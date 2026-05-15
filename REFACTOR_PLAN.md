# SoupX → Python + Rust 重构方案

> **原项目**: [constantAmateur/SoupX](https://github.com/constantAmateur/SoupX) (R package, GigaScience 2020)
> **目标**: 将 core 算法用 Rust 重写，Python API 集成到 scanpy/AnnData 生态，实现 8–15× 端到端加速和 90%+ 内存节省

---

## 目录

1. [原项目分析](#1-原项目分析)
2. [技术决策：为什么是 Rust 而非 Go](#2-技术决策为什么是-rust-而非-go)
3. [分层架构设计](#3-分层架构设计)
4. [Rust Core 模块设计](#4-rust-core-模块设计)
5. [Python API 设计](#5-python-api-设计)
6. [性能预估与量化分析](#6-性能预估与量化分析)
7. [内存分析与 OOM 风险评估](#7-内存分析与-oom-风险评估)
8. [开发路线图](#8-开发路线图)
9. [风险与缓解](#9-风险与缓解)

---

## 1. 原项目分析

### 1.1 是什么

SoupX 是 R 包，用于去除液滴法单细胞 RNA-seq（10X Genomics）中的 **ambient mRNA 污染**（"soup"）。液滴法的背景污染普遍占总表达量的 2–50%，会产生批次效应影响下游分析。

### 1.2 算法管线（4 步）

```
load10X() / SoupChannel()
    ↓ 创建 SoupChannel 对象，自动 estimateSoup()
setClusters()
    ↓ 提供聚类信息（必需，否则效果极差）
autoEstCont()
    ↓ 自动估计全局污染比例 ρ
adjustCounts()
    ↓ 三种校正方法，输出净化后的计数矩阵
```

### 1.3 核心数据结构

```r
SoupChannel:
├── $tod          # Table of Droplets（所有液滴矩阵，内存大头 85%+）
├── $toc          # Table of Counts（仅细胞液滴）
├── $soupProfile  # data.frame: est(基因在soup中的比例), counts
├── $metaData     # data.frame: nUMIs, clusters, rho, DR坐标
├── $fit          # 拟合结果（GLM 或后验密度）
└── $DR           # 降维坐标
```

### 1.4 六大核心算法

#### (a) `estimateSoup` — 估计 soup profile

**目标**: 获取 ambient RNA 中每个基因的比例。

**方法**: 取 UMI 总数在 [0, 100] 的空液滴，假设其全部为背景 RNA：
```
soupProfile$est[gene_i] = sum(空液滴中 gene_i 的 counts)
                         / sum(空液滴中所有基因的 total counts)
```

#### (b) `quickMarkers` — 快速标志基因发现

**目标**: 为每个簇找特异性基因（tf-idf 方法）。

**算法**:
1. 二值化稀疏矩阵（counts > 0.9 → 表达）
2. 对每个簇 × 每个基因计算 **tf-idf score**：
   ```
   TF  = 基因在簇中的表达频率
   IDF = log(总细胞数 / 表达该基因的细胞数)
   score = TF × IDF
   ```
3. 超几何检验 + BH 校正，取每个簇的 top N

**复杂度**: O(n_clusters × n_genes)，天然适合并行。

#### (c) `estimateNonExpressingCells` — 找出不表达某基因集的细胞

**目标**: 给定一组基因（如 HB 基因），找出确定不表达它们的细胞，用于 ρ 估计。

**算法**:
1. 对每个基因集 × 每个细胞，计算预期背景计数：
   ```
   exp = total_UMSs * maxContam * sum(soupProfile[geneSet])
   ```
2. Poisson 检验：`p = P(X ≥ obs-1 | λ = exp)`
3. BH 校正
4. 簇级决策：如果簇中任何细胞的 p 值显著 → 整个簇排除

#### (d) `calculateContaminationFraction` — 手动模式估计 ρ

**前提**: 用户提供在部分细胞中确定不表达的基因集。

**算法**:
1. 用 `estimateNonExpressingCells` 过滤可用的细胞
2. 对过滤后的细胞，拟合 **Poisson GLM**（log-link + offset）：
   ```
   counts ~ Poisson(λ)
   log(λ) = log(expected_soup_counts) + log(ρ)
   ```
3. `ρ = exp(intercept)` 即全局污染比例

#### (e) `autoEstCont` — 自动模式估计 ρ（核心）

**核心洞察**: 多个不同簇独立给出 n 个 ρ 估计值 → 真值附近密度聚集。

**算法步骤**:
1. `quickMarkers` 找标志基因，过滤 tfidf > tfidfMin 且 soup 表达高
2. 对每个标志基因 × 每个簇，`estimateNonExpressingCells` 过滤
3. 对每个可用估计建立 **Gamma 后验分布**：
   ```
   先验：Gamma(k, θ)    (mode=priorRho=0.05, sd=priorRhoStdDev=0.10)
   后验：Gamma(k + obsCnt, scale = θ / (1 + θ × expCnt))
   ```
4. 将所有后验密度聚合（按基因加权），取 contaminationRange 内的 MAP 点

**瓶颈**: `sapply(rhoProbes, ...)` 在 1000 个探针点上逐点算密度，完全串行。

#### (f) `adjustCounts` — 计数校正（最核心、最复杂）

**三种方法**:

| 方法 | 原理 | 速度 | 精确度 |
|------|------|:---:|:---:|
| **subtraction** | 从每个基因迭代减去期望背景计数 | 中 | 高 |
| **soupOnly** | Poisson 检验鉴定纯污染基因整基因删除 | 快 | 中 |
| **multinomial** | 显式最大化多项似然，贪心交换分子 | 最慢 | 最高 |

**subtraction 细节**:
```
对每个细胞:
  expSoupCnts = nUMIs × ρ
  alloc(expSoupCnts, column, soupProfile)
    ↓ 按权重分配到各基因，受限于 observed counts
  迭代直到剩余 < tol
```

**`alloc` 子算法**（关键工具函数）:
```
输入: tgt(目标总量), bucketLims(各基因上限), ws(权重)
1. 按 bucketLims/ws 排序（确定"填满顺序"）
2. 依次填满每个桶直到达到上限
3. 剩余量按重新归一化的权重分配
输出: 每个基因应减去的背景计数（≤ 观察值）
```

**soupOnly 细节**:
1. 对每个非零 entry 计算 Poisson p 值
2. 按 cell 内 p 值排序，累积求和
3. Fisher 联合 p 值（χ²(4)），pCut 阈值过滤
4. *关键瓶颈*：全局排序 `order(-(out@j+1), p)`

**multinomial 细节**:
```
对每个细胞:
  初始化 fit = 减法结果
  while True:
    找最大 delInc（增加似然的基因）和 delDec（减少似然的基因）
    贪心交换
    直到 delInc + delDec ≤ 0
```

---

## 2. 技术决策：为什么是 Rust 而非 Go

| 维度 | Rust | Go |
|------|------|-----|
| **Python 互操作** | PyO3 + maturin，直接暴露 numpy 内存，零拷贝 | cgo + 手动 FFI，每个矩阵传递都要复制 |
| **稀疏矩阵** | `sprs` 直接操作 CSR/CSC，对标 scipy 内部格式 | 无成熟库，需手写 |
| **数值计算** | `ndarray` + `sprs` + `statrs`，生态丰富 | 手写或 C 绑定 |
| **并行化** | `rayon` 无锁并行迭代器，极简 | goroutine 但数据共享复杂 |
| **包分发** | `maturin build --release` → 单 wheel，pip install 即可 | 需要 C 编译链 + Go 工具链 |
| **内存控制** | 零成本抽象 + 确定性析构 | GC 暂停 + 内存不可预测 |

**结论**: Rust 是唯一合理的选择。Go 在数值计算 + Python 互操作的场景下没有竞争力。

---

## 3. 分层架构设计

```
┌─────────────────────────────────────────────────────┐
│  用户 API（纯 Python）                                │
│                                                      │
│  soupx.load_10x(adata, cellranger_dir)               │
│  soupx.set_clusters(adata, key="leiden")             │
│  soupx.auto_est_cont(adata)                          │
│  soupx.adjust_counts(adata, method="subtraction")    │
│  soupx.decontaminate(adata)          # 一步到位       │
│  soupx.plot_marker_distribution(adata, gene_list)    │
│  soupx.plot_marker_map(adata, "HB")                  │
└──────────────┬──────────────────────────────────────┘
               │
┌──────────────▼──────────────────────────────────────┐
│  Python Glue Layer（soupx/_bridge.py）               │
│                                                      │
│  - AnnData → numpy CSR/CSC arrays                    │
│  - 调用 Rust 核心（零拷贝传递 numpy buffer）           │
│  - Cluster label → vector mapping                    │
│  - 10X CellRanger 目录解析                            │
│  - 可视化封装（matplotlib/plotly）                     │
└──────────────┬──────────────────────────────────────┘
               │  PyO3 (zero-copy numpy)
┌──────────────▼──────────────────────────────────────┐
│  Rust Core Engine（Python 扩展模块）                   │
│                                                      │
│  ┌──────────────────────────────────────────────┐   │
│  │           Public API (PyO3 #[pyfunction])     │   │
│  │  quick_markers, estimate_contam,              │   │
│  │  auto_est_cont, adjust_counts, expand_to_cells│   │
│  └──────────────────┬───────────────────────────┘   │
│                     │                                 │
│  ┌──────────────────▼───────────────────────────┐   │
│  │              Engine Layer                     │   │
│  │  ┌───────────────┐  ┌────────────────────┐   │   │
│  │  │ Marker Engine │  │ Contam Estimator   │   │   │
│  │  │ (tf-idf +     │  │ (Poisson GLM +     │   │   │
│  │  │  hypergeomet.)│  │  Gamma posterior)  │   │   │
│  │  └───────────────┘  └────────────────────┘   │   │
│  │  ┌──────────────────────────────────────┐    │   │
│  │  │    Count Adjustment Engine            │    │   │
│  │  │    subtraction / soupOnly / multinom. │    │   │
│  │  └──────────────────────────────────────┘    │   │
│  └──────────────────┬───────────────────────────┘   │
│                     │                                 │
│  ┌──────────────────▼───────────────────────────┐   │
│  │             Core Utilities                    │   │
│  │  ┌─────────┐ ┌──────────┐ ┌──────────────┐  │   │
│  │  │ Sparse  │ │  Stats   │ │   Allocator  │  │   │
│  │  │  Matrix │ │ (Poisson,│ │  (alloc +    │  │   │
│  │  │  Ops    │ │  Gamma,  │ │   expand)    │  │   │
│  │  │         │ │  Hyperg.)│ │              │  │   │
│  │  └─────────┘ └──────────┘ └──────────────┘  │   │
│  └──────────────────────────────────────────────┘   │
└─────────────────────────────────────────────────────┘
```

---

## 4. Rust Core 模块设计

### 4.1 Crate 结构

```
soupx-core/
├── Cargo.toml
├── src/
│   ├── lib.rs              # PyO3 入口，暴露 #[pyfunction] 到 Python
│   ├── sparse.rs           # CSR/CSC 零拷贝适配层
│   ├── stats.rs            # Poisson, Gamma, Hypergeometric 分布函数
│   ├── soup_profile.rs     # estimate_soup（可保留 Python，极简单）
│   ├── markers.rs          # quick_markers: tf-idf + hypergeometric + BH
│   ├── non_expressing.rs   # estimate_non_expressing_cells
│   ├── contamination.rs    # calculate_contamination + auto_est_cont
│   ├── adjustment/
│   │   ├── mod.rs          # adjust_counts 总入口
│   │   ├── subtraction.rs
│   │   ├── soup_only.rs    # soupOnly 方法（CSC 天然列排序优化）
│   │   └── multinomial.rs  # multinomial 方法（BinaryHeap 优化）
│   ├── alloc.rs            # alloc 分配算法（分级策略）
│   └── expand.rs           # expand_clusters 簇→单细胞展开
```

### 4.2 关键依赖

```toml
[dependencies]
pyo3 = { version = "0.22", features = ["extension-module"] }
numpy = "0.22"           # Python numpy 互操作
sprs = "0.11"            # 高性能稀疏矩阵 (CSR/CSC)
ndarray = "0.16"         # 密集多维数组
statrs = "0.16"          # 统计分布 (也可能自写，取决于 statrs 精度)
rayon = "1.10"           # 数据并行
anyhow = "1"             # 错误处理
thiserror = "1"          # 结构化错误类型
```

### 4.3 关键技术点

#### (a) 数据桥接与 GIL 安全策略

**分级处理：tod 与 toc 区别对待。**

`soup_profile` 估计只需要空液滴子集（UMI ∈ [0,100]）的简单求和，Python/numpy 完全可胜任，无需传给 Rust。而 `toc`（细胞计数矩阵）是所有重计算的输入，且体积可控。

| 数据 | 大小（典型） | 处理方式 | 理由 |
|------|:---:|------|------|
| `tod`（全液滴） | ~2.5 GB | Python 层估计 soup profile，**释放后不再持有** | 只在 load_10x 中用一次，无需 clone |
| `toc`（细胞计数） | ~0.15 GB | 入口处 clone 到 Rust 自有 `spars::CsMat`，之后 Rust 全权持有 | 数据量小，clone 代价 <0.1s，换 GIL 安全 |

**GIL 释放策略**：全部重计算函数（markers、contam、adjust）在 Rust 入口处拿到自有 CSR 后立即调用 `py.allow_threads()` 释放 GIL，允许 Python 端并发。绝不跨 `allow_threads` 边界持有 numpy buffer 指针。

**为什么不用零拷贝视图**：
- `CsMatView` 指向 Python numpy buffer，生命周期依赖 Python 端
- 一旦 Rust 释放 GIL 并行计算，Python GC 可能回收底层数组 → use-after-free
- 约束"永不释放 GIL"则锁死 Python 端，浪费多线程能力
- `toc` clone 代价（~0.15 GB, ~0.1s）远小于收益（完整 GIL 释放 + 零 UB 风险）

```rust
use numpy::{PyArray1, PyReadonlyArray1};
use sprs::CsMat;

/// 安全路径（推荐默认）：从 numpy CSR 复制构造自有 sprs::CsMat
/// clone 代价：50K cells × 25K genes × 5% sparsity ≈ 0.15 GB, ~0.1s
pub fn csr_from_numpy_owned(
    data: &PyArray1<f64>,
    indices: &PyArray1<i64>,
    indptr: &PyArray1<i64>,
    shape: (usize, usize),
) -> CsMat<f64> {
    // 读取 numpy buffer 并复制到自有 CsMat（在 GIL 持有期间完成）
}

/// 零拷贝路径（仅限短生命周期，不释放 GIL，文档标注 unsafe）
/// 适用场景：单次简单操作，不跨越 py.allow_threads
pub unsafe fn csr_from_numpy_view<'py>(
    data: &'py PyArray1<f64>,
    indices: &'py PyArray1<i64>,
    indptr: &'py PyArray1<i64>,
    shape: (usize, usize),
) -> sprs::CsMatView<'py, f64> {
    // 仅在持有 GIL 时使用，绝不与 allow_threads 共存
}
```

#### (b) SoupOnly 的 CSC 优化（最大单项优化，15–25×）

```rust
/// R 原版：
///   order(-(out@j+1), p)     # 全局排序，O(N log N)
///   split(s, out@j[o]+1)     # 按 cell 分割排序结果
///
/// Rust 版（CSC 格式）：
///   对每一列（= 一个细胞）：
///     该列的 entries 天然连续存储，无需排序
///     用 nth_element 找 pCut 分界点
///     只输出保留的 entries
///
/// 时间从 O(N log N) 降为 O(N + n_cells × log k)
```

#### (c) Alloc 的分级优化

```rust
/// R 原版：
///   order(bucketLims/ws)  # 全排序，O(k log k)
///
/// Rust 版（分级策略）：
///   1. 早期短路：所有桶已满或全为零 → O(k) 检测后直接返回
///   2. 少量分配（tgt 很小）：partial_sort 前 m 个桶 → O(k log m)
///   3. 一般情况：全排序 → O(k log k)
///
/// 注意：quickselect 只能定位第 m 统计量，不能保留填充顺序；
/// alloc 需要按 bucketLims/ws 升序依次填满，必须用 partial_sort 或全排序。
```

#### (d) Multinomial 的 BinaryHeap 优化

```rust
/// R 原版：
///   while True:
///     delInc = log(ps[increasable]) - log(fit[increasable]+1)   # O(k) 扫描
///     wInc = which(increasable)[which.max(delInc)]
///     delDec = ...
///     wDec = which(decreasable)[which.max(delDec)]
///
/// Rust 版：
///   维护两个 BinaryHeap: incHeap (max-heap of delInc), decHeap (max-heap of delDec)
///   每次迭代: pop 堆顶 → O(log k)，更新后 push 回 → O(log k)
///   总复杂度: O(iter × log k) vs 原版 O(iter × k)
```

### 4.4 并行化策略

| 函数 | 并行维度 | 工具 | GIL 策略 | 预期线程利用率 |
|------|---------|------|:---:|:---:|
| `quickMarkers` | 按簇（n_clusters 路并行） | `rayon::par_iter` | 释放 GIL（自有数据） | 高（簇数 >> CPU 核数时） |
| `estimateNonExpressing` | 按基因集（n_geneSets 路并行） | `rayon::par_iter` | 释放 GIL（自有数据） | 中（通常 5–20 个基因集） |
| `autoEstCont` 密度计算 | 按 rhoProbe 点（1000 路并行） | `rayon::par_iter` | 释放 GIL（自有数据） | 极高（1000 个独立计算） |
| `adjustCounts` | 按细胞/簇（n_cells 路并行） | `rayon::par_iter`（按 indptr 分块） | 释放 GIL（自有数据） | 极高（细胞数 >> CPU 核数） |
| `expandClusters` | 按簇（n_clusters 路并行） | `rayon::par_iter` | 释放 GIL（自有数据） | 高 |
| `alloc` 内部 | 不并行（每个 call 太小，~几十个基因） | 串行 | — | — |

---

## 5. Python API 设计

### 5.1 顶层 API

```python
import soupx
import scanpy as sc

# === 完整工作流（推荐） ===
adata = sc.read_10x_mtx("path/to/cellranger/outs/")

# 1. 加载空液滴数据（从 10X raw_gene_bc_matrices 读取）
soupx.load_10x(
    adata,
    cellranger_dir="path/to/cellranger/outs/",
    soup_range=(0, 100),        # UMI 范围，同 R 版的 soupRange
    keep_droplets=False,        # 是否保留完整液滴矩阵
)

# 2. 正常聚类（scanpy 标准流程）
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata)
soupx.set_clusters(adata, key="leiden")

# 3. 自动估计污染比例
soupx.auto_est_cont(adata)

# 4. 校正计数（输出到 adata.layers["soupx_corrected"]）
soupx.adjust_counts(adata, method="subtraction", round_to_int=True)

# === 一步到位 ===
adata = soupx.decontaminate(
    adata,
    clusters="leiden",
    method="subtraction",
    round_to_int=True,
)

# === 高级用法：手动指定污染基因集 ===
gene_list = {"HB": ["HBB", "HBA2", "HBD"],
             "IG": ["IGHG1", "IGHG2", "IGKC", "IGLC2"]}
rho_est = soupx.estimate_contamination(adata, gene_list)
print(f"Estimated contamination: {rho_est:.2%}")

# === 可视化 ===
soupx.plot_marker_distribution(adata, gene_list)
soupx.plot_marker_map(adata, "HB")
soupx.plot_change_map(adata, "HB")  # 需要先运行 adjust_counts

# === 检查结果 ===
print(adata.obs["soupx_rho"].describe())
print(adata.uns["soupx_fit"])          # 拟合详情
print(adata.varm["soup_profile"].head())  # soup 表达谱
```

### 5.2 AnnData 结构约定

| AnnData 位置 | 内容 | 类型 | 说明 |
|-------------|------|------|------|
| `adata.X` | 细胞计数矩阵 | `scipy.sparse.csr_matrix` | 标准 AnnData 格式 |
| `adata.raw.X` | 所有液滴矩阵 | `scipy.sparse.csr_matrix` | 由 `load_10X` 设置 |
| `adata.obs["soupx_rho"]` | 每个细胞的 ρ | `float` | 由 `auto_est_cont` / `estimate_contamination` 设置 |
| `adata.obs["clusters"]` | 聚类标签 | `str`（或自定义 key） | 由 `set_clusters` 指定 key |
| `adata.uns["soupx_fit"]` | 拟合详情 | `dict` | 包括 prior, posterior, rhoEst, markersUsed 等 |
| `adata.varm["soup_profile"]` | 基因 × [est, counts] | `DataFrame` | soup 表达谱 |
| `adata.layers["soupx_corrected"]` | 净化后的计数矩阵 | `scipy.sparse.csr_matrix` | `adjust_counts` 的输出 |
| `adata.obsm["X_umap"]` 或 `X_tsne` | 降维坐标 | `DataFrame` | 可视化需要 |

### 5.3 函数签名

```python
def load_10x(
    adata: AnnData,
    cellranger_dir: str | Path,
    soup_range: tuple[int, int] = (0, 100),
    keep_droplets: bool = False,
) -> None:
    """从 10X CellRanger 输出目录加载所有液滴数据。"""

def set_clusters(
    adata: AnnData,
    key: str = "leiden",
) -> None:
    """设置用于背景估计的聚类标签。"""

def auto_est_cont(
    adata: AnnData,
    *,
    tfidf_min: float = 1.0,
    soup_quantile: float = 0.90,
    max_markers: int = 100,
    contamination_range: tuple[float, float] = (0.01, 0.8),
    prior_rho: float = 0.05,
    prior_rho_stddev: float = 0.10,
    force_accept: bool = False,
    verbose: bool = True,
) -> None:
    """自动估计全局污染比例 ρ。（对应 R 版 autoEstCont）"""

def estimate_contamination(
    adata: AnnData,
    gene_list: dict[str, list[str]],
    *,
    maximum_contamination: float = 1.0,
    fdr: float = 0.05,
    force_accept: bool = False,
) -> float:
    """使用手动指定的基因集估计 ρ。（对应 calculateContaminationFraction）"""

def adjust_counts(
    adata: AnnData,
    *,
    method: Literal["subtraction", "soup_only", "multinomial"] = "subtraction",
    round_to_int: bool = False,
    p_cut: float = 0.01,
    verbose: bool = True,
) -> AnnData:
    """去除背景污染，将净化后的计数矩阵写入 adata.layers["soupx_corrected"]。

    原地修改 adata 并返回它以支持链式调用（与 scanpy 惯例一致）。
    """

def decontaminate(
    adata: AnnData,
    *,
    clusters: str = "leiden",
    method: Literal["subtraction", "soup_only", "multinomial"] = "subtraction",
    round_to_int: bool = True,
    **auto_est_kwargs,
) -> AnnData:
    """便捷接口：设置聚类 → 自动估计 ρ → 校正计数。

    前提：已通过 load_10x() 加载所有液滴数据，且已完成 scanpy 聚类（如 sc.tl.leiden）。
    原地修改 adata 并返回以支持链式调用。
    """

def quick_markers(
    adata: AnnData,
    *,
    n: int = 10,
    fdr: float = 0.01,
    express_cut: float = 0.9,
) -> pd.DataFrame:
    """快速标志基因发现（tf-idf + 超几何检验）。"""
```

---

## 6. 性能预估与量化分析

### 6.1 方法论说明

**重要**: R 原版中统计分布函数（`ppois`, `phyper`, `dgamma`, `qgamma`）已调用 C/Fortran 编译代码，这部分 Rust 无法超越。Rust 的加速来自其他地方。

### 6.2 分模块预估

典型数据集：50K 细胞、2M 空液滴、25K 基因、稀疏度 ~5%

| 模块 | R 瓶颈 | Rust 优化点 | 预估加速 | 理由 |
|------|--------|-----------|:---:|------|
| **`quickMarkers`** | `dgTMatrix` 遍历 + `lapply`/`sapply` + 串行 | rayon 并行 + CSR 列切片 + 避免 data.frame 构造 | **8–15×** | 统计计算同等，赢在并行和结构 |
| **`estimateNonExpressingCells`** | `do.call(rbind, lapply(...))` 重复矩阵分配 + 串行 | rayon 双维并行（基因集×簇）+ 原地计算 | **6–10×** | 主要是并行 |
| **`autoEstCont`** 密度计算 | `sapply(rhoProbes, ...)` 1000 个点全串行 | rayon 在 rhoProbes 维并行 + SIMD dgamma | **10–20×** | Embarrassingly parallel，最适合加速 |
| **`adjustCounts` subtraction** | 逐细胞 `lapply` + 逐基因 `alloc()` → `order()` O(k log k) | rayon 逐列并行 + quickselect O(k) | **8–15×** | 算法降复杂度 + 并行 |
| **`adjustCounts` soupOnly** | **全局排序** `order(-(out@j+1), p)` + split | CSC 天然列分组 → **免排序** + nth_element | **15–25×** | **最大单项优化** |
| **`adjustCounts` multinomial** | while 循环中每轮 O(k) 线性扫描 | BinaryHeap 优先队列 O(log k) + 并行 | **20–40×** | 算法复杂度降级 |

### 6.3 端到端预估

| 数据集 | R 原版 | Rust 预估 | 加速比 |
|--------|:---:|:---:|:---:|
| 10K cells | ~45s | ~5s | 9× |
| 50K cells | ~240s (4 min) | ~26s | 9× |
| 200K cells | ~900s (15 min) | ~80s (1.3 min) | 11× |
| 1M cells | ~1500s (25 min) | ~120s (2 min) | 12× |

**保守估计：端到端 8–15 倍加速，数据越大加速越明显。**

---

## 7. 内存分析与 OOM 风险评估

### 7.1 原版 R 内存剖面（50K cells, 2M droplets, 25K genes）

```
对象                        大小         占比
────────────────────────────────────────────
sc$tod (dgCMatrix)         ~2.5 GB       85%
sc$toc (dgCMatrix)         ~0.15 GB       5%
sc$metaData               ~0.01 GB      <1%
sc$soupProfile            ~0.01 GB      <1%
────────────────────────────────────────────
稳态基线                   ~2.7 GB
────────────────────────────────────────────
adjustCounts 中间态:
  tmp$toc (簇聚合 CSR)      ~0.01 GB
  expandClusters 合并        ~0.3 GB
  Tsparse 操作缓冲           ~0.5 GB
────────────────────────────────────────────
峰值                       ~3.5 GB
```

### 7.2 Rust 版内存剖面（相同数据）

```
对象                        大小         说明
────────────────────────────────────────────────
tod CSR view               ~0 GB        零拷贝读取 numpy buffer
toc CSR view               ~0 GB        零拷贝
metaData Vecs              ~0.01 GB
soupProfile Vec            ~0.005 GB
────────────────────────────────────────────────
稳态基线                   ~0.015 GB    比 R 省 ~99%
────────────────────────────────────────────────
adjustCounts 中间态:
  簇聚合 CSR                ~0.01 GB
  并行列缓冲区              ~0.01 GB × N_threads
  expandClusters 结果        ~0.15 GB
────────────────────────────────────────────────
峰值                       ~0.25 GB     比 R 省 ~93%
```

### 7.3 关键结论

| 维度 | R 原版 | Rust 版 |
|------|:---:|:---:|
| 自身内存占用（不含数据） | ~1.0 GB 峰值 | ~0.25 GB 峰值 |
| 内存节省 | — | **90–95%** |
| 同台机器可处理数据上限 | 100K cells (16GB) | **500K–800K cells** (16GB) |

**OOM 硬限制**: `tod` 矩阵本身（2M 液滴 × 25K 基因 × 5% 稀疏 × 12 bytes/entry = ~30 GB 密集等效，但稀疏压缩后 ~2.5 GB）。如果 `tod` 本身 > 机器 RAM 的 60%，无论什么语言都会 OOM。Rust 的优势是**消除额外 90% 的自身开销**，把限制推迟到硬数据边界。

### 7.4 超大数据的应对方案

| 方案 | 代价 | 适用范围 |
|------|------|---------|
| 采样空液滴（取 1/10 估计 soup profile） | 精度轻微下降 | 一般场景可接受 |
| mmap 读取 `tod`（内存映射文件） | IO 变慢 | >5M 液滴 |
| 按通道分批处理 | 无法跨通道共享信息 | 多通道数据 |
| 基因预过滤（去除在所有空液滴中计数 <10 的基因） | 丢失少量低表达基因 | 可配置 |

---

## 8. 开发路线图

### 里程碑

| 周次 | 里程碑 | 可交付物 |
|:---:|--------|---------|
| 1 | 项目脚手架 | Cargo.toml + PyO3 骨架 + `maturin build` 通过 + CI |
| 2 | 稀疏矩阵桥接 + 统计函数库 | CSR/CSC 零拷贝 + Poisson/Gamma/Hypergeo + BH → 单元测试对账 R 输出 |
| 3 | `quickMarkers` Rust 实现 | Python 调 Rust，对账 R 版 quickMarkers 结果 |
| 4 | `estimateNonExpressingCells` + `calculateContaminationFraction` | 手动 ρ 估计可用 |
| 5–6 | `autoEstCont` 完整实现 | 后验密度聚合 + Gamma 先验 |
| 7 | `adjustCounts` 三种方法 + `expandClusters` | 核心校正完成 |
| 8 | Python 高层 API + AnnData 集成 + `load10X` | 可用 pip install 安装使用 |
| 9 | 可视化 + 文档 + 示例 notebook | `plotMarkerDistribution/Map/ChangeMap` 对标 R 版 |
| 10 | 系统对账 + 性能 benchmark + PyPI 发布 | 与原版在标准数据集上逐层对账（ρ, 校正矩阵） |

### 对账策略

与 R 原版逐层对账，确保结果一致：

```
Level 1: soup profile      → expect identical（简单求和，应完全相同）
Level 2: quickMarkers      → expect tfidf identical, qval within 1e-8
Level 3: ρ estimate        → expect within 1e-6
Level 4: adjusted counts   → expect >99.9% entries identical（浮点路径差异）
```

---

## 9. 风险与缓解

| 风险 | 概率 | 影响 | 缓解 |
|------|:---:|:---:|------|
| statrs 统计函数精度不够（与 R 的 C 实现相比） | 中 | 中 | 对账 → 必要时自写纯 Rust 统计函数（对标 Rmath） |
| sprs 与 scipy CSR 内部格式不兼容 | 低 | 高 | 零拷贝适配层是重中之重，首周验证 |
| `autoEstCont` 密度计算的数值稳定性 | 中 | 中 | 使用 log-dgamma 避免下溢；与 R 对账每种场景 |
| rayon 并行导致结果不可重现（浮点累加顺序不同） | 中 | 低 | 提供 `n_threads=1` 模式用于测试；文档说明 |
| pyo3/maturin 版本升级破坏兼容性 | 低 | 中 | 锁定依赖版本 + CI 覆盖率矩阵 |
| 10X CellRanger 输出格式变化 | 低 | 低 | scanpy 社区已处理，复用 `sc.read_10x_mtx` |
| 用户不接受 AnnData 约定 | 低 | 低 | 明确文档 + 提供 `soupx.SoupChannel` 类（兼容原 R 用户心智模型） |
