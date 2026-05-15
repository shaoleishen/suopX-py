# REFACTOR_PLAN.md 审核意见

> 审核日期：2026-05-14
> 审核人：bioshen
> 总体评价：**计划方向正确，细节需打磨。核心修正 3 项，优化建议 14 项。**

---

## 总体评价

**优点**：Rust 选型论证充分，分层架构清晰合理，内存/性能分析有量化依据，风险识别诚实。

**核心问题**：有几个技术细节存在偏差，以及一些工程实践的缺失。以下按严重程度排列。

---

## 需要修正的问题

### 1. Alloc 的 Quickselect 优化逻辑存疑（原文 4.3c）

> "找填满分界点只需要第 m 顺序统计量，O(k) … 对大多数 trivial 基因（全 soup 或零 soup），直接跳过"

`alloc` 的核心逻辑是**按 `bucketLims/ws` 升序依次填满每个桶**，并非只找一个分界点：

```
1. 按 bucketLims/ws 排序（确定"填满顺序"）
2. 依次填满每个桶直到达到上限
3. 剩余量按重新归一化的权重分配
```

你需要的是完整的有序序列，而非仅第 m 个顺序统计量。

- quickselect（`nth_element`）只定位第 m 小的元素，**不保留前 m 个的排序关系**，无法直接替代
- 如果你的意图是：对 trivial case 早期短路（全零桶直接返回零向量，全满桶直接返回 `bucketLims`），这与 quickselect 是两个独立概念，不应混淆

**建议**：

| 场景 | 策略 | 复杂度 |
|------|------|:---:|
| 所有桶已满 / 所有桶为零 | 早期短路，O(k) 检测 | O(k) |
| `tgt` 很小，只有少量桶被填充 | `partial_sort` 前 m 个桶 | O(k log m) |
| 一般情况 | 全排序 | O(k log k) |

文档中应修正 "O(k) quickselect 版" 的说法，改为 "O(k) 早期短路 + O(k log k) 全排序"。

---

### 2. Decontaminate 签名与文档矛盾（原文 5.1 + 5.3）

```python
def decontaminate(adata, *, clusters="leiden", method="subtraction", **auto_est_kwargs):
    """一步完成：加载 → 聚类 → 估计 ρ → 校正"""
```

说"加载 → 聚类"，但参数里**没有 `cellranger_dir`，也没有聚类逻辑**。实际流程是：

```python
# 用户在外部完成这些步骤：
soupx.load_10x(adata, cellranger_dir, ...)  # 加载
sc.tl.leiden(adata)                          # 聚类
# 然后 decontaminate 只做：
soupx.set_clusters(adata, key=clusters)     # 设置聚类
soupx.auto_est_cont(adata, **auto_est_kwargs)  # 估计 ρ
soupx.adjust_counts(adata, method=method)    # 校正
```

**建议**：文档字符串改为：

> "便捷接口：设置聚类 → 自动估计 ρ → 校正计数。**前提：已通过 load_10x() 加载数据，已完成 scanpy 聚类。**"

---

### 3. Zero-copy 桥接的 GIL 安全风险（原文 4.3a）

原文描述了零拷贝构造 `CsMatView`，但**完全未涉及 GIL（全局解释器锁）管理**。这是整个项目最脆弱的环节。

**风险场景**：

```
1. Rust 通过 unsafe 构造 CsMatView<'a>，指向 Python numpy 的 data/indices/indptr buffer
2. Rust 调用 py.allow_threads(|| { rayon 并行计算 }) 释放 GIL
3. 此时如果 numpy 数组在 Python 端失去引用（GC），Rust 持有的指针悬垂
   → use-after-free → Undefined Behavior
```

**建议在架构中明确补充**：

#### (a) 两种视图策略

```rust
// 策略 A：零拷贝视图（不安全，需约束调用者）
// 适用场景：短生命周期调用，不释放 GIL，不跨越 py.allow_threads
pub unsafe fn csr_view_from_numpy<'py>(...) -> CsMatView<'py, f64> {
    // 仅在持有 GIL 时构造和使用
}

// 策略 B：安全克隆（推荐作为默认）
// 适用场景：释放 GIL 进行长时间计算
pub fn csr_owned_from_numpy(py: Python, csr: &PyAny) -> CsMat<f64> {
    // 从 numpy 读取并复制到自有 sprs::CsMat
    // 复制代价 ~0.1-0.2 GB（50K cells × 25K genes × 5% sparsity）
    // 相对全年节省 90%+ 内存的优势，这点复制代价可接受
}
```

#### (b) GIL 释放策略表

| 函数 | 使用的数据 | 是否释放 GIL | 理由 |
|------|-----------|:---:|------|
| `quick_markers` | 自有 CSR（已 clone） | ✅ 释放 | 纯 Rust 数据 |
| `estimate_non_expressing` | 自有 CSR | ✅ 释放 | 纯 Rust 数据 |
| `auto_est_cont` 密度扫描 | 自有 Vec | ✅ 释放 | 纯 Rust 数据 |
| `adjust_counts` subtraction | 自有 CSR | ✅ 释放 | 纯 Rust 数据 |
| `adjust_counts` soupOnly | 自有 CSR | ✅ 释放 | 纯 Rust 数据 |
| `adjust_counts` multinomial | 自有 CSR | ✅ 释放 | 纯 Rust 数据 |

**推荐方案**：在 `lib.rs` 入口处统一 clone 一次 CSR（O(N) 复制），之后全部释放 GIL 计算。**牺牲 ~0.15 GB 内存换来完整的 GIL 释放和零 UB 风险。** 这比原文中 `unsafe` 零拷贝视图的方案更安全、更简单。

#### (c) 如果坚持零拷贝

如果必须零拷贝，至少需要：
1. 用 `Python::with_gil(|py| { ... })` 包裹所有持有视图的计算，绝不释放 GIL
2. 这意味着**整个计算期间锁死 GIL**，Python 端无法并发
3. 文档显式标注 `#[pyfunction]` 持有 GIL 的行为

---

### 4. `adjust_counts` 返回值语义不明确（原文 5.3）

```python
def adjust_counts(adata, ...) -> AnnData:
    """去除背景污染，返回净化后的计数矩阵。"""
```

文档注释说"返回净化后的计数矩阵"，但签名返回 `AnnData`。实际语义有两种可能：

| 语义 | 行为 | 优缺点 |
|------|------|--------|
| **In-place 修改 + 返回自身** | 写入 `adata.layers["soupx_corrected"]`，返回同一个对象 | 可链式调用，但容易误以为返回了新对象 |
| **返回新 AnnData** | 复制 adata，修改副本的 `.X`，返回副本 | 语义清晰但内存翻倍 |

**建议**：明确选择一种，并在文档中标注。推荐 in-place + 返回自身（与 scanpy 惯例一致，如 `sc.pp.normalize_total`）：

```python
def adjust_counts(
    adata: AnnData,
    *,
    method: Literal["subtraction", "soup_only", "multinomial"] = "subtraction",
    round_to_int: bool = False,
    p_cut: float = 0.01,
    verbose: bool = True,
) -> AnnData:
    """去除背景污染，将净化后的计数矩阵写入 adata.layers["soupx_corrected"]。
    
    原地修改 adata 并返回它以支持链式调用。
    """
```

---

## 设计与工程优化建议

### 5. `load_10x` 默认行为浪费 I/O 和内存（原文 5.1）

```python
soupx.load_10x(adata, cellranger_dir, keep_droplets=False)  # 默认
```

`keep_droplets=False` 是默认，意味着加载全部 2M+ 液滴矩阵，估计完 soup profile 后立刻丢弃。这浪费了：

- 读取 2M 液滴的 I/O 时间（~数秒到数十秒）
- 在内存中短暂持有 ~2.5 GB 的 `tod` 矩阵

**建议增加 `soup_profile_only` 参数**：

```python
def load_10x(
    adata: AnnData,
    cellranger_dir: str | Path,
    soup_range: tuple[int, int] = (0, 100),
    keep_droplets: bool = False,
    soup_profile_only: bool = True,   # 新增：只读空液滴，不加载全矩阵
) -> None:
```

- `soup_profile_only=True`：只读取 UMI ∈ soup_range 的空液滴，直接估计 soup profile，不构造 `tod`
- `keep_droplets=True`：全量加载并保留（高级用户需要）

这可将 I/O 从全矩阵降为空液滴子集，通常节省 80%+。

---

### 6. 缺少非 10X 数据的通用 API 入口

原文只在风险表中提到"提供 `SoupChannel` 类"作为缓解措施（第 9 节），但未在 API 设计中体现。这是社区采用率的关键——非 10X 平台（inDrops, Drop-seq, Smart-seq2 等）的用户也需要这个工具。

**建议在 Python API 中增加**：

```python
class SoupChannel:
    """SoupX 核心对象的 Python 版本，兼容非 10X 数据源。"""
    
    def __init__(
        self,
        toc: sparse.csr_matrix,    # 细胞计数矩阵 (genes × cells)
        tod: sparse.csr_matrix,    # 全部液滴矩阵 (genes × all_droplets)
        *,
        soup_profile: np.ndarray | None = None,
        clusters: np.ndarray | None = None,
    ):
        ...
    
    def auto_est_cont(self, **kwargs) -> float:
        ...
    
    def adjust_counts(self, method="subtraction", **kwargs) -> sparse.csr_matrix:
        ...
```

面向 AnnData 的高层 API（`soupx.load_10x` / `soupx.decontaminate`）可以考虑内部把 AnnData 转换到 `SoupChannel` 再调用 Rust。

---

### 7. 并行化章节遗漏 GIL 释放策略（原文 4.4）

原文的并行化表缺少 GIL 策略列，但这是 PyO3 + rayon 组合中最关键的工程决策。

**建议在 4.4 节增加 GIL 策略列**（与建议 #3 联动）：

| 函数 | 并行维度 | 工具 | GIL 策略 | 预期利用率 |
|------|---------|------|:---:|:---:|
| `quickMarkers` | 按簇 | `rayon::par_iter` | 释放 GIL（自有数据） | 高 |
| `estimateNonExpressing` | 按基因集 | `rayon::par_iter` | 释放 GIL（自有数据） | 中 |
| `autoEstCont` 密度 | 按 rhoProbe | `rayon::par_iter` | 释放 GIL（自有数据） | 极高 |
| `adjustCounts` | 按细胞 | `rayon::par_bridge` | 释放 GIL（自有数据） | 极高 |
| `expandClusters` | 按簇 | `rayon::par_iter` | 释放 GIL（自有数据） | 高 |

同时标注 `par_bridge` 有非零开销，如果矩阵能预先按行分块，`par_chunks` 更优。

---

### 8. 开发路线图偏乐观（原文 8）

| 周次 | 原文 | 建议 |
|:---:|------|------|
| 1 | 项目脚手架 + CI | 建议 **1.5 周**：PyO3 + maturin 首次配置的坑不少，特别涉及跨平台 wheel |
| 2 | 稀疏桥接 + 全统计库 + 对账 R | **拆为两周**：这是最硬核的工作。第一周桥接，第二周统计函数 |
| 5-6 | autoEstCont 完整实现 | 保留 2 周 |
| 新增 | **Buffer 周** | 在第 5-6 周后插入 1 周 buffer，应对预料之外的数值精度问题 |
| **总计** | **10 周** | **11 周**，更现实 |

---

### 9. 对账策略需增加 Level 0（原文 8）

原文 4 级对账从 soup profile 开始，缺少统计分布函数的单元级对账。这是 statrs 风险的**实际缓解手段**。

**建议增加**：

```
Level 0: 统计分布函数    ppois / dpois / phyper / dgamma / pgamma
                        与 R 的原生函数逐值对账
                        - 常用参数范围（λ=0.1~1000, df=1~100）
                        - 极端值（λ→0, p→0, p→1）
                        - 容差 1e-8
                        输出对账报告，决定是否自写统计函数

Level 1: soup profile    → 完全相同
Level 2: quickMarkers    → tfidf identical, qval within 1e-8
Level 3: ρ estimate      → within 1e-6
Level 4: adjusted counts → >99.9% entries identical
```

---

### 10. statrs 精度风险的实操缓解（原文 9）

原文在风险表写了"必要时自写纯 Rust 统计函数"。**建议补充具体方案**：

- **Week 2 第一步**：编写 `statrs_vs_R_stats.py` 对账脚本，覆盖常用参数空间 + 极端值
- **如果 statrs 达标**：直接使用，记录验证报告
- **如果 statrs 不达标**：不从头实现，优先考虑：
  1. 绑定 [libR-sys](https://github.com/extendr/libR-sys) 直接调用 R 的 C 数学库（Rmath），精度 100% 一致
  2. 移植 [dqrng](https://www.dqrng.org/) 的分布实现（纯 C++，已验证）
  3. 最后才考虑从零自写

---

### 11. 模块文件树结构修正（原文 4.1）

```rust
// 原文：
src/
├── adjustment.rs       // 这不能同时是文件又包含子模块
│   ├── subtraction.rs
│   ├── soup_only.rs
│   └── multinomial.rs
```

Cargo 的模块规则是：文件 + 同名目录不可共存。应改为：

```rust
src/
├── adjustment/
│   ├── mod.rs          // adjust_counts 总入口
│   ├── subtraction.rs
│   ├── soup_only.rs
│   └── multinomial.rs
```

或者在 `adjustment.rs` 中内联所有子逻辑（不推荐，文件会很长）。

---

### 12. `par_bridge` 可优化为 `par_chunks`（原文 4.4）

`adjustCounts` 用 `rayon::par_bridge` 处理 CSR 逐列迭代器，但 `par_bridge` 会在内部加 Mutex，有非零开销。

**优化**：如果能根据 CSR 的 `indptr` 将列预先分成等量的大块（chunks），就可以用 `rayon::par_chunks` 获得更干净的并行：

```rust
let n_cells = indptr.len() - 1;
let chunk_size = (n_cells / rayon::current_num_threads()).max(1);
indptr.par_chunks(chunk_size).for_each(|chunk| { ... });
```

对大型数据集（50K+ cells）此优化有意义，可避免重复 partition 开销。

---

### 13. AnnData 保护逻辑（原文 5.2）

`adata.raw` 可能为 `None`（用户跳过 `load_10x` 或使用通用入口）。所有依赖 `adata.raw.X` 的函数应在入口处检查：

```python
def _check_raw_data(adata: AnnData) -> None:
    """确保 adata.raw 存在且包含合法计数矩阵。"""
    if adata.raw is None:
        raise ValueError(
            "adata.raw is None. Call soupx.load_10x() first, "
            "or pass the droplet matrix via soupx.SoupChannel()."
        )
    if not sparse.issparse(adata.raw.X):
        raise ValueError("adata.raw.X must be a scipy sparse matrix.")
```

---

### 14. multinomial 加速比预估偏高（原文 6.2）

原文标注 multinomial 加速比 **20-40×**，但 BinaryHeap 的常数因子较大（每次 pop/push 涉及 log k 次比较和交换），且循环迭代次数完全取决于数据收敛速度。

**建议**：保守下调为 **15-25×**（与 soupOnly 同档）。实际 benchmark 后再修正。

---

### 15. 超大数据的基因预过滤（原文 7.4）

原文列出了 4 种超大方案，但未提及：**基因预过滤其实已在 R 版的内部实现中**（`geneGiven` 参数允许限制使用的基因子集）。应确认 Rust 版是否保留这一行为，并保持一致。

---

### 16. `quick_markers` 内部与公开 API 的分层（原文 5.3）

`quick_markers` 既作为公开 API（返回 DataFrame 给用户），也作为 `autoEstCont` 的内部调用。这两个场景的数据需求不同：

```python
# 公开 API → 接收 AnnData
def quick_markers(adata: AnnData, *, n: int = 10, ...) -> pd.DataFrame:
    ...

# 内部调用 → 直接接收 CSR 矩阵（避免重复 AnnData → numpy 转换）
def _quick_markers_inner(
    toc: sparse.csr_matrix,
    clusters: np.ndarray,
    soup_profile: np.ndarray,
    *,
    tfidf_min: float = 1.0,
    ...
) -> dict[int, list[tuple[int, float]]]:  # cluster_id → [(gene_idx, score)]
    ...
```

建议明确分开，公开 API 做薄封装。

---

### 17. 增加 Logging & Progress 章节

对于可能耗时数分钟的操作（`autoEstCont`、`adjustCounts`），Python 用户期望看到进度。建议增加一个章节说明：

```python
# 方案 A: tqdm 集成
from tqdm import tqdm

def auto_est_cont(adata, *, show_progress=True, **kwargs):
    if show_progress:
        # 在 Rust 函数中通过回调报告进度
        ...

# 方案 B: logging + verbose 控制
import logging
logger = logging.getLogger("soupx")

def auto_est_cont(adata, *, verbose=True, **kwargs):
    if verbose:
        logger.info("Running quickMarkers...")
        ...
        logger.info(f"Found {n_markers} marker genes")
        logger.info("Calculating posterior densities...")
        ...
        logger.info(f"Estimated rho = {rho:.4f}")
```

推荐用 tqdm + logging 双轨（tqdm 控制台进度条，logging 记录详细日志）。

---

## 小建议汇总

| # | 建议摘要 | 类型 |
|---|---------|:---:|
| 修正 1 | alloc 的 quickselect 声明需要修正 | 🔴 必须 |
| 修正 2 | decontaminate 文档字符串与实际流程对齐 | 🔴 必须 |
| 修正 3 | 补全 GIL 安全策略（推荐 clone 而非零拷贝） | 🔴 必须 |
| 修正 4 | adjust_counts 明确 in-place + 返回自身语义 | 🟡 建议 |
| 5 | load_10x 增加 soup_profile_only 模式 | 🟡 建议 |
| 6 | 增加 SoupChannel 通用 API | 🟡 建议 |
| 7 | 并行化表增加 GIL 策略列 | 🟡 建议 |
| 8 | 开发路线图延长至 11 周 | 🟢 可选 |
| 9 | 对账增加 Level 0（统计函数单测） | 🟡 建议 |
| 10 | statrs 风险具体化（优先 libR-sys） | 🟡 建议 |
| 11 | 修正模块文件树结构（adjustment/mod.rs） | 🟢 可选 |
| 12 | par_bridge → par_chunks 优化 | 🟢 可选 |
| 13 | 增加 adata.raw 为 None 的保护 | 🟡 建议 |
| 14 | multinomial 加速比下调至 15-25× | 🟢 可选 |
| 15 | 确认基因预过滤行为与 R 版一致 | 🟢 可选 |
| 16 | quick_markers 内部分层 | 🟢 可选 |
| 17 | 增加 Logging & Progress 章节 | 🟡 建议 |

---

## 结论

计划**方向完全正确**，无需推翻重来。核心修正（#1-4）应先解决，再进入实现。其余建议可根据开发节奏灵活取舍。
