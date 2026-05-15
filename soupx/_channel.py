"""SoupChannel — low-level API for non-10X data sources.

Mirrors the R SoupX::SoupChannel class. Useful when the user has
count data not from 10X CellRanger (e.g., inDrops, Drop-seq, Smart-seq2).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sparse


class SoupChannel:
    """Core SoupX object for arbitrary count matrices.

    Accepts raw droplet and cell matrices directly, bypassing the
    10X CellRanger parsing layer.

    Parameters
    ----------
    toc : sparse.csr_matrix
        Table of Counts — cell-only count matrix (genes × cells).
    tod : sparse.csr_matrix
        Table of Droplets — all droplet count matrix (genes × droplets).
    soup_profile : np.ndarray, optional
        Pre-computed soup proportions. If None, estimated from empty droplets.
    clusters : np.ndarray, optional
        Cluster labels per cell.
    """

    def __init__(
        self,
        toc: sparse.csr_matrix,
        tod: sparse.csr_matrix,
        *,
        soup_profile: np.ndarray | None = None,
        clusters: np.ndarray | None = None,
    ):
        self.toc = toc
        self.tod = tod
        self.soup_profile = soup_profile
        self.clusters = clusters

        if soup_profile is None:
            self.estimate_soup()

    def estimate_soup(self, soup_range: tuple[int, int] = (0, 100)) -> None:
        """Estimate soup profile from empty droplets."""
        # TODO: implement
        raise NotImplementedError("SoupChannel.estimate_soup not yet implemented")

    def auto_est_cont(self, **kwargs) -> float:
        """Estimate contamination fraction ρ automatically."""
        # TODO: implement
        raise NotImplementedError("SoupChannel.auto_est_cont not yet implemented")

    def adjust_counts(self, method: str = "subtraction", **kwargs) -> sparse.csr_matrix:
        """Remove contamination from cell counts.

        Returns corrected count matrix (does not modify toc in-place).
        """
        # TODO: implement
        raise NotImplementedError("SoupChannel.adjust_counts not yet implemented")

    @property
    def n_genes(self) -> int:
        return self.toc.shape[0]

    @property
    def n_cells(self) -> int:
        return self.toc.shape[1]

    @property
    def n_droplets(self) -> int:
        return self.tod.shape[1]
