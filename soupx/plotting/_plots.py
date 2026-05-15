"""Plotting functions for soupX diagnostics."""

from __future__ import annotations

from typing import Any
from anndata import AnnData


def plot_marker_distribution(
    adata: AnnData,
    gene_list: dict[str, list[str]],
    **kwargs: Any,
) -> None:
    """Plot the observed vs expected distribution of marker genes.

    Mirrors R SoupX::plotMarkerDistribution.
    """
    # TODO: implement
    raise NotImplementedError("plot_marker_distribution not yet implemented")


def plot_marker_map(
    adata: AnnData,
    gene_set: str,
    **kwargs: Any,
) -> None:
    """Plot expression of marker genes on a DR embedding.

    Mirrors R SoupX::plotMarkerMap.
    """
    # TODO: implement
    raise NotImplementedError("plot_marker_map not yet implemented")


def plot_change_map(
    adata: AnnData,
    gene_set: str,
    **kwargs: Any,
) -> None:
    """Plot expression change (before vs after adjustment) on a DR embedding.

    Mirrors R SoupX::plotChangeMap.
    """
    # TODO: implement
    raise NotImplementedError("plot_change_map not yet implemented")
