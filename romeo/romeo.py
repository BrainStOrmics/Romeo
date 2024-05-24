import scanpy as sc
import numpy as np
import pandas as pd
import sys
import logging
from typing import Dict, Any

from .utile import find_markers, adata_add_metadata

def romeo(
    adata,
    groupby: str = 'group',
    key_layer: str = None,
    normalize: bool = True,
    key_added: str = 'romeo',
    merge_mode: str = 'outer',
    angular_consistency: float = 0.1,
    min_positivity_rate: float = 0.0
) -> None:
    """
    Find marker genes based on the consistency between gene expression levels 
    and cell positivity rates.

    Args:
        adata: Annotated data matrix or a list containing multiple anndata matrix.
        groupby: The key of cell groups in adata.obs or a list containing multiple 
            groupbys for each slice. Defaults to group.
        key_layer: The key from `adata.layers` whose value will be used or a list 
            for containing multiple key_layer for each slice. If None, the adata.X 
            will be used. Defaults to None.
        normalize: Normalize the count matrix by sc.pp.log1p(). Defaults to True.
        key_added: The key in `adata.uns` where information is saved. Defaults to romeo.
        merge_mode: The merge mode for multiple slices, intersection (inner) or union (outer).
            Defaults to outer.
        angular_consistency: The weight of angular consistency. Defaults to 0.1.
        min_positivity_rate: The minimum cell positivity rate in group. Defaults to 0.0.

    Returns
    -------
    names: Structured array to be indexed by group id storing the gene names.
    scores: Structured array to be indexed by group id storing scores for each gene.

    Examples
    --------
    >>> import scanpy as sc
    >>> import romeo
    >>> adata_counts = sc.datasets.pbmc3k()
    >>> adata = sc.datasets.pbmc3k_processed()
    >>> adata.layers['counts'] = adata_counts[adata.obs_names, adata.var_names].X.copy()
    >>> romeo.romeo(adata, 'louvain', key_layer='counts')
    >>> romeo.romeo_markers_dotplot(adata, groupby='louvain', dotplot_kwargs={'cmap': 'Spectral_r'})
    """
    if type(adata) != list:
        adata = [adata]

    pd_romeo_markers = find_markers(
        adata=adata,
        groupby=groupby,
        key_layer=key_layer,
        normalize=normalize,
        merge_mode=merge_mode,
        angular_consistency=angular_consistency,
        min_positivity_rate=min_positivity_rate
    )

    adata_add_metadata(
        adata=adata,
        markers=pd_romeo_markers,
        key_added=key_added,
    )

def romeo_markers(
    adata,
    groupby: str = 'annotation',
    key_added: str = 'romeo',
    top: int = 3,
    return_marker_file: str = None,
) -> Dict[Any, str]:
    """
    Get the top n genes for each group based on the romeo score.

    Args:
        adata: Annotated data matrix.
        groupby: The key of cell groups in adata.obs. Defaults to annotation.
        key_added: The key in `adata.uns` where information is saved. Defaults to romeo.
        top: The `top` genes with the most lowest romeo score in each group. Defaults to 3.
        return_marker_file: The file to save marker gene list. Defaults to None.

    Returns:
        Dict[Any, str]: A dict for the top genes, labels and names.
    """
    if key_added in adata.uns:
        pd_conep_markers = adata.uns[key_added]
        if (len(pd_conep_markers)/len(np.unique(adata.obs[groupby])) != len(adata.var_names)):
            pd_conep_markers = pd_conep_markers[
                (~np.isin(pd_conep_markers['names'], 
                    np.setdiff1d(pd_conep_markers['names'], adata.var_names))
                ) &
                (~np.isin(pd_conep_markers['labels'], 
                    np.setdiff1d(pd_conep_markers['labels'], adata.obs[groupby]))
                )
            ]
        pd_conep_markers['labels'] = pd.Categorical(
            pd_conep_markers['labels'], 
            categories=adata.obs[groupby].cat.categories
        )

        pd_top_genes = pd_conep_markers.sort_values(['labels', 'scores'], 
            ascending=[True, True]).groupby('labels').head(top)

        if (return_marker_file != None):
            pd_top_genes.to_csv(return_marker_file, index=None)
            
        return pd_top_genes.groupby('labels')['names'].apply(list).to_dict()
    else:
        logging.basicConfig(format='%(asctime)s %(message)s')
        logging.error(f"The {key_added} does not exists in .uns, process conep first.")
        sys.exit(-1)

def romeo_markers_dotplot(
    adata,
    groupby: str = 'group',
    key_added: str = 'romeo',
    top: int = 3,
    return_marker_file: str = None,
    return_fig_file: str = None,
    dotplot_kwargs: dict = {},
) -> None:
    """
    Make a dot plot for marker genes by using sc.pl.dotplot.

    Args:
        adata: Annotated data matrix.
        groupby: The key of cell groups in adata.obs. Defaults to group.
        key_added: The key in `adata.uns` where information is saved. Defaults to romeo.
        top: The `top` genes with the most highest romeo score in each group. Defaults to 3.
        return_marker_file: The file to save marker gene list. Defaults to None.
        return_fig_file: The file to save dotplot figure. Defaults to None.
        dotplot_kwargs: Other parameters for sc.pl.dotplot. Defaults to {}.
    
    Returns:
        If return_fig_file is not None, return a dotplot figure.
    """
    var_names = romeo_markers(
        adata=adata,
        groupby=groupby,
        key_added=key_added,
        top=top,
        return_marker_file=return_marker_file,
    )

    if return_fig_file != None:
        sc.pl.dotplot(
            adata=adata,
            var_names=var_names,
            groupby=groupby,
            return_fig=True,
            **dotplot_kwargs,
        ).savefig(return_fig_file)
    else:
        sc.pl.dotplot(
            adata=adata,
            var_names=var_names,
            groupby=groupby,
            **dotplot_kwargs,
        )
  