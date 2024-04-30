import pandas as pd
import scanpy as sc
import sys
import logging
from scipy import sparse
from typing import Any, List

from .misc import multiple_slices, single_slice, parms_to_list

def find_markers(
    adata,
    groupby: str = 'group',
    key_layer: str = None,
    normalize: bool = True,
    merge_mode: str = 'outer',
    angular_consistency: float = 0.1,
    min_positivity_rate: float = 0.0
) -> pd.DataFrame:
    """
    Find the marker genes based on the consistency between gene expression and 
    cell positivity rates.

    Args:
        adata: A list containing multiple anndata matrix for each slice.
        groupby: The key of cell groups in adata.obs or a list containing multiple 
            groupbys for each slice. Defaults to group.
        key_layer: The key from `adata.layers` whose value will be used or a list 
            for containing multiple key_layer for each slice. If None, the adata.X 
            will be used. Defaults to None.
        normalize: Normalize the count matrix by sc.pp.log1p(). Defaults to True.
        merge_mode: merge multiple slices based on intersection (inner) or union (outer).
            Defaults to outer.
        angular_consistency: The weight of angular consistency. Defaults to 0.1.
        min_positivity_rate: The minimum cell positivity rate in group. Defaults to 0.0.

    Returns:
        pd.DataFrame: a `pandas.DataFrame` contains labels, names and scores
    """
    logging.basicConfig(format='%(asctime)s %(message)s')

    n_adata = len(adata)
    list_groupby = parms_to_list(groupby, n_adata)
    list_key_layer = parms_to_list(key_layer, n_adata)

    if n_adata == 1:
        pd_markers = single_slice(
            adata=adata[0],
            groupby=list_groupby[0],
            key_layer=key_layer,
            normalize=normalize,
            angular_consistency=angular_consistency,
            min_positivity_rate=min_positivity_rate
        )
    else:
        pd_markers = multiple_slices(
            list_adata=adata,
            list_groupby=list_groupby,
            list_key_layer=list_key_layer,
            normalize=normalize,
            angular_consistency=angular_consistency,
            min_positivity_rate=min_positivity_rate,
            merge_mode=merge_mode
        )
    pd_markers['labels'] = pd_markers['labels'].astype(type(adata[0].obs[list_groupby[0]][0]))

    return pd_markers

def adata_add_metadata(
    adata,
    markers,
    key_added: str = 'conep',
) -> None:
    """
    Add conep markers to .uns[key_added].

    Args:
        adata: A list containing annotated data matrix.
        markers: The marker gene list detected by conep.
        key_added: The key in `adata.uns` where information is saved. Defaults to conep.

    Returns:
        None
    """

    for i in range(len(adata)):
        adata[i].uns[key_added] = markers
        
