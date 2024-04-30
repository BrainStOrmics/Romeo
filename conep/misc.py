from typing import Tuple, List
import numpy as np
import scanpy as sc
from scipy import sparse
from sklearn.metrics.pairwise import euclidean_distances, cosine_similarity
import pandas as pd

Max_consistency_score = 100

def adata_blocks(
    adata, 
    groupby: str = 'group'
) -> Tuple[bool]:
    """
    Classify adata based on groupby.

    Args:
        adata: Annotated data matrix.
        groupby: The key of cell groups in adata.obs. Defaults to 'group'.

    Returns:
        Tuple: Cell info for each group.
    """
    blocks = np.zeros(
        (len(adata.obs[groupby].cat.categories), adata.obs[groupby].values.size), 
        dtype=bool
    )
    
    for i, key in enumerate(adata.obs[groupby].cat.categories):
        blocks[i] = adata.obs[groupby].cat.categories[i] == adata.obs[groupby].values

    return blocks

def group_stat(
    count_matrix,
    expr_min, 
    expr_max
) -> np.ndarray:
    """
    Matrix statistics.

    Args:
        count_matrix: a DataFrame containing gene expression for each cell.
        expr_min: The minimum value for scaling.
        expr_max: The maximum value for scaling.

    Returns:
        np.ndarray: expression levels and positive rates for each gene in each group. 
    """
    if sparse.issparse(count_matrix):
        expression_levels = (np.mean(count_matrix.A, axis=0)-expr_min)/(expr_max-expr_min)
        positive_rates = count_matrix.getnnz(axis=0) / count_matrix.shape[0]
    else:
        expression_levels = (np.mean(count_matrix, axis=0)-expr_min)/(expr_max-expr_min)
        positive_rates = np.count_nonzero(count_matrix, axis=0) / count_matrix.shape[0]
    
    return expression_levels, positive_rates

def convert_to_mds_and_rts(
    arr
) -> np.ndarray:
    """
    Convert to mean difference and ratio of the array.

    Args:
        arr: A np.ndarray.

    Returns:
        np.ndarray: mean difference and ratio.
    """
    arr_sum = np.sum(arr, axis=0)
    arr_sum = np.where(arr_sum==0, 1e-5, arr_sum)
    arr_mean = np.mean(arr, axis=0)

    return arr-arr_mean, arr/arr_sum

def single_slice_stat(
    adata, 
    groupby: str = 'group', 
    key_layer: str = None,
    normalize: bool = True,
    min_positivity_rate: float = 0.0
) -> np.ndarray:
    """
    Calculate the postivity rates and expression levels for each slice.

    Args:
        adata: Annotated data matrix.
        groupby: The key of cell groups in adata.obs. Defaults to group.
        key_layer: The key from `adata.layers` whose value will be used.
            If None, the adata.X will be used. Defaults to None.
        normalize: Normalize the count matrix by sc.pp.log1p(). Defaults to True.
        min_positivity_rate: The minimum cell positivity rate in group. Defaults to 0.0.
    
    Returns:
        np.ndarray: positivity rates for each gene in all groups [n_groups, n_genes]
    """
    blocks = adata_blocks(adata, groupby=groupby)
    n_groups, n_genes = blocks.shape[0], adata.X.shape[1]
    expression_levels = np.zeros((n_groups, n_genes))
    positive_rates = np.zeros((n_groups, n_genes))

    if key_layer in adata.layers:
        tmp_X = adata.X.copy()
        adata.X = adata.layers[key_layer].copy()
    if normalize:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    expr_min, expr_max = adata.X.min(), adata.X.max()
    expression_levels_mean, positive_rates_mean = group_stat(adata.X, expr_min, expr_max)

    for i, block in enumerate(blocks):
        count_matrix_block = adata.X[block]
        expression_levels[i], positive_rates[i] = group_stat(count_matrix_block, expr_min, expr_max)

    expression_levels_mean_diff, expression_levels_ratio = convert_to_mds_and_rts(expression_levels)
    positive_rates_mean_diff, positive_rates_ratio = convert_to_mds_and_rts(positive_rates)

    # convert to a 4 columns array = [n_groups * n_genes, 4]
    gene_info = np.zeros((n_groups*n_genes, 4))
    for i, block in enumerate(blocks):
        gene_info[i*n_genes:(i+1)*n_genes] = np.concatenate((
            expression_levels_mean_diff[i][:, None],
            expression_levels_ratio[i][:, None],
            positive_rates_mean_diff[i][:, None],
            positive_rates_ratio[i][:, None]), axis=1
        )

    # recover the adata.X if necessary
    if key_layer in adata.layers:
        adata.X = tmp_X

    # filter genes with lowly positivity rates
    if min_positivity_rate > 0:
        group_names = adata.obs[groupby].cat.categories
        gene_names = adata.var_names
        filtered = np.where(positive_rates < min_positivity_rate)
        genes_filtered = pd.DataFrame({'labels': group_names[filtered[0]], 'names': gene_names[filtered[1]]})
        genes_filtered.set_index(['labels', 'names'], inplace=True)
    else:
        genes_filtered = pd.DataFrame()

    return gene_info, genes_filtered

def score(
    X, 
    angular_consistency: float = 0.1
) -> np.ndarray:
    """
    Calculate the consistency score.

    Args:
        X: A 4 columns array where each row is a gene.
        angular_consistency: The weight of angular consistency. Defaults to 0.1.

    Returns:
        np.ndarray: Consistency score for each gene in each group.
    """
    Y = np.ones(X.shape[1]).reshape(1, -1)
    consistency_score = euclidean_distances(X, Y) + \
        np.arccos(cosine_similarity(X, Y)) * angular_consistency
    
    return np.round(consistency_score.flatten(), decimals=3)

def output(
    scores,
    names,
    labels,
    genes_filtered
) -> pd.DataFrame:
    """
    Generate the output of conep.

    Args:
        scores: A np.ndarray for scores.
        names: A list for gene names.
        labels: A list of group for each gene.
        genes_filtered: Genes to be filtered.

    Returns:
        pd.DataFrame: A `pandas.DataFrame` contains labels, names and scores.
    """
    if len(genes_filtered) == 0:
        return pd.DataFrame({'labels': labels, 'names': names, 'scores': scores})
    else:
        pd_genes = pd.DataFrame({'labels': labels, 'names': names, 'scores': scores})
        pd_genes.set_index(['labels', 'names'], inplace=True)
        pd_genes.loc[genes_filtered.index] = Max_consistency_score
 
        return pd_genes.reset_index()

def parms_to_list(
    arr,
    n
) -> List:
    """
    Repeat parameters based on the length.

    Args:
        arr: A list for parameter.
        n: The length of returned parameters.

    Returns:
        List: A list of returned parameters.
    """
    if arr == None:
        return [None] * n
    else:
        if type(arr) != list:
            arr = [arr]
        if len(arr) < n:
            return repeat_n(arr[0], n)
        else:
            return arr[:n]

def repeat_n(
    arr,
    n: int = 1
) -> List:
    """
    Repeat each element n times.

    Args:
        arr: An array or list to be repeated.
        n: Repeat times. Defaults to 1.

    Returns:
        List: The repeated list.
    """

    return [item for s in arr for item in [s]*n]

def single_slice(
    adata,
    groupby: str = 'group',
    key_layer: str = None,
    normalize: bool = True,
    angular_consistency: float = 0.1,
    min_positivity_rate: float = 0.0
) -> pd.DataFrame:
    """
    Detect the marker genes for single slice.

    Args:
        adata: Annotated data matrix.
        groupby: The key of cell groups in adata.obs. Defaults to group.
        key_layer: The key from `adata.layers` whose value will be used.
            If None, the adata.X will be used. Defaults to None.
        normalize: Normalize the count matrix by sc.pp.log1p(). Defaults to True.
        angular_consistency: The weight of angular consistency. Defaults to 0.1.
        min_positivity_rate: The minimum cell positivity rate in group. Defaults to 0.0.

    Returns:
        pd.DataFrame: A `pandas.DataFrame` contains labels, names and scores.
    """
    groups = adata.obs[groupby].cat.categories
    np_gene_info, genes_filtered = single_slice_stat(
        adata=adata,
        groupby=groupby,
        key_layer=key_layer,
        normalize=normalize,
        min_positivity_rate=min_positivity_rate
    )

    return output(
        scores=score(np_gene_info, angular_consistency), 
        names=np.array(tuple(adata.var_names)*len(groups)),
        labels=repeat_n(groups,adata.shape[1]),
        genes_filtered=genes_filtered
    )

def multiple_slices(
    list_adata,
    list_groupby,
    list_key_layer,
    normalize: bool = True,
    angular_consistency: float = 0.1,
    min_positivity_rate: float = 0.0,
    merge_mode: str = 'outer',
) -> pd.DataFrame:
    """
    Detect the marker genes for multiple slices.

    Args:
        list_adata: A list containing annotated data matrix for each slice.
        list_groupby: A list containing the key of cell groups in adata.obs 
            for each slice.
        list_key_layer: A list containing the key from `adata.layers` 
            whose value will be used.
        normalize: Normalize the count matrix by sc.pp.log1p(). Defaults to True.
        angular_consistency: The weight of angular consistency. Defaults to 0.1.
        min_positivity_rate: The minimum cell positivity rate in group. Defaults to 0.0.
        merge_mode: Merge multiple slices based on intersection (inner) or 
            union (outer). Defaults to outer.

    Returns:
        pd.DataFrame: A `pandas.DataFrame` contains labels, names and scores.
    """
    pd_genes = pd.DataFrame()
    genes_filtered = pd.DataFrame()

    for i, adata in enumerate(list_adata):
        slice_groups = adata.obs[list_groupby[i]].cat.categories
        np_gene_info, slice_genes_filtered = single_slice_stat(
            adata=adata,
            groupby=list_groupby[i],
            key_layer=list_key_layer[i],
            normalize=normalize,
            min_positivity_rate=min_positivity_rate
        )
        slice_names=np.array(tuple(adata.var_names)*len(slice_groups))
        slice_labels=repeat_n(slice_groups, adata.shape[1])
        pd_slice = pd.concat([
            pd.DataFrame(np_gene_info),
            pd.DataFrame({'names': slice_names, 'labels':slice_labels})], axis=1)
        pd_slice.set_index(['names', 'labels'], inplace=True)
        pd_genes = pd.concat([pd_genes, pd_slice], axis=1, join=merge_mode)
        genes_filtered = pd.concat([genes_filtered, slice_genes_filtered], axis=1, join=merge_mode)
    
    return output(
        scores=score(pd_genes.to_numpy(), angular_consistency),
        names=pd_genes.index.get_level_values(0),
        labels=pd_genes.index.get_level_values(1),
        genes_filtered=genes_filtered
    )
