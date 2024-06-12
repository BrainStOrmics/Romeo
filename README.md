# A Robust Marker Gene Selection by Harmonizing Expression Levels and Positivity Rates in Single-Cell Resolution Transcriptome
[![python~=3.8](https://img.shields.io/badge/python-3.8-brightgreen)](https://www.python.org/)
[![License: GPL3.0](https://img.shields.io/badge/License-GPL3.0-yellow)](https://opensource.org/license/gpl-3-0/)

Cell type classification is a crucial stage in single-cell and spatial transcriptome analysis. Romeo, which stands for **RO**bust **M**arker identifier with **E**xpression level and positive rati**O**, is a valuable tool for precisely and reliably identifying marker genes. Its efficiency in identifying biologically significant marker genes quickly makes it an indispensable asset for enhancing our comprehension of cellular diversity and function.

![image](assets/workflow_of_romeo.png)

# Installation

```
pip install git+https://github.com/BrainStOrmics/Romeo.git
```

# Usage

The [Romeo tutorials](https://github.com/BrainStOrmics/Romeo/tree/main/tutorials/romeo_tutorials.ipynb) provides a quick-start guide based on the pbmc3k dataset.
## parameters

options | description
---- | ----
adata | Annotated data matrix or a list containing multiple anndata matrix.
groupby | The key of cell groups in adata.obs or a list containing multiple groupbys for each slice. Defaults to group.
key_layer | The key from `adata.layers` whose value will be used or a list for containing multiple key_layer for each slice. If None, the adata.X will be used. Defaults to None.
normalize | Normalize the count matrix by sc.pp.log1p(). Defaults to True.
key_added | The key in `adata.uns` where information is saved. Defaults to romeo.
merge_mode | The merge mode for multiple slices, intersection (inner) or union (outer). Defaults to outer.
angular_consistency | The weight of angular consistency. Defaults to 0.1.
min_positivity_rate | The minimum cell positive ratio in each group. Defaults to 0.0.


# Enviroments
- python>=3.8.0
- numpy>=1.18.0
- pandas>=0.25.1
- scanpy>=1.9.0
- scipy>=1.9.0
- anndata>=0.8.0
- scikit-learn>=0.19.0

# Question

For questions about the code and tutorial, please contact Qianhua ZHU, zhuqianhua@genomics.cn.

# Citation
If Romeo is useful for your research, please consider citing.