#!/usr/bin/env python 

import pandas as pd 
import scanpy as sc 
import re 

def rank_entity_groups(adata, groupby, method = 'wilcoxon', use_raw = False, copy = False, **kwargs) -> None: 
    """ 
    Ranks entities (e.g. metabolites) by the impact they have on the differentiation of a group from another. If

    INPUT
    ----- 
    adata (:class:`AnnData`)
        AnnData Object
    groupby (str)
        Group/Property that is used to define different states (e.g. Cluster of UMAP or experimental condition). Any text from adata.obs_keys() is in principle valid.
    method (str) 
        Method for ranking
    use_raw (bool)
        Whether to use raw values or not
    copy
        Whether to copy the AnnData object 
    kwargs
        Valid arguments of sc.tl.rank_gene_groups()
    """
    if copy: 
        anndata = copy.deepcopy(adata) 
        return  sc.tl.rank_genes_groups(anndata, groupby, method, use_raw, kwargs)
    else: 
        sc.tl.rank_genes_groups(adata, groupby, method, use_raw, kwargs)
    

def parse_ranked_ions(df:pd.DataFrame, pcutoff:None|float = None, changecutoff:None|float = None, ion_column:str = 'names', molecular_formula_column:str = 'molecular_formula', adduct_column:str = 'adduct', logchange_column:str = 'logfoldchanges', pval_column:str = 'pvals_adj'): 
    """ 
    df (:class:`pd.DataFrame`)
        Dataframe which is supposed to be parsed
    pcutoff (float)
        p value cutoff. 
    changecutoff (None|float)
        Cutoff for change in abundance (default: None/Not performed)
    ion_column (str)
        Name of column that contains ion formulas (default: names) 
    logchange_column (str)
        Name of column that contains logfold changes of metabolites (default: logfoldchanges) 
    pval_column (str)
        Name of column that contains p values of changes of metabolites (default: pvals_adj) 
   
    EXAMPLE
    -------
    >>> df
                   names     scores  logfoldchanges          pvals      pvals_adj
    121    C10H10N2O3S-H   3.120556       -0.096448   1.805101e-03   2.758008e-03
    241  C11H11ClN2O2+Cl  -2.799665       -6.435749   5.115568e-03   7.651858e-03
    17         C11H9N3-H  29.209204        0.578870  1.481648e-187  2.776140e-186
    >>> parse_ranked_genes(df)
    
    """

    # Extract molecular formulas 
    df[molecular_formula_column] = df[ion_column].str.extract('(.*)[\+|\-].*') 
    df[adduct_column] = df[ion_column].str.extract('.*[\+|\-](.*)') 


    if pcutoff is not None: 
        df = df.loc[df[pval_column] < pcutoff]

    if changecutoff is not None: 
        df = df.loc[df[logchange_column] > changecutoff]
    
    return df 



def ion_to_molformula(ion_formula:str) -> str:
    """ 
    Returns molecular formula from ion formulas of the form
    MOLFORMULA[+/-]ION

    EXAMPLE 
    -------
    >>> print(ion_to_molformula(C5H5N5-H))
    ... C5H5N5

    >>> df['molecular_formula'] = df['ion'].apply(lambda x: ion_to_molformula(x))
    """
    return re.findall('(.*)[\+|\-].*', ion_formula)[0]
