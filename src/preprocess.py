#!/usr/bin/env python 

def parse_ranked_genes(df, pcutoff:float =1e-3, changecutoff:None|float = None, ion_column:str = 'names', logchange_column:str = 'logfoldchanges', pval_column:str = 'pvals_adj'): 
    """ 
    df (:class:`pd.DataFrame`)
        Dataframe which is supposed to be parsed
    pcutoff (float)
        p value cutoff (default 0.001). If no cutoff is desired, set pcutoff = 1 
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
    df['molecular_formula'] = df[ion_column].str.extract('(.*)[\+|\-].*') 

    df = df.loc[df[pval_column] < pcutoff]

    if changecutoff is not None: 
        df = df.loc[df[pval_column] > changecutoff]
    
    return df 

