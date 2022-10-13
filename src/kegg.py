#!/usr/bin/env python 

import pandas as pd 
import tqdm

def kegg_compounds(): 
    """ 
    Gets IDs of all KEGG compounds with their name 

    RETURNS 
    -------
    :class:`pd.DataFrame`
        Pandas Dataframe with column of KEGG ID and name of compound 
    """
    # Get compound list from https://rest.kegg.jp/list/compound
    df = pd.read_csv('https://rest.kegg.jp/list/compound', sep = '\t', names = ['id', 'name'])
    # ids are returned in the form cpd:ID (e.g. cpd:C00931). Remove leading cpd:
    df['id'] = df['id'].apply(lambda x: x.split(':')[1])
    # Only store first name 
    df['name'] = df['name'].apply(lambda x: x.split(';')[0])
    
    return df 

    
def keggID_from_molformula(molecular_formula:str):
    """ 
    Uses KEGG API to return a dataframe that stores all KEGG IDs whose molecular formulas match a query. 

    INPUT
    ----- 
    molecular_formula (`str`)
        Molecular formula of molecule of interest

    RETURNS 
    ------- 
    :class:`pd.DataFrame`

    EXAMPLE 
    ------- 
    >>> keggID_from_molformula('H2O')
       id molecular_formula
    0  C00001               H2O
    """

    # KEGG API 
    request = f'https://rest.kegg.jp/find/compound/{molecular_formula}/formula'

    # Use web-retrieval tool of pandas and read in as dataframe 
    df = pd.read_csv(request, sep='\t', names = ['id', 'molecular_formula'])

    # Only consider entries with exact match of molecular formula  
    df = df.loc[df['molecular_formula'] == molecular_formula]

    # ids are returned in the form cpd:ID (e.g. cpd:C00931). Remove leading cpd:
    df['id']= df['id'].apply(lambda x: x.split(':')[1])
    
    return df 


def kegg_batch_retrieval(molecular_formulas:List[str]): 
    """ 
    Retrieves information for several molecular formulas that are passed as list and returns a single dataframe 

    EXAMPLE 
    -------
    >>> molformulas = ['H2O', 'C6H12O6', 'C6H8O7']
    >>> df = kegg_batch_retrieval(molformulas)
    >>> df.sample(3)
            id molecular_formula
    23  C01906           C6H12O6
    46  C21050           C6H12O6
    7   C00267           C6H12O6
    
    """
    from tqdm import tqdm
    
    df_list = []
    for i, mol in zip(tqdm(range(len(molecular_formulas))), molecular_formulas): 
        df_list.append(keggID_from_molformula(mol)) 

    return pd.concat(df_list, axis=0)
