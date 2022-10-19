#!/usr/bin/env python 

import pandas as pd 
from tqdm import tqdm
from typing import List


def all_compounds(): 
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

    
def id_from_molformula(molecular_formula:str, molecular_formula_column = 'molecular_formula'):
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
    df = pd.read_csv(request, sep='\t', names = ['id', molecular_formula_column])

    # Only consider entries with exact match of molecular formula  
    df = df.loc[df['molecular_formula'] == molecular_formula]

    # ids are returned in the form cpd:ID (e.g. cpd:C00931). Remove leading cpd:
    df['id']= df['id'].apply(lambda x: x.split(':')[1])
    
    return df 


def batch_retrieval(molecular_formulas:List[str],  molecular_formula_column = 'molecular_formula'): 
    """ 
    Retrieves information for several molecular formulas that are passed as list and returns a single dataframe 

    EXAMPLE 
    -------
    >>> molformulas = ['H2O', 'C6H12O6', 'C6H8O7']
    >>> df = kegg_batch_retrieval(molformulas)
    >>> df.sample(3)
            id molecular_formula
    1   C00001           H2O
    46  C21050           C6H12O6
    7   C00267           C6H12O6
    
    """    
    df_list = []
    for mol in tqdm(molecular_formulas, total = len(molecular_formulas)): 
        df_list.append(id_from_molformula(mol,  molecular_formula_column)) 

    return pd.concat(df_list, axis=0).reset_index(drop = True)

def calculate_redundancy(df:pd.DataFrame, molecular_formula_column:str = 'molecular_formula'): 
    """ 
    Calculates number of molecules with same molecular formula. 
    """
    counts = df.value_counts(molecular_formula_column).reset_index().rename(columns = {0:'redundancy'})
    return df.merge(counts, on = molecular_formula_column)  
    

def annotate_data(df:pd.DataFrame, molecular_formula_column = 'molecular_formula'): 
    """ 
    Returns fully annotated dataframe. This function concatenates the functions `kegg.batch_retrieval`, `kegg.all_compounds` with some basic pandas functionalities for convenience. 
    """
    molecular_formulas = df[molecular_formula_column].to_list()
    df_keggIDs = batch_retrieval(molecular_formulas)
    df_all_compounds = all_compounds()

    df_annotation = pd.merge(df_keggIDs, df_all_compounds, how = 'left', on = 'id')

    return pd.merge(df, df_annotation, how = 'outer', on = molecular_formula_column) 
    

