#!/usr/bin/env python 

import requests
from typing import List 

def make_ipath_selection(ids:List[str], colors:None|str|List[str] = None, sizes:None|int|List[int] = None, highlight:None|List[dict] = None, save:None|str = None): 
    """ 
    Returns input for ipath3 for a number of ids. Coloring and widths can be specified. Specific nodes or pathways can be highlighted in the keyword highlight.  

    INPUT
    ----- 

    ids (List[str]):
        List of KEGG IDs that are supposed to be highlighted in the map
    colors (None|str|List[int]):
        Colors of the highlighted ids. Valid color formats are HEX (#XXXXXX), RGB (RGB(XXX, YYY, ZZZ)) and CYMK. If None, everything is highlighted red, if str, the respective color is used for all ids. 
        If list, every id is colored according to the specified color. 
    sizes (None|int|List[int]): 
        Size of highlighted ids (px). If None, everything is set to size 10px, if int, the respective size is used for all ids. 
        If list, every id is set to the respective size. 
    highlight (None|Dict): 
        Highlight specific KEGG IDs in map (e.g. pathways) in json-format [{keggID: <KEGG ID>, color:<color>, width:<width>}]
        Glycolysis: 00010 (hsa00010)
        TCA cycle: 00020
        Pentose Phosphate Pathway: 00030 (hsa00030)
        Fatty acid degradation: 00061
        Fatty acid degradation: 00071

    
    """
    # Set defaults 
    if colors is None: 
        colors = ['#FF0000']*len(ids)
    if type(colors) is str:
        colors = [colors]*len(ids)
    if type(colors) is list and len(colors) != len(ids): 
        raise ValueError('Colors must have same length as ids')

    if sizes is None: 
        sizes = [10]*len(ids)
    if type(sizes) is int:
        sizes = [sizes]*len(ids)
    if type(sizes) is list and len(sizes) != len(ids): 
        raise ValueError('Sizes must have same length as ids')



    # Make main seleciton based on input lists 
    ipath_selection = ''
    for ipathID, color, size in zip(ids, colors, sizes): 
        ipath_selection += f'{ipathID} {color} W{size}\n'
    

    if type(highlight) is list: 
        for item in highlight: 
            ipathID, color, size = item.get('keggID'), item.get('color'), item.get('size')

            if ipathID is None: 
                raise ValueError('All ids must be defined')
            if color is None: 
                color = '#cccccc'
            if size is None: 
                size = 10

            ipath_selection += f'{ipathID} {color} W{size}\n'

    # Return/save output
    if save is None: 
        return ipath_selection
    else: 
        with open(save, 'w') as f: 
            f.write(ipath_selection)
            return f'Saved in {save}'
    

def post(selection = '', 
               default_opacity = 1, 
               default_edge_width = 3, 
               default_node_radius = 7, 
               keep_colors = 0, 
               default_color = '#cccccc', 
               background_color = '#ffffff', 
               whole_pathways = 0, 
               whole_modules = 0, 
               query_reactions = 0, tax_filter = 9606, metabolic_map = 'metabolic', export_type = 'SVG') -> bytes:
    """
    Posts input to [ipath3 server](https://pathways.embl.de/tools.cgi) (HTTPS:POST server: https://pathways.embl.de/mapping.cgi) and returns image with highlighted pathways.
    Keywords are the same as in the online version.  
    Since the utility is currently disabled on the website, this method cannot be used at the moment. 

    INPUT
    ----- 
    selection (str)
        Selection of highlighted entities in pathway map in ipath3. The selection can have an arbitrary number of rows. 
        Each row has the form 
        `<ID/KEGG ID> <color (HEX #XXXXXX/RGB RGB(X,Y,Z)/CYMK)> <width px>`
    default_opacity (float {0..1})
        Default opacity/alpha value of nodes (default is 1)
    default_edge_width (int)
        Default width of edges in graph (represent reactions)  
    default_node_radius (int)
        Default size of nodes in graph (represent metabolites)
    keep_colors (int {0,1})
        Whether to keep default colors of pathways provided by ipath3 (disabled per default)
    background_color (str)
        Color of background in HEX (#XXXXXX), RGB (RGB(X,Y,Z)) or CYMK code 
    whole_pathways (int, {0,1})
        If enabled, any pathway with at least one matching edge or compound will be highlighted (disabled per default). 
    whole_modules (int, {0,1})
        If enabled, any KEGG module with at least one matching edge or compound will be highlighted (disabled per default).
    query_reactions (int, {0,1})
        If enabled, compound presence within each edges reactions will also be checked (disabled per default).
    tax_filter (int)
        An NCBI tax ID or KEGG 3 letter species code can be provided. Only pathways present in selected species will be included in the map. 
        Per default, only human metabolic pathways are displayed (NCBI species ID: `9606`). Note that either `selection` or `tax_filter` have to be specified. 
        Human (Homo sapiens): 9606, Mouse (Mus musculus): 10090
    metabolic_map (str {'metabolic', 'secondary', 'microbial', 'antibiotic'})
        Corresponds to setting map in ipath3. Select the overview map to use for the initial customization (default: `metabolic`)
    export_type (str {svg})
        Select the graphical file format for the generated map. Only SVG is available at the moment.

    """

    url = 'https://pathways.embl.de/mapping.cgi'

    post = {'selection': selection, 
        'default_opacity': default_opacity, 
        'default_edge_width': default_edge_width,
        'default_node_radius': default_node_radius,
        'keep_colors': keep_colors,
        'default_color': default_color,           
        'background_color': background_color, 
        'whole_pathways': whole_pathways,
        'whole_modules': whole_modules,
        'query_reactions': query_reactions,
        'tax_filter': tax_filter, 
        'metabolic_map': metabolic_map, 
        'export_type': export_type}
        
    ipath3_request = requests.post(url, json = post)

    # if ipath3_request.status_code == 400: 
    #     raise ValueError('Bad Request')

    return ipath3_request.content



def calculate_coverage(): 
    """ 
    Calculates the fraction of KEGG IDs corresponding to a specific molecular formula that were actually found in a metabolic map. 
    """
    pass 