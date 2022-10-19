import pandas as pd
import numpy as np 
import re
from typing import List

def convert_to_rgb(val, colors, maxval = 1, minval = 0): 
    """
    Creates a color gradient for values that are between a maximum and minimum value
    https://stackoverflow.com/questions/20792445/calculate-rgb-value-for-a-range-of-values-to-create-heat-map 

    """
    EPSILON = 0.001  # Smallest difference.

    # Relative position of value in range minval-maxval 
    
    i_f = (val-minval)/(maxval-minval) * (len(colors)-1)
    
    if not np.isnan(i_f): 
        i, f = int(i_f // 1), i_f % 1
    else: return [100, 100, 100]

    if f < EPSILON: 
        return colors[i]     
    else: 
        # Interpolate between colors
        (r1,g1,b1), (r2,g2,b2) = colors[i], colors[i+1]
        r, g, b = r1 + f*(r2-r1), g1+f*(g2-g1), b1+f*(b2-b1)
        return [int(r), int(g), int(b)]



def color(df:pd.DataFrame, column:str, colors:List[tuple], normalization = None, colortype = 'HEX') -> list: 
    """ 
    Returns list of colors on a linear color scale of argument `colors` according to values in `column` of the dataframe `df`. 
    Different normalization methods are applicable. Colors can be returned 

    INPUT
    ----- 
    df (`class`:pd.DataFrame)
        Pandas DataFrame that contains a column whose values should be colorcoded
    column (str)
        Name of column 
    colors (`List[tuple]`)
        List of tuples of colors in RGB format (no HEX support yet)
    normalization (str, {zscore, log, quantilXX})
        Normalization procedure 
        zscore: Normalizes according to z-score of values 
        .. math::
            Z = (X_i - < X >) / \sigma(X)

        log: Normalizes according to log of values
        quantilXX: Normalizes according to quantils of values were XX quantils are considered. 
        E.g. quantil2 classifies values into 0 (50%) and 1 (50%)
    """
    values = df[column]


    if normalization == 'zscore': 
        values = (values - values.mean())/values.std()
    
    if normalization == 'log': 
        values = np.log(values)

    if re.match('quantil\d*', str(normalization)): 
        nrQuantils = int(re.findall('quantil(\d*)', normalization)[0])
        quantils = [numerator/nrQuantils for numerator in range(nrQuantils)]
        perc = np.quantile(values, quantils)
        values = pd.Series(np.digitize(values, perc))


    minval, maxval = values.min(), values.max()
    color_list_rgb = values.apply(lambda x: convert_to_rgb(x, colors, maxval, minval))

    if colortype == 'RGB': 
        return [f'RGB({rgb[0]},{rgb[1]},{rgb[2]})' for rgb in color_list_rgb] 
    
    else: return [f'#{rgb[0]:02x}{rgb[1]:02x}{rgb[2]:02x}' for rgb in color_list_rgb]


def size(df:pd.DataFrame, column:str, minsize:int = 2, maxsize:int = 20, normalization = None) -> list: 
    """ 
    Returns list of sizes on a linear color scale from size_min to size_max

    INPUT
    ----- 
    df (`class`:pd.DataFrame)
        Pandas DataFrame that contains a column whose values should be colorcoded
    column (str)
        Name of column 
    size_min (int)
        Minimal size in px
    size_max (int)
        Maximal size in px
    normalization (str, {zscore, log, quantilXX})
        Normalization procedure 
        zscore: Normalizes according to z-score of values 
        .. math::
            Z = (X_i - < X >) / \sigma(X)
        log: Normalizes according to log of values
        quantilXX: Normalizes according to quantils of values were XX quantils are considered. 
        E.g. quantil2 classifies values into 0 (50%) and 1 (50%)
    """
    import re 
    values = df[column]


    if normalization == 'zscore': 
        values = (values - values.mean())/values.std()
    
    if normalization == 'log': 
        values = np.log(values)

    if re.match('quantil\d*', str(normalization)): 
        nrQuantils = int(re.findall('quantil(\d*)', normalization)[0])
        quantils = [numerator/nrQuantils for numerator in range(nrQuantils)]
        perc = np.quantile(values, quantils)
        values = pd.Series(np.digitize(values, perc))


    minval, maxval = values.min(), values.max()

    fractions = [(val - minval)/(maxval - minval) for val in values]
    sizes = [int(fraction*(maxsize-minsize) + minsize) for fraction in fractions]
    
    return sizes 
    

    