#!/usr/bin/env python

import scipy.stats
import pandas as pd

def inv_ranknorm(unnorm_table, meth = 'first'):
    '''
    similar to the r function used here: https://github.com/cartercompbio/TIMEgermline/blob/master/scripts/gwas/rank-norm.sh
    
    Inputs:
        unnorm_table: pandas dataframe with numerical values for all columns
        meth: option (default 'first') method for handling ties in rank method; 'average', 'first', 'last', 'dense'
    Outputs:
        norm_table: inverse rank normalized values with same column and row index as input unnorm_table
        
    test showed the biggest deviation from R version was <10E-14
    test showed "method='first'" yielded distributions centered on zero, even with many ties (e.g. sparse with multiple zero valued samples)
    '''
    norm_table = scipy.stats.norm.ppf((unnorm_table.rank(method = meth)-0.5) / (~unnorm_table.isna()).sum())
    norm_table = pd.DataFrame(norm_table, index = unnorm_table.index, columns = unnorm_table.columns)
    return norm_table

def inv_ranknorm_ser(unnorm_ser, meth = 'first'):
    norm_table = scipy.stats.norm.ppf((unnorm_ser.rank(method = meth)-0.5) / (~unnorm_ser.isna()).sum())
    norm_table = pd.Series(norm_table, index = unnorm_ser.index)
    return norm_table    
