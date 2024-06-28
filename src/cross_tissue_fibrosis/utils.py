### environment

import pandas as pd
import numpy as np

### functions

def pandas_rename_columns(df, keys, values):
    """
    Rename the columns and add empty columns if keys are not found
    """
    # make sure that there is no self-renaming
    i = 0
    while i < len(keys):
        if i < len(keys) and keys[i] == values[i]:
            keys = np.delete(keys, i)
            values = np.delete(values, i)
        else:
            i = i + 1
    # remove duplicated names
    df = df.drop(list(set(values) & set(df.columns)), axis=1)
    # add missing keys if they do not exist yet
    for i in range(len(keys)):
        if not keys[i] in df.columns:
            df[keys[i]] = None
    # do the actual mapping
    mapping = {keys[i]: values[i] for i in range(len(keys))}
    df.rename(columns = mapping, inplace = True) # Labels not contained in a dict / Series will be left as-is
    return df


def R_match_idx(v, w):
    """
    Mimics the R match() function on numpy arrays
    """
    match = [None] * len(v)
    for i in range(len(v)):
        if v[i] in w:
            match[i] = np.where(w == v[i])[0]
    return np.concatenate(match)


def R_match(w, index):
    """
    Update the indexes after R match
    """
    mapping = w
    mapping.index = index
    return mapping

# TODO: add SYMBOL to ENSEMBL mapping

