
### environment

import logging
import warnings

import os
import pandas as pd
import numpy as np
import sklearn
# from sklearn.metrics import balanced_accuracy_score
from sklearn import linear_model

### global settings

MAX_ITER = 1000
# FDR = 0.1
FDR = 1 # ONLY FOR TESTING

### functions

# TODO: do also ridge or lasso
def da_logit_from_celltype(in_table, C):
    """Description.

    Parameters:
    ------------
    in_table 
        file containing the table with milo output + soft cell-type annotation per neighbourhood
    C
        inverse of the regularisation weight 
    """
    M = pd.read_csv(in_table, sep = "\t", header = 0, index_col = 0)
    X = M.iloc[:,7:len(M.columns)] # remove the output of milo
    nhoods = len(X.index)
    ncelltypes = len(X.columns)

    b_train = []
    for x in range(nhoods):
        if M["logFC"].iloc[x] > 0 and M["SpatialFDR"].iloc[x] < FDR:
            b_train.append(1)
        elif M["logFC"].iloc[x] < 0 and M["SpatialFDR"].iloc[x] < FDR:
            b_train.append(-1)
        else:
            b_train.append(0)
    y = np.array(b_train).astype(int) # bool: expanded (1) or not (0)

    # logistic regression
    pipe_lr = linear_model.LogisticRegression(max_iter = MAX_ITER, C = C, multi_class = "multinomial").fit(X, y)
    coef = pipe_lr.coef_
    names = X.columns
    
    return coef, names


def de_logit_from_celltype(in_table, C):
    """Description.

    Parameters:
    ------------
    in_table 
        file containing the table with miloDE output + soft cell-type annotation per neighbourhood
    C
        inverse of the regularisation weight 
    """
    M = pd.read_csv(in_table, sep = "\t", header = 0, index_col = 0)
    nhood_names = np.unique(M['Nhood'])
    genes = np.unique(M['gene'])
    XX = M.iloc[:,8:len(M.columns)] # remove the output of miloDE
    ncelltypes = len(XX.columns)
    
    b_train = []
    for x in range(len(M.index)):
        if M["logFC"].iloc[x] > 0 and M["pval_corrected_across_genes"].iloc[x] < FDR and M["pval_corrected_across_nhoods"].iloc[x] < FDR:
            b_train.append(1)
        elif M["logFC"].iloc[x] < 0 and M["pval_corrected_across_genes"].iloc[x] < FDR and M["pval_corrected_across_nhoods"].iloc[x] < FDR:
            b_train.append(-1)
        else:
            b_train.append(0)
    yy = np.array(b_train).astype(int) # bool: expanded (1) or not (0)

    # logistic regression for each gene
    names = XX.columns
    XX = XX.to_numpy()
    coef_matrix = np.full([len(genes),ncelltypes], np.nan)
    for i in range(len(genes)):
        is_gene = (M["gene"].to_numpy() == genes[i])
        is_not_nan = (np.logical_not(np.isnan(M["logFC"].to_numpy())))
        idx = np.where(is_gene & is_not_nan)
        X = XX[idx]
        y = yy[idx]
        if len(np.unique(y)) > 1:
            pipe_lr = linear_model.LogisticRegression(max_iter = MAX_ITER, C = C, multi_class = "multinomial").fit(X, y)
            coef_matrix[i,:] = pipe_lr.coef_
    coef_df = pd.DataFrame(coef_matrix, index = genes)
    coef_df.columns = names
    
    return coef_df


def write_da_logit_from_celltype(dirname, coef, names, C):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    filename = os.path.join(dirname, "logit_" + str(C) + ".tsv")
    df = pd.DataFrame(coef)
    df.columns = names
    df.to_csv(filename, sep = "\t")


def write_de_logit_from_celltype(dirname, coef_df, C):
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    filename = os.path.join(dirname, "logit_" + str(C) + ".tsv")
    coef_df.to_csv(filename, sep = "\t")


