#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 10:20:30 2021

@author: qiaojixu
"""

import numpy as np
# from matplotlib import pyplot as plt
import pandas as pd
# import os
import scipy.cluster.hierarchy as sch
import itertools
# import sys



# def plot_corr(df,size):
#     '''Plot a graphical correlation matrix for a dataframe.

#     Input:
#         df: pandas DataFrame
#         size: vertical and horizontal size of the plot'''
    
#     # %matplotlib inline
#     # import matplotlib.pyplot as plt

#     # Compute the correlation matrix for the received dataframe
#     # corr = df.corr()
#     corr = df
#     # Plot the correlation matrix
#     fig, ax = plt.subplots(figsize=(size, size))
#     # cax = ax.matshow(corr, cmap='Greys')
#     plt.xticks(range(len(corr.columns)), corr.columns, rotation=90);
#     plt.yticks(range(len(corr.columns)), corr.columns);
    
#     # Add the colorbar legend
#     # cbar = fig.colorbar(cax, ticks=[-1, 0, 1], aspect=40, shrink=.8)

def get_GroupID(matrixCo, numC):

    matrixCo = pd.DataFrame(matrixCo)
    # matrixCo =sys.argv[1]
    
    matrixCo.index =  matrixCo.columns
    
    X = matrixCo.values
    d = sch.distance.pdist(X)
    L = sch.linkage(d, method='complete')
    # ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    ind = sch.fcluster(L, numC, 'maxclust')
    # print(ind)
    columns = [matrixCo.columns.tolist()[i] for i in list(np.argsort(ind))]
    matrixCo = matrixCo.reindex(columns, axis=1).reindex(columns, axis=0)
    
    unique, counts = np.unique(ind, return_counts=True)
    counts = dict(zip(unique, counts))
    print(counts)
    
    
    # plot_corr(matrixCo, 18)
    
    ##create cluster label for each contig
    groupID = []
    for item in unique:
        # print(item, counts.get(item))
        reps = itertools.repeat(item, counts.get(item))
        groupID += reps
    
    
    ##output of our data
    outidx = pd.DataFrame(columns)
    outidx.columns = ["X"]
    outidx["cluster id"] = groupID
    return outidx
    # outidx.to_excel("newcluster_trn1.xlsx", index=False, header=True)



