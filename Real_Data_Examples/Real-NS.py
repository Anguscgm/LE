#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 12:09:09 2022

@author: tianxili
"""

import sys

import pickle  
import NSHeaderReal as NS
import numpy as np
import scipy.io
mat = scipy.io.loadmat('data/'+sys.argv[1]+'-A.mat')
# P = mat['P']
A_orig = mat['A']
method = "NS"
rho_seq = [.5,.3,.2,.1]

M = 50

for rho in rho_seq:
    result = {}
    auc_array =[]
    pr_auc_array=[]
    tpr_array = []
    fpr_array = []
    precision_array= []
    recall_array = []
    runtime_array = []
    #density_array =[]
    #eigen_array = []
    #cc_array = []
    # norm_array = []
    # rel_err_array = []
    # rel_err_full_array = []
    
    for j in range(0,M):
        # A_orig, Q = NS.generator(P, deg)    
        auc_NS, fpr_NS, tpr_NS, PR_AUC, precision, recall, runtime = NS.ego_Adj_NS_evaluation(A_orig,rho)
        auc_array.append(auc_NS)
        pr_auc_array.append(PR_AUC)
        tpr_array.append(tpr_NS)
        fpr_array.append(fpr_NS)
        precision_array.append(precision)
        recall_array.append(recall)
        runtime_array.append(runtime)
        #density_array.append(density_NS)
        #eigen_array.append(eigen_NS)
        #cc_array.append(cc_NS)
        # norm_array.append(norm_NS)
        # rel_err_array.append(rel_err_NS)
        # rel_err_full_array.append(rel_err_full_NS)
    
    #result = {'auc_array':auc_array,'pr_auc_array':pr_auc_array, 'tpr_array': tpr_array, 'fpr_array': fpr_array,'precision_array': precision_array, 'recall_array': recall_array,'density_array':density_array,'eigen_array':eigen_array,'cc_array':cc_array, 'norm_array':norm_array, 'rel_err_array':rel_err_array}
    result = {'auc_array':auc_array,'pr_auc_array':pr_auc_array, 'tpr_array': tpr_array, 'fpr_array': fpr_array,
    'precision_array': precision_array, 'recall_array': recall_array,'runtime_array': runtime_array}
     #print(result)
    np.save(method+"_Results/"+sys.argv[1]+"-"+method+"-rho="+str(rho)+".npy",result)
                







