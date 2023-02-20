
from __future__ import division
from scipy import linalg
import scipy.sparse as sparse
from scipy.sparse import csr_matrix
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.metrics import roc_auc_score
from sklearn.metrics import confusion_matrix
from sklearn.ensemble import RandomForestClassifier
from scipy.optimize import nnls
import random 
#import randnet
import time
import timeit

from itertools import combinations

import numpy as np
import scipy.io as sio
from scipy.spatial.distance import cdist
from scipy.linalg import eigh
from scipy import sparse
from scipy.stats.mstats import spearmanr, kendalltau, pearsonr
from sklearn.metrics import auc, roc_curve, precision_recall_curve,average_precision_score
from sklearn.decomposition import TruncatedSVD
#import matplotlib.pyplot as plt



def NS(R, idx):
    n = R.shape[1]
    r = R.shape[0]
    rho = r/n
    A = np.zeros((n, n))
    for i, k in enumerate(idx):        
        A[k, :] = R[i]
        A[:, k] = R[i]
    A2 = R.T @ R * (1 / n)
    D = cdist(A2, A2, 'chebyshev') 
    K = D < np.percentile(D, np.sqrt(np.log(n) / n) * 100, 0)
    Q = A @ (K * (1 / (np.sum(K, 0) + 1e-10)))
    # return (Q + Q.T) * 0.5/rho
    return (Q + Q.T) * 0.5/(1-(1-rho)**2)



def global_CC(A_orig):
    A = A_orig
    np.fill_diagonal(A,0)
    A2 = np.matmul(A,A)
    A3 = np.matmul(A2,A)
    val_cc = np.trace(A3)/np.sum(A2)
    return val_cc

   
def genDist(n, d):
    X = np.random.normal(0, 1, (n, d))
    P = 1 / (1 + np.exp(cdist(X, X)))
    return P


def genProduct(n, d):
    X = np.random.beta(0.5, 1, (n, d)) 
    P = X @ X.T * (1 / d) 
    return P


def genSBM(n, d):
    S = np.ones((d, d)) * 0.3 / d
    np.fill_diagonal(S, np.arange(1, d) / (d + 1))
    B = np.equal.outer(np.arange(d), np.repeat(np.arange(d), n // d))
    P = B.T @ S @ B
    print("Density", np.mean(P))   
    return P
    

def AUC(Aout, Qout):
    fpr, tpr, _ = roc_curve(Aout, Qout, pos_label=1)
    return auc(fpr, tpr)
    

def evaluate(Aout, Pout, Qout, results, method_name, i, k, rho):
    if "pauc" in results:
        fpr, tpr, _ = roc_curve(Aout, Qout, pos_label=1)
        results["pauc"][method_name][i, k] = (auc(fpr, tpr) if np.isnan(results["pauc"][method_name][i, k]) else 
                                             max(results["pauc"][method_name][i, k], auc(fpr, tpr)))
    if "spea" in results or "ktau" in results:
        arr = list(zip(Qout, Pout))
        np.random.shuffle(arr)
        q, p = zip(*arr[:300])
        if "spea" in results:
            s = spearmanr(p, q)[0]        
            results["spea"][method_name][i, k] = 0 if np.isnan(s) else s
        if "ktau" in results:
            s = kendalltau(p, q)[0]
            results["ktau"][method_name][i, k] = 0 if np.isnan(s) else s               
    if "rmse" in results:
        results["rmse"][method_name][i, k] = np.linalg.norm(Pout - Qout) / np.linalg.norm(Pout)
    
    if "adj" in results:
        results["adj"][method_name][i, k] = np.linalg.norm(Pout - Qout / rho) / np.linalg.norm(Pout)
    print(method_name,
          [(measure, results[measure][method_name][i, k]) for measure in results])

def eigen_evaluate(A,A_impute):
    w_P, U_P = eigh(A, eigvals=(A.shape[0] - 1, A.shape[0] - 1))
    w_Q, U_Q = eigh(A_impute, eigvals=(A_impute.shape[0] - 1, A_impute.shape[0] - 1))
    c_star_P = np.max(np.abs(U_P))
    cent_P = np.sum(c_star_P-np.abs(U_P))
    c_star_Q = np.max(np.abs(U_Q))
    cent_Q = np.sum(c_star_Q-np.abs(U_Q))
    return np.abs(cent_P-cent_Q)/cent_P        
  

def evaluate_new(Aout, Pout, Pout_mat, Qout, Qout_mat, results, method_name, i, k, rho,A_impute,A):
    if "pauc" in results:
        print("check pauc")
        fpr, tpr, _ = roc_curve(Aout.A1, Qout.A1, pos_label=1)
        results["pauc"][method_name][i, k] = (auc(fpr, tpr) if np.isnan(results["pauc"][method_name][i, k]) else 
                                             max(results["pauc"][method_name][i, k], auc(fpr, tpr)))
    if "spea" in results or "ktau" in results:
        print("check spea or ktau")
        arr = list(zip(Qout, Pout))
        np.random.shuffle(arr)
        q, p = zip(*arr[:300])
        if "spea" in results:
            print("check spea")
            s = spearmanr(p, q)[0]        
            results["spea"][method_name][i, k] = 0 if np.isnan(s) else s
        if "ktau" in results:
            print("check ktau")
            s = kendalltau(p, q,method='asymptotic')[0]
            #s = kendalltau(p, q)[0]
            results["ktau"][method_name][i, k] = 0 if np.isnan(s) else s               
    if "rmse" in results:
        results["rmse"][method_name][i, k] = np.linalg.norm(Pout - Qout) / np.linalg.norm(Pout)    
    if "density" in results:
        results["density"][method_name][i, k] = np.abs(np.mean(A_impute)-np.mean(A))/np.mean(A)
    if "eigen_centrality" in results:
        w_P, U_P = eigh(A, eigvals=(A.shape[0] - 1, A.shape[0] - 1))
        w_Q, U_Q = eigh(A_impute, eigvals=(A_impute.shape[0] - 1, A_impute.shape[0] - 1))
        c_star_P = np.max(np.abs(U_P))
        cent_P = np.sum(c_star_P-np.abs(U_P))
        c_star_Q = np.max(np.abs(U_Q))
        cent_Q = np.sum(c_star_Q-np.abs(U_Q))
        results["eigen_centrality"][method_name][i, k] = np.abs(cent_P-cent_Q)/cent_P        
    print(method_name,
          [(measure, results[measure][method_name][i, k]) for measure in results])

    
def sparseReader(filename):
    with open(filename, 'r') as data:
        data.readline()
        _, m, n, *_ = data.readline().split()
        n = int(n)
        m = int(m)
        print("#nodes = {0}, #edges = {1}".format(n, m))
        A = sparse.dok_matrix((n, n), dtype=np.int)
        for line in data:
            if line[0] != '%':
                s, t, *_ = line.split()
                i = int(s) - 1
                j = int(t) - 1
                A[i, j] = 1
                A[j, i] = 1
    print("Data loaded")
    return n, A


def IIDSampler(A, rho):
    n = A.shape[0]
    sample = []
    out_set = []
    rho_prime = 2 * rho - rho ** 2
    Omega = np.zeros(A.shape)
    for i in range(n - 1):
        for j in range(i + 1, n):
            if np.random.uniform(0, 1) < rho_prime:
                Omega[i, j] += 1
                Omega[j, i] += 1
                sample.append((i, j))
            else:
                out_set.append((i, j))                
    return Omega, sample, list(zip(*out_set))


def egoSampler(A, rho):
    n = A.shape[0]
    nsample = int(n * rho)
    idx = np.random.choice(range(n), nsample, replace=False)
    out_set = np.array([i for i in np.arange(n) if i not in set(idx)])
    return A[idx, :], idx, out_set
     

def getOutSet(X, out_ego, out_iid):
    return X[out_ego][:, out_ego].ravel(), X[out_iid[0], out_iid[1]].ravel()


def generator(P, deg_avg):
    n = P.shape[0]
    Q = P * (n  * deg_avg / np.sum(P))
    A = np.tril(np.random.uniform(size=(n, n)) < Q, -1) * 1
    A = A + A.T    
    return A, Q


def stack_impute(A,P_hat_mat,out_ego):
    n = A.shape[0]
    in_ego = np.setxor1d(range(n),out_ego)
    R11 = A[in_ego][:,in_ego]
    R12 = A[in_ego][:,out_ego]
    R21 = A[out_ego][:,in_ego]
    R22 = P_hat_mat
    R1x = np.concatenate((R11,R12),1)
    R2x = np.concatenate((R21,R22),1)
    A_impute = np.concatenate((R1x,R2x),0)
    return A_impute




def adj_to_nodes_edges(A):
    
    """ 
    This function change adjacency matrix to list of nodes and edges.

    Input and Parameters:
    -------
    A: the adjacency matrix

    Returns:
    -------
    nodes: node list of the given network
    edges: edge list of the given network

    Examples:
    -------
    >>> nodes, edges = adj_to_nodes_edges(A)
    """
    
    num_nodes = A.shape[0]
    nodes = range(num_nodes)
    edges = np.where(np.triu(A,1))
    row = edges[0]
    col = edges[1]
    edges = np.vstack((row,col)).T
    return nodes, edges


def auc_evaluate(Aout,Qout):
    print("check pauc")
    #fpr, tpr, _ = roc_curve(Aout.A1, Qout.A1, pos_label=1)
    precision, recall, _ = precision_recall_curve(Aout.A1, Qout.A1, pos_label=1)
    #plt.plot(precision, recall, 'o', color='black')
    #pAUC = auc(fpr, tpr) 
    pAUC = auc(recall, precision) 
    print("AUC: " +str(np.round(pAUC,2)))
    return pAUC, precision, recall


def auc_evaluate_orig(Aout,Qout):
    print("check pauc")
    fpr, tpr, _ = roc_curve(Aout.A1, Qout.A1, pos_label=1)
    #precision, recall, _ = precision_recall_curve(Aout.A1, Qout.A1, pos_label=1)
    #plt.plot(precision, recall, 'o', color='black')
    pAUC = auc(fpr, tpr) 
    #pAUC = auc(recall, precision) 
    print("AUC: " +str(np.round(pAUC,2)))
    return pAUC, fpr, tpr

   
        
    

""" 
    This function converts an edgelist to the adjacency matrix
""" 
def edgelist2Adj(edges_orig):
    num_nodes = int(np.max(edges_orig)) + 1
    row = np.array(edges_orig)[:,0]
    col = np.array(edges_orig)[:,1]

    data_aux = np.ones(len(row))
    A_orig = csr_matrix((data_aux,(row,col)),shape=(num_nodes,num_nodes))
    A_orig = sparse.triu(A_orig,1) + sparse.triu(A_orig,1).transpose()
    A_orig[A_orig>0] = 1 
    A_orig = A_orig.todense()
    return A_orig




def ego_Adj_NS_evaluation(A_orig,Q,rho=0.5): 
    
    """ 
    This function evaluates Neighborhood Smoothing (NS) performance on a network (A_orig)

    Input and Parameters:
    -------
    A_orig: the adjacency matrix of the original network
    rho: ego-sample proportion
    Returns:
    -------

     """
    
    #### load the original netowrk A_orig
    A_orig = np.asmatrix(A_orig)
    #### construct the holdout and training matriced from the original matrix
    
    P = A_orig
    Omega, sample_iid, out_iid = IIDSampler(A_orig, rho)
    R, sample_ego, out_ego = egoSampler(A_orig, rho)
    Pout_ego, Pout_iid = getOutSet(P, out_ego, out_iid)
    Aout_ego, Aout_iid = getOutSet(A_orig, out_ego, out_iid)
    
    
    print("Run NS")
    start_time = timeit.default_timer()
    Phat = np.asmatrix(NS(R, sample_ego))
    runtime = timeit.default_timer() - start_time
    P_hat_mat = Phat[out_ego][:, out_ego]
    P_hat = P_hat_mat.ravel()
    A_impute = stack_impute(A_orig,P_hat_mat,out_ego)
        
    # PR_AUC, precision, recall = auc_evaluate(Aout_ego, P_hat)
    auc, fpr, tpr = auc_evaluate_orig(Aout_ego, P_hat)

    #pred_cc = global_CC(A_impute)    
    #true_cc = global_CC(A_orig)    
    #cc_NS = np.abs(pred_cc-true_cc)/true_cc
 
    #density_NS = np.abs(np.mean(A_impute)-np.mean(A_orig))/np.mean(A_orig)
    #eigen_NS = eigen_evaluate(A_orig,A_impute)
    
    norm_NS = np.linalg.norm(Q[out_ego][:, out_ego].ravel()-P_hat,ord='fro')
    rel_err_NS = norm_NS/np.linalg.norm(Q[out_ego][:, out_ego],ord='fro')
    rel_err_full_NS = np.linalg.norm(Q-Phat,ord='fro')/np.linalg.norm(Q,ord='fro')
    mse = np.linalg.norm(Q[out_ego][:, out_ego].ravel()-P_hat,ord='fro')**2/len(out_ego)**2
    mse_full = np.linalg.norm(Q-Phat,ord='fro')**2/A_orig.shape[0]**2

    #return auc, fpr, tpr,density_NS, eigen_NS, PR_AUC, precision, recall, cc_NS, norm_NS, rel_err_NS
    return auc, fpr, tpr, norm_NS, rel_err_NS, rel_err_full_NS, mse, mse_full, runtime





def ego_LP_NS_evaluation(edges_orig,rho=0.5): 
    
    """ 
    This function evaluates Neighborhood Smoothing (NS) performance on a network (edge_orig) 

    Input and Parameters:
    -------
    A_orig: the adjacency matrix of the original network
    rho: ego-sample proportion
    Returns:
    -------

     """
    
    #### load the original netowrk A_orig
    A_orig = edgelist2Adj(edges_orig)
    #### construct the holdout and training matriced from the original matrix
    
    P = A_orig
    Omega, sample_iid, out_iid = IIDSampler(A_orig, rho)
    R, sample_ego, out_ego = egoSampler(A_orig, rho)
    Pout_ego, Pout_iid = getOutSet(P, out_ego, out_iid)
    Aout_ego, Aout_iid = getOutSet(A_orig, out_ego, out_iid)
    
    
    print("Run NS")
    Phat = np.asmatrix(NS(R, sample_ego))
    P_hat_mat = Phat[out_ego][:, out_ego]
    P_hat = P_hat_mat.ravel()
    A_impute = stack_impute(A_orig,P_hat_mat,out_ego)
        
    PR_AUC, precision, recall = auc_evaluate(Aout_ego, P_hat)
    auc, fpr, tpr = auc_evaluate_orig(Aout_ego, P_hat)

    pred_cc = global_CC(A_impute)    
    true_cc = global_CC(A_orig)    
    cc_NS = np.abs(pred_cc-true_cc)/true_cc
 
    density_NS = np.abs(np.mean(A_impute)-np.mean(A_orig))/np.mean(A_orig)
    eigen_NS = eigen_evaluate(A_orig,A_impute)
    


    return auc, fpr, tpr,density_NS, eigen_NS, PR_AUC, precision, recall, cc_NS


def nonuniform_generator(P, deg_avg):
    n = P.shape[0]
    P_ind = np.argsort(P.sum(axis=0))[::-1]
    P = P[P_ind, :]
    P = P[:, P_ind]
    Q = P * (n  * deg_avg / np.sum(P))
    A = np.tril(np.random.uniform(size=(n, n)) < Q, -1) * 1
    A = A + A.T    
    return A, Q


def nonuniformSampler(A, rho, n_ratio, sample_ratio):
    n = A.shape[0]
    n_part = np.round(n*n_ratio).astype(int)
    if sum(n_part)!=n:
        n_part[1] += n-sum(n_part)
    sample_part = np.minimum(np.ceil(n_part*rho*sample_ratio).astype(int),n_part) # Just insurance
    idx = np.array([])
    temp = 0
    for i in range(len(n_ratio)):
        idx = np.append(idx,temp+np.random.choice(range(n_part[i]), sample_part[i], replace=False))
        temp += n_part[i]
    idx = idx.astype(int)
    nsample = int(sum(sample_part))
    # idx = np.random.choice(range(n), nsample, replace=False)
    out_set = np.array([i for i in np.arange(n) if i not in set(idx)])
    return A[idx, :], idx, out_set


def ego_nonuniform_Adj_NS_evaluation(A_orig,Q,rho=0.5,n_ratio=np.array([.33,.34,.33]),sample_ratio=np.array([1.5,1,.5])): 
    
    """ 
    This function evaluates Neighborhood Smoothing (NS) performance on a network (A_orig)

    Input and Parameters:
    -------
    A_orig: the adjacency matrix of the original network
    rho: ego-sample proportion
    Returns:
    -------

     """
    
    #### load the original netowrk A_orig
    A_orig = np.asmatrix(A_orig)
    #### construct the holdout and training matriced from the original matrix
    
    P = A_orig
    Omega, sample_iid, out_iid = IIDSampler(A_orig, rho)
    # R, sample_ego, out_ego = egoSampler(A_orig, rho)
    R, sample_ego, out_ego = nonuniformSampler(A_orig, rho, n_ratio, sample_ratio)
    Pout_ego, Pout_iid = getOutSet(P, out_ego, out_iid)
    Aout_ego, Aout_iid = getOutSet(A_orig, out_ego, out_iid)
    
    
    print("Run NS")
    start_time = timeit.default_timer()
    Phat = np.asmatrix(NS(R, sample_ego))
    runtime = timeit.default_timer() - start_time
    P_hat_mat = Phat[out_ego][:, out_ego]
    P_hat = P_hat_mat.ravel()
    A_impute = stack_impute(A_orig,P_hat_mat,out_ego)
        
    # PR_AUC, precision, recall = auc_evaluate(Aout_ego, P_hat)
    auc, fpr, tpr = auc_evaluate_orig(Aout_ego, P_hat)

    #pred_cc = global_CC(A_impute)    
    #true_cc = global_CC(A_orig)    
    #cc_NS = np.abs(pred_cc-true_cc)/true_cc
 
    #density_NS = np.abs(np.mean(A_impute)-np.mean(A_orig))/np.mean(A_orig)
    #eigen_NS = eigen_evaluate(A_orig,A_impute)
    
    norm_NS = np.linalg.norm(Q[out_ego][:, out_ego].ravel()-P_hat,ord='fro')
    rel_err_NS = norm_NS/np.linalg.norm(Q[out_ego][:, out_ego],ord='fro')
    rel_err_full_NS = np.linalg.norm(Q-Phat,ord='fro')/np.linalg.norm(Q,ord='fro')
    mse = np.linalg.norm(Q[out_ego][:, out_ego].ravel()-P_hat,ord='fro')**2/len(out_ego)**2
    mse_full = np.linalg.norm(Q-Phat,ord='fro')**2/A_orig.shape[0]**2

    #return auc, fpr, tpr,density_NS, eigen_NS, PR_AUC, precision, recall, cc_NS, norm_NS, rel_err_NS
    return auc, fpr, tpr, norm_NS, rel_err_NS, rel_err_full_NS, mse, mse_full, runtime


def nonuniformbyA_generator(P, deg_avg):
    n = P.shape[0]
    Q = P * (n  * deg_avg / np.sum(P))
    A = np.tril(np.random.uniform(size=(n, n)) < Q, -1) * 1
    A = A + A.T
    A_ind = np.argsort(A.sum(axis=0))[::-1]
    Q = Q[A_ind, :]
    Q = Q[:, A_ind]
    A = A[A_ind, :]
    A = A[:, A_ind]
    return A, Q