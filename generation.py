# -*- coding:utf-8 -*-
# generation
import numpy as np
import scipy.sparse as sparse
import sys

from hyperg import HyperG

def gen_hg(X,with_feature,tad,hyperedge,hyperedge_flat):
    """
    :param X: numpy array, shape = (n_samples, n_features)
    :param with_feature: bool, optional(default=False)
    :return: instance of HyperG
    """
    # tad = np.loadtxt("/home/zsc/study/biye/output/chr19_50kb/optics_minsample_4.txt")
    # hyperedge = np.loadtxt('/home/zsc/study/biye/output/chr19_50kb/hyperedge_all.txt')
    # hyperedge_flat = np.loadtxt('/home/zsc/study/biye/output/chr19_50kb/hyperedge_flat.txt')
    n_nodes = tad.shape[0]
    n_edges = hyperedge.shape[0]

    node_idx = hyperedge_flat[:,0]
    edge_idx = hyperedge_flat[:,1]

    values = np.ones(node_idx.shape[0])

    H = sparse.coo_matrix((values, (node_idx, edge_idx)), shape=(n_nodes, n_edges))
    w = np.ones(n_edges)

    if with_feature:
        return HyperG(H, w = w, X=X)

    return HyperG(H,w = w)