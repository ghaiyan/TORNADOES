import numpy as np
import pyBigWig
import os
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import scipy.sparse as sparse
import sys
from hyperg import HyperG
step_length=50000 #resolution

# %load learning.py
from sklearn.cluster import k_means
from sklearn.cluster import SpectralClustering
from sklearn.utils import check_symmetric, check_random_state
from scipy.linalg import eigh

from hyperg import HyperG
from scipy.stats import pearsonr,spearmanr
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.metrics import jaccard_score, pairwise_distances
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.metrics.pairwise import manhattan_distances
from scipy.stats import pearsonr,spearmanr
def normalization(data):
    _range = np.max(data)-np.min(data)
    return (data-np.min(data)) / _range

def gen_hg(X,with_feature,tad,hyperedge,hyperedge_flat):
    """
    :param X: numpy array, shape = (n_samples, n_features)
    :param with_feature: bool, optional(default=False)
    :return: instance of HyperG
    """
    
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

def spectral_hg_partitioning(hg, n_clusters, assign_labels='kmeans', n_components=None, random_state=None, n_init=10):
    """
    :param hg: instance of HyperG
    :param n_clusters: int,
    :param assign_labels: str, {'kmeans', 'discretize'}, default: 'kmeans'
    :param n_components: int, number of eigen vectors to use for the spectral embedding
    :param random_state: int or None (default)
    :param n_init: int, number of time the k-means algorithm will be run
    with different centroid seeds.
    :return: numpy array, shape = (n_samples,), labels of each point
    """

    assert isinstance(hg, HyperG)
    assert n_clusters <= hg.num_nodes()

    random_state = check_random_state(random_state)

    if n_components is None:
        n_components = n_clusters

    L = hg.laplacian().toarray()
    L = check_symmetric(L)

    eigenval, eigenvec = eigh(L)
    embeddings = eigenvec[:, :n_components]
    embeddings = np.concatenate((embeddings, feature), axis=1)
    print(embeddings.shape)
    if assign_labels == 'kmeans':
        _, labels, _ = k_means(embeddings, n_clusters = n_clusters, random_state=random_state,
                               n_init=n_init)
    else:
        labels = SpectralClustering(n_clusters, random_state=random_state).fit(embeddings).labels_

    return labels
for chr in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]:
    cellline = 'IMR90'  
    #hic = np.load('/data/ghy_data/'+cellline+'/50kb/chr'+str(chr)+'_50000.npy')#Hi-c data
    #print(hic.shape)
    subcomparts = np.loadtxt('output/'+cellline+'/caspian_subcompart/'+cellline+'_chr'+str(chr)+'_50000_subcomparts.txt')
    tad = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')#此处用四列的TAD数据

    outpath = "output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge.txt"
    i = 0
    j = 0
    nodes = []

    with open(outpath, "w") as out:
        while i<len(subcomparts):
            start = int(subcomparts[i][0])
            end = int(subcomparts[i][2])
            while j<len(tad):
                if not(tad[j][2] + 3< start or tad[j][0] - end >3) :
                    nodes.append(j+1)
                j = j+1
            if len(nodes)>0:
                out.write("\t".join((str(nodes[0]),  str(nodes[-1]))) + "\n")
            nodes = []
            j = 0
            i = i+1
        out.close()
    hyperedge = np.loadtxt("output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge.txt")
    
    outpath1 = "output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge_all.txt"
    i = 0
    hyperedge2 = np.array([[0,0]])
    if hyperedge[0][0]>1:
        hyperedge2 = np.array([[1,int(hyperedge[0][0])-1]])

    while i < len(hyperedge)-1 :
        if int(hyperedge[i][1])+1 < int(hyperedge[i+1][0]):
            hyperedge2 = np.append(hyperedge2,np.array([[int(hyperedge[i][1])+1,int(hyperedge[i+1][0])-1]]),axis=0)
        i= i+1
    hyperedge2 = np.delete(hyperedge2, 0,0)
    hyperedge = np.append(hyperedge,hyperedge2,axis=0)
    hyperedge = np.sort(hyperedge,axis=0)
    np.savetxt(outpath1, np.c_[hyperedge],fmt='%d',delimiter='\t')

    hyperedge3 = np.loadtxt("output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge_all.txt")
    outpath2 = "output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge_flat.txt"

    res = []
    temp = []
    i = 0

    while i<len(hyperedge3):
        temp = list(range(int(hyperedge[i][0])-1,int(hyperedge[i][1])))
        for _ in temp:
            res.append([_, i])
        temp = []
        i = i+1

    res = np.array(res)
    print(res.shape)

    np.savetxt(outpath2, res, fmt='%d',delimiter='\t')
  
    A1_outpath = 'output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_A1.txt'
    B1_outpath = 'output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B1.txt'
    A2_outpath = 'output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_A2.txt'
    B2_outpath = 'output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B2.txt'
    B3_outpath = 'output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B3.txt'
    with open('output/'+cellline+'/SNIPER-sub-compartment/'+cellline+'_track_hg38.bed') as nkFile:
        lines = [line.strip().split() for line in nkFile]
    group = np.array(lines)
    start = group[:,1]
    end = group[:,2]
    AB_KR = group[:,3]
    with open(A2_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'A2' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'A2' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(A1_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'A1' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'A1' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(B1_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'B1' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'B1' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(B2_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'B2' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'B2' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1

    with open(B3_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'B3' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'B3' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1


    # 处理子区室结果

    A1_outpath = 'output/'+cellline+'/sub-compartment/chr'+str(chr)+'/A1.txt'
    B1_outpath = 'output/'+cellline+'/sub-compartment/chr'+str(chr)+'/B1.txt'
    A2_outpath = 'output/'+cellline+'/sub-compartment/chr'+str(chr)+'/A2.txt'
    B2_outpath = 'output/'+cellline+'/sub-compartment/chr'+str(chr)+'/B2.txt'
    with open('output/'+cellline+'/sub-compartment/chr'+str(chr)+'/chr'+str(chr)+'_sub_compartments.bed') as nkFile:
        lines = [line.strip().split() for line in nkFile]
    group = np.array(lines)
    start = group[:,1]
    end = group[:,2]
    AB_KR = group[:,3]
    with open(A2_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'A.2' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'A.2' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(A1_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'A.1' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'A.1' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(B1_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'B.1' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'B.1' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1
    with open(B2_outpath, "w") as out:
        left = -1
        for i in range(len(AB_KR)):
            if AB_KR[i][:3] == 'B.2' and left==-1:
                left = i
            elif left>=0 and (AB_KR[i][:3] != 'B.2' or str(AB_KR[i])=="nan") :
                out.write("\t".join((str(int((int(start[left])-1)/50000)+1),str(int(int(end[i-1])/50000)))) + "\n")
                left=-1

    #epi process, bigwig -> npy
    for input_epi in ["H3K4me3","POLR2A","H3K9me3","H3K27ac","CTCF"]:
        input_file="/data/ghy_data/"+cellline+"/"+input_epi+".bigWig"
        bw = pyBigWig.open(input_file)
        output_folder="epi/output/"+cellline+"/"
        chrom="chr"+str(chr)
        peaks = np.zeros(int(np.ceil(bw.chroms(chrom)/step_length)))
        print(peaks.shape)
        write_start, length, sum_pos = 0, 0, 0
        for entry in bw.intervals(chrom):
            start, end, pos = entry
            old_length = length
            length += end - start
            if length < step_length:
                sum_pos += (end - start) * pos
                continue
            while length + step_length>= 0:
                sum_pos += (step_length - old_length) * pos
                old_length = 0
                peaks[int(write_start/step_length)-1] = sum_pos / step_length
                write_start += step_length
                length -= step_length
                sum_pos = 0
            sum_pos += length * pos
        print(peaks[0])
        np.save(os.path.join(output_folder, '{}_{}b_{}'.format(chrom, step_length,input_epi)), peaks)

    data_CTCF = np.load('epi/output/'+cellline+'/chr'+str(chr)+'_50000b_CTCF.npy')
    data_H3K4me3 = np.load('epi/output/'+cellline+'/chr'+str(chr)+'_50000b_H3K4me3.npy')
    data_H3K27ac = np.load('epi/output/'+cellline+'/chr'+str(chr)+'_50000b_H3K27ac.npy')
    data_POLR2A = np.load('epi/output/'+cellline+'/chr'+str(chr)+'_50000b_POLR2A.npy')
    data_H3K9me3 = np.load('epi/output/'+cellline+'/chr'+str(chr)+'_50000b_H3K9me3.npy')

    print(data_CTCF.shape)
    print(data_H3K4me3.shape)
    print(data_H3K27ac.shape)
    print(data_POLR2A.shape)
    print(data_H3K9me3.shape)
    shape_num=min(data_CTCF.shape[0],data_H3K4me3.shape[0],data_H3K27ac.shape[0],data_POLR2A.shape[0],data_H3K9me3.shape[0])
    feature = np.dstack((data_CTCF[:shape_num], data_H3K4me3[:shape_num], data_H3K27ac[:shape_num], data_POLR2A[:shape_num], data_H3K9me3[:shape_num]))
    print(feature.shape)
    print(feature[0])
    # 特征归一化
    # feature_norm = []
    # for j in feature[0]:
    #     feature_norm.append(normalization(j))
    # np.savetxt("./output/feature_norm/chr"+str(i)+'.txt',
    #            feature_norm, fmt="%f", delimiter=' ')
    np.savetxt('epi/output/'+cellline+'/chr'+str(chr)+'.txt',feature[0], fmt="%f", delimiter=' ')
    ##计算每个TAD的特征值 ##TAD内节点特征值的均值
    def feature_edge(feature_path,out_path,chr):
        dataAll = []
        single_edge_feature = []
        feature_19 = np.loadtxt(feature_path)
        feature_shape=feature_19.shape[0]
        data_19 = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')

        f=open('output/'+cellline+'/50kb_TAD/chr'+str(chr)) #此处用二列的TAD数据
        for line in f.readlines():
            dataSet = []
            #for i in np.arange(int(line.split()[0]), int(line.split()[1])+1):
            for i in np.arange(int(line.split()[0]), feature_shape+1):
                dataSet.append(feature_19[i-1])
            single_TAD_feature = np.mean(dataSet, axis=0)
            dataAll.append(single_TAD_feature)
        f.close

        feature_norm = []
        for j in dataAll:
            feature_norm.append(normalization(j))

        #feature_norm = np.array(feature_norm, dtype=float)
        #feature_norm=feature_norm.astype(float)
        np.savetxt(out_path, feature_norm, fmt="%f")
        # np.savetxt("./output/feature_edge/chr19/"+'5861_PCA.txt',
        #            feature_PCA, fmt="%f", delimiter=' ')

    ##计算每个TAD的特征值 ##TAD内节点特征值的均值
    feature_edge('epi/output/'+cellline+'/chr'+str(chr)+'.txt','epi/output/'+cellline+'/chr'+str(chr)+'_TAD.txt',chr)
    feature = np.loadtxt('epi/output/'+cellline+'/chr'+str(chr)+'_TAD.txt')
    print(feature.shape)
    hyperedge_flat = np.loadtxt( "output/"+cellline+"/hypergraph/chr"+str(chr)+"_hyperedge_flat.txt")

    hg = gen_hg(None,False,tad,hyperedge3,hyperedge_flat)
    #change n_clusters to control the number of types
    for ncluster in [2,3,4,5]:
        clusters = spectral_hg_partitioning(hg, n_clusters = ncluster, assign_labels='kmeans', n_components=None, random_state=None, n_init=10)
        print(clusters)

        #use fanc tool to call A/B compartment
        A = np.loadtxt('output/'+cellline+'/AB/chr'+str(chr)+'_A.txt')
        B = np.loadtxt('output/'+cellline+'/AB/chr'+str(chr)+'_B.txt')
        # TAD = np.loadtxt('/home/zsc/study/biye/output/chr19_50kb/optics_new.txt')

        def getClass_AB():
            res = np.ones(len(tad))
            for i in range(len(tad)):
                start = tad[i][0]
                end = tad[i][1]
                for j in range(len(A)):
                    s = A[j][0]
                    e = A[j][1]
                    if not (start>e or end<s):
                        res[i] = 0
            return res
        cos_sim_AB = cosine_similarity(np.vstack((clusters, getClass_AB())))[0,1]
        print(clusters,getClass_AB())
        pccs = pearsonr(clusters, getClass_AB())
        print(cos_sim_AB)
        print(pccs)

        #A1,A2,B1,B2
        A1 = np.loadtxt('output/'+cellline+'/sub-compartment/chr'+str(chr)+'/A1.txt')
        B1 = np.loadtxt('output/'+cellline+'/sub-compartment/chr'+str(chr)+'/B1.txt')
        A2 = np.loadtxt('output/'+cellline+'/sub-compartment/chr'+str(chr)+'/A2.txt')
        B2 = np.loadtxt('output/'+cellline+'/sub-compartment/chr'+str(chr)+'/B2.txt')

        #A1,A2,B1,B2,B3
        A1_1 = np.loadtxt('output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_A1.txt')
        A2_1 = np.loadtxt('output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_A2.txt')
        B1_1 = np.loadtxt('output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B1.txt')
        B2_1 = np.loadtxt('output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B2.txt')
        B3_1 = np.loadtxt('output/'+cellline+'/SNIPER-sub-compartment/chr'+str(chr)+'_B3.txt')


        length=[0]*len(tad)
        def getlen(fir,sec):
            dist = 0
            start = fir[0]
            end = fir[1]
            s = sec[0]
            e = sec[1]
            if end>s and start<s:
                dist = end-s+1
            if end>e and start<e:
                dist = end-e+1
            return dist
        def getClass_sub():
            res = [-1]*len(tad)
            for i in range(len(tad)):
                start = tad[i][0]
                end = tad[i][1]

                for j in range(len(A1)):
                    s = A1[j][0]
                    e = A1[j][1]
                    if not(start>=e or end<=s):
                        res[i] = 0
                        if start>s and end<e:
                            length[i] = end-start+1
                        else:
                            length[i] = getlen(tad[i],A1[j])

                for j in range(len(A2)):
                    s = A2[j][0]
                    e = A2[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 1
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 1
                                length[i]=getlen(tad[i],A2[j])
                            else:
                                if getlen(tad[i],A2[j])>length[i]:
                                    length[i] = getlen(tad[i],A2[j])
                                    res[i] = 1
                for j in range(len(B1)):
                    s = B1[j][0]
                    e = B1[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 2
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 2
                                length[i]=getlen(tad[i],B1[j])
                            else:
                                if getlen(tad[i],B1[j])>length[i]:
                                    length[i] = getlen(tad[i],B1[j])
                                    res[i] = 2
                for j in range(len(B2)):
                    s = B2[j][0]
                    e = B2[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 3
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 3
                                length[i]=getlen(tad[i],B2[j])
                            else:
                                if getlen(tad[i],B2[j])>length[i]:
                                    length[i] = getlen(tad[i],B2[j])
                                    res[i] = 3
            return res
        
        def getClass_sub5():
            res = [-1]*len(tad)
            for i in range(len(tad)):
                start = tad[i][0]
                end = tad[i][1]

                for j in range(len(A1_1)):
                    s = A1_1[j][0]
                    e = A1_1[j][1]
                    if not(start>=e or end<=s):
                        res[i] = 0
                        if start>s and end<e:
                            length[i] = end-start+1
                        else:
                            length[i] = getlen(tad[i],A1_1[j])

                for j in range(len(A2_1)):
                    s = A2_1[j][0]
                    e = A2_1[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 1
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 1
                                length[i]=getlen(tad[i],A2_1[j])
                            else:
                                if getlen(tad[i],A2_1[j])>length[i]:
                                    length[i] = getlen(tad[i],A2_1[j])
                                    res[i] = 1
                for j in range(len(B1_1)):
                    s = B1_1[j][0]
                    e = B1_1[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 2
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 2
                                length[i]=getlen(tad[i],B1_1[j])
                            else:
                                if getlen(tad[i],B1_1[j])>length[i]:
                                    length[i] = getlen(tad[i],B1_1[j])
                                    res[i] = 2
                for j in range(len(B2_1)):
                    s = B2_1[j][0]
                    e = B2_1[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 3
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 3
                                length[i]=getlen(tad[i],B2_1[j])
                            else:
                                if getlen(tad[i],B2_1[j])>length[i]:
                                    length[i] = getlen(tad[i],B2_1[j])
                                    res[i] = 3
                for j in range(len(B3_1)):
                    s = B3_1[j][0]
                    e = B3_1[j][1]
                    if not(start>=e or end<=s):
                        if start>=s and end<=e:
                                res[i] = 3
                                length[i] = end-start+1
                        else:
                            if res[i] ==-1:
                                res[i] = 3
                                length[i]=getlen(tad[i],B3_1[j])
                            else:
                                if getlen(tad[i],B3_1[j])>length[i]:
                                    length[i] = getlen(tad[i],B3_1[j])
                                    res[i] = 3
            return res
        pccs = pearsonr(clusters, getClass_sub())
        cos_sim_subAB = cosine_similarity(np.vstack((clusters, getClass_sub())))[0,1]
        print("cos_sim of A1,A2,B1,B2",cos_sim_subAB)

        pccs = pearsonr(clusters, getClass_sub5())
        cos_sim_subAB5 = cosine_similarity(np.vstack((clusters, getClass_sub5())))[0,1]
        print("cos_sim of A1,A2,B1,B2,B3",cos_sim_subAB5)

        if ncluster==2:
            one_outpath = 'output/'+cellline+'/hypergraph/two/chr'+str(chr)+'_first.txt'
            two_outpath = 'output/'+cellline+'/hypergraph/two/chr'+str(chr)+'_second.txt'

            with open(one_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 0:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(two_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 1:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")

        if ncluster==3:
            one_outpath = 'output/'+cellline+'/hypergraph/three/chr'+str(chr)+'_first.txt'
            two_outpath = 'output/'+cellline+'/hypergraph/three/chr'+str(chr)+'_second.txt'
            three_outpath = 'output/'+cellline+'/hypergraph/three/chr'+str(chr)+'_three.txt'

            with open(one_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 0:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(two_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 1:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(three_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 2:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
        if ncluster==4:
            one_outpath = 'output/'+cellline+'/hypergraph/four/chr'+str(chr)+'_first.txt'
            two_outpath = 'output/'+cellline+'/hypergraph/four/chr'+str(chr)+'_second.txt'
            three_outpath = 'output/'+cellline+'/hypergraph/four/chr'+str(chr)+'_three.txt'
            four_outpath = 'output/'+cellline+'/hypergraph/four/chr'+str(chr)+'_four.txt'
            with open(one_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 0:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(two_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 1:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(three_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 2:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")
            with open(four_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 3:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")
        
        if ncluster==5:
            one_outpath = 'output/'+cellline+'/hypergraph/five/chr'+str(chr)+'_first.txt'
            two_outpath = 'output/'+cellline+'/hypergraph/five/chr'+str(chr)+'_second.txt'
            three_outpath = 'output/'+cellline+'/hypergraph/five/chr'+str(chr)+'_three.txt'
            four_outpath = 'output/'+cellline+'/hypergraph/five/chr'+str(chr)+'_four.txt'
            five_outpath = 'output/'+cellline+'/hypergraph/five/chr'+str(chr)+'_five.txt'
            with open(one_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 0:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(two_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 1:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3])))) + "\n")
            with open(three_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 2:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")
            with open(four_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 3:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")
            with open(five_outpath, "w") as out:
                for i in range(len(clusters)):
                    if clusters[i] == 4:
                        out.write("\t".join((str(int(tad[i,0])),str(int(tad[i][1])),str(int(tad[i][2])),str(int(tad[i][3]))))+ "\n")
        
        #similarity with A/B compartment
        first_class = []
        second_class = []
        third_class=[]
        fourth_class=[]
        fifth_class=[]
        def getCluster(i,j):
            for i in range(len(clusters)):
                if clusters[i]==i:
                    first_class.append(i)
                if clusters[i]==j:
                    second_class.append(i)
            return first_class,second_class
        res = getClass_AB()
        pccs,pvalue = pearsonr(clusters, res)
        def getScore(i,j):
            first_class,second_class = getCluster(i,j)
            res1 = list(res[first_class])+list(res[second_class])
            classes = list(clusters[first_class])+list(clusters[second_class])

            jaccard_sim = jaccard_score(classes, res1,average='micro')
            cos_sim = cosine_similarity(np.vstack((classes, res1)))[0,1]
            euclidean_dist = euclidean_distances(np.vstack((classes, res1)))[0,1]
            manhattan_dist = manhattan_distances(np.vstack((classes, res1)))[0,1]
            per,pvalue = pearsonr(classes, res1)
            return jaccard_sim,cos_sim,euclidean_dist,manhattan_dist,per

        a = getScore(0,1)
        b = getScore(0,2)
        c = getScore(1,2)
        jaccard_sim = max(a[0],b[0],c[0])
        cos_sim = max(a[1],b[1],c[1])
        euclidean_dist = min(a[2],b[2],c[2])
        manhattan_dist = min(a[3],b[3],c[3])
        per = c[4]
        print("Jaccard 相似系数: %0.3f"%(jaccard_sim))
        print("余弦相似度:", cos_sim)
        print("欧几里得距离:", euclidean_dist)
        print("曼哈顿距离:", manhattan_dist)
        print("pearson:",per)
        similarityofAB="output/"+cellline+"/similarityofAB.txt"
        
        
        with open(similarityofAB, "a+") as f:
            f.write("\t".join(("chr",str(chr),str(ncluster),str(jaccard_sim),str(cos_sim),str(euclidean_dist),str(manhattan_dist),str(per)))+'\n')
            f.close()
        
        #similarity with A1/A2/B1/B2 compartment
        print("similarity with A1/A2/B1/B2 compartment:")
        res = getClass_sub()
        pccs,pvalue = pearsonr(clusters, res)
        

        jaccard_sim = jaccard_score(clusters, res,average='micro')
        cos_sim = cosine_similarity(np.vstack((clusters, res)))[0,1]
        euclidean_dist = euclidean_distances(np.vstack((clusters, res)))[0,1]
        manhattan_dist = manhattan_distances(np.vstack((clusters, res)))[0,1]
        per,pvalue = pearsonr(clusters, res)
        print("Jaccard 相似系数: %0.3f"%(jaccard_sim))
        print("余弦相似度:", cos_sim)
        print("欧几里得距离:", euclidean_dist)
        print("曼哈顿距离:", manhattan_dist)
        print("pearson:",per)
        similarityofA1B1="output/"+cellline+"/similarityofA1B1.txt"
        with open(similarityofA1B1, "a+") as f:
            f.write("\t".join(("chr",str(chr),str(ncluster),str(jaccard_sim),str(cos_sim),str(euclidean_dist),str(manhattan_dist),str(per)))+'\n')
            f.close()

         #similarity with A1/A2/B1/B2/B3 compartment
        print("similarity with A1/A2/B1/B2/B3 compartment:")
        res5 = getClass_sub5()
        pccs,pvalue = pearsonr(clusters, res5)
        

        jaccard_sim5 = jaccard_score(clusters, res5,average='micro')
        cos_sim5 = cosine_similarity(np.vstack((clusters, res5)))[0,1]
        euclidean_dist5 = euclidean_distances(np.vstack((clusters, res5)))[0,1]
        manhattan_dist5 = manhattan_distances(np.vstack((clusters, res5)))[0,1]
        per5,pvalue = pearsonr(clusters, res5)
        print("Jaccard 相似系数: %0.3f"%(jaccard_sim5))
        print("余弦相似度:", cos_sim5)
        print("欧几里得距离:", euclidean_dist5)
        print("曼哈顿距离:", manhattan_dist5)
        print("pearson:",per5)
        similarityofA1B1_5="output/"+cellline+"/similarityofA1B1_5.txt"
        with open(similarityofA1B1_5, "a+") as f:
            f.write("\t".join(("chr",str(chr),str(ncluster),str(jaccard_sim5),str(cos_sim5),str(euclidean_dist5),str(manhattan_dist5),str(per5)))+'\n')
            f.close()