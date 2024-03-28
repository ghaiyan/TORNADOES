import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import heapq
import argparse
import math
import matplotlib as mpl
import hdbscan
import datetime
import sklearn.metrics as metrics
mpl.rcParams['pdf.fonttype']=42
mpl.rcParams['ps.fonttype']=42

cellline="IMR90"

#22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1
#gm12878:22,20,21,19,15,14,13,12,11,10,9,8,7,6,4,3,2,1
#K562:1,2,3,4,5,6,8,10,11,12,14,15,16,17,18,19,22
#IMR90:1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22
for chr in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,22]:
    for clusternum in ["two","three","four","five"]:
        if clusternum=="two":
            for type in ["first","second"]:
                feature_path='epi/output/'+cellline+'/chr'+str(chr)+'.txt'
                dataAll = []
                single_edge_feature = []
                feature_19 = np.loadtxt(feature_path)
                #data_19 = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')
                f=open("output/"+cellline+"/hypergraph/"+clusternum+"/chr"+str(chr)+"_"+type+".txt")
                for line in f.readlines():
                    start = int(line.split()[0])
                    end = int(line.split()[2])
                    # if start>=48 and end<=196:
                    dataSet = []
                    for i in np.arange(start, end+1):
                        dataSet.append(feature_19[i-1])
                    single_TAD_feature = np.sum(dataSet, axis=0)
                    dataAll.append(single_TAD_feature)
                f.close
                v_H3K27ac = np.sum(dataAll, axis=0)[0]
                v_H3K4me3 = np.sum(dataAll, axis=0)[1]
                v_CTCF = np.sum(dataAll, axis=0)[2]
                v_POLR2A = np.sum(dataAll, axis=0)[3]
                v_H3K9me3 = np.sum(dataAll, axis=0)[4]
                densityFile="output/"+cellline+"/factors-density-metric-5.txt"
                with open(densityFile,'a+') as f:
                    f.write("\t".join(("chr",str(chr),str(clusternum),str(type),str(v_H3K27ac),str(v_H3K4me3),str(v_CTCF),str(v_POLR2A),str(v_H3K9me3)))+'\n')
                    f.close()
        if clusternum=="three":
                for type in ["first","second","three"]:
                    feature_path='epi/output/'+cellline+'/chr'+str(chr)+'.txt'
                    dataAll = []
                    single_edge_feature = []
                    feature_19 = np.loadtxt(feature_path)
                    #data_19 = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')
                    f=open("output/"+cellline+"/hypergraph/"+clusternum+"/chr"+str(chr)+"_"+type+".txt")
                    for line in f.readlines():
                        start = int(line.split()[0])
                        end = int(line.split()[2])
                        # if start>=48 and end<=196:
                        dataSet = []
                        for i in np.arange(start, end+1):
                            dataSet.append(feature_19[i-1])
                        single_TAD_feature = np.sum(dataSet, axis=0)
                        dataAll.append(single_TAD_feature)
                    f.close
                    v_H3K27ac = np.sum(dataAll, axis=0)[0]
                    v_H3K4me3 = np.sum(dataAll, axis=0)[1]
                    v_CTCF = np.sum(dataAll, axis=0)[2]
                    v_POLR2A = np.sum(dataAll, axis=0)[3]
                    v_H3K9me3 = np.sum(dataAll, axis=0)[4]
                    densityFile="output/"+cellline+"/factors-density-metric-5.txt"
                    with open(densityFile,'a+') as f:
                        f.write("\t".join(("chr",str(chr),str(clusternum),str(type),str(v_H3K27ac),str(v_H3K4me3),str(v_CTCF),str(v_POLR2A),str(v_H3K9me3)))+'\n')
                        f.close()
        if clusternum=="four":
                for type in ["first","second","three","four"]:
                    feature_path='epi/output/'+cellline+'/chr'+str(chr)+'.txt'
                    dataAll = []
                    single_edge_feature = []
                    feature_19 = np.loadtxt(feature_path)
                   #data_19 = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')
                    f=open("output/"+cellline+"/hypergraph/"+clusternum+"/chr"+str(chr)+"_"+type+".txt")
                    for line in f.readlines():
                        start = int(line.split()[0])
                        end = int(line.split()[2])
                        # if start>=48 and end<=196:
                        dataSet = []
                        for i in np.arange(start, end+1):
                            dataSet.append(feature_19[i-1])
                        single_TAD_feature = np.sum(dataSet, axis=0)
                        dataAll.append(single_TAD_feature)
                    f.close
                    v_H3K27ac = np.sum(dataAll, axis=0)[0]
                    v_H3K4me3 = np.sum(dataAll, axis=0)[1]
                    v_CTCF = np.sum(dataAll, axis=0)[2]
                    v_POLR2A = np.sum(dataAll, axis=0)[3]
                    v_H3K9me3 = np.sum(dataAll, axis=0)[4]
                    densityFile="output/"+cellline+"/factors-density-metric-5.txt"
                    with open(densityFile,'a+') as f:
                        f.write("\t".join(("chr",str(chr),str(clusternum),str(type),str(v_H3K27ac),str(v_H3K4me3),str(v_CTCF),str(v_POLR2A),str(v_H3K9me3)))+'\n')
                        f.close()
        if clusternum=="five":
                for type in ["first","second","three","four","five"]:
                    feature_path='epi/output/'+cellline+'/chr'+str(chr)+'.txt'
                    dataAll = []
                    single_edge_feature = []
                    feature_19 = np.loadtxt(feature_path)
                   #data_19 = np.loadtxt('output/'+cellline+'/50kb_TAD/chr'+str(chr)+'.txt')
                    f=open("output/"+cellline+"/hypergraph/"+clusternum+"/chr"+str(chr)+"_"+type+".txt")
                    for line in f.readlines():
                        start = int(line.split()[0])
                        end = int(line.split()[2])
                        # if start>=48 and end<=196:
                        dataSet = []
                        for i in np.arange(start, end+1):
                            dataSet.append(feature_19[i-1])
                        single_TAD_feature = np.sum(dataSet, axis=0)
                        dataAll.append(single_TAD_feature)
                    f.close
                    v_H3K27ac = np.sum(dataAll, axis=0)[0]
                    v_H3K4me3 = np.sum(dataAll, axis=0)[1]
                    v_CTCF = np.sum(dataAll, axis=0)[2]
                    v_POLR2A = np.sum(dataAll, axis=0)[3]
                    v_H3K9me3 = np.sum(dataAll, axis=0)[4]
                    densityFile="output/"+cellline+"/factors-density-metric-5.txt"
                    with open(densityFile,'a+') as f:
                        f.write("\t".join(("chr",str(chr),str(clusternum),str(type),str(v_H3K27ac),str(v_H3K4me3),str(v_CTCF),str(v_POLR2A),str(v_H3K9me3)))+'\n')
                        f.close()