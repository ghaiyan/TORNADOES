#!/bin/bash
juicer='/data2/ghy_data/ghy_data/juicer'
cellline='IMR90'
rawdatapath=/data1/ghy_data/IMR90
norm_kr_contact=/data1/ghy_data/IMR90/1Mb
for j in {1..22}
do
java -jar $juicer/juicer_tools_1.22.01.jar  dump observed KR $rawdatapath/4DNFIH7TH4MF.hic  ${j} ${j}   BP 1000000 ${norm_kr_contact}/chr${j}_kr_1000000.hic
#java -jar /home/software/juicer/CPU/juicer_tools.jar  eigenvector -p KR inter_30.hic ${j} BP 100000 ${1}_chr${j}_pc1_100k.txt
#java -jar /home/software/juicer/CPU/juicer_tools.jar  dump observed NONE inter_30.hic  ${j} ${j}   BP 50000 ${1}_chr${j}_raw_50k.txt
done