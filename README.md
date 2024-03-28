# A Method for chromatin domain partitioning based on hypergraph clustering

![image](https://github.com/ghaiyan/TORNADOES/assets/27917393/18d3985b-ba0a-4070-813d-70339cdbd22c)

Figure 1. Workflow of the TORNADOES method. (a) Workflow of analyzing the chromatin domains. (b) Workflow of generating the clusters based on hypergraph learning.
## datasets: Hi-C and ChIP-seq data.
Cell line	Data type	Data linker

**GM12878	:**

Hi-C	https://data.4dnucleome.org/files-processed/4DNFI1UEG1HD/

H3K27ac ChIP-seq	https://www.encodeproject.org/files/ENCFF087YCU/@@download/ENCFF087YCU.bigWig
https://www.encodeproject.org/files/ENCFF367KIF/@@download/ENCFF367KIF.bed.gz

H3K4me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF480KNX/@@download/ENCFF480KNX.bigWig
https://www.encodeproject.org/files/ENCFF188SZS/@@download/ENCFF188SZS.bed.gz

 H3K9me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF533NIQ/@@download/ENCFF533NIQ.bigWig
https://www.encodeproject.org/files/ENCFF682WIQ/@@download/ENCFF682WIQ.bed.gz

 CTCF ChIP-seq	https://www.encodeproject.org/files/ENCFF680XUD/@@download/ENCFF680XUD.bigWig
https://www.encodeproject.org/files/ENCFF827JRI/@@download/ENCFF827JRI.bed.gz

 POLR2A ChIP-seq	https://www.encodeproject.org/files/ENCFF532WJA/@@download/ENCFF532WJA.bigWig
https://www.encodeproject.org/files/ENCFF886PSD/@@download/ENCFF886PSD.bed.gz

 RAD21 ChIP-seq	https://www.encodeproject.org/files/ENCFF940IMA/@@download/ENCFF940IMA.bed.gz

 SMC3 ChIP-seq	https://www.encodeproject.org/files/ENCFF837YJA/@@download/ENCFF837YJA.bed.gz

 EZH2 ChIP-seq	https://www.encodeproject.org/files/ENCFF339YTO/@@download/ENCFF339YTO.bed.gz

 H3K36me3	https://www.encodeproject.org/files/ENCFF268HMO/@@download/ENCFF268HMO.bed.gz

 H3K4me1	https://www.encodeproject.org/files/ENCFF321BVG/@@download/ENCFF321BVG.bed.gz
 
**H1-hESC:**

Hi-C	https://data.4dnucleome.org/files-processed/4DNFI2TK7L2F/

 H3K27ac ChIP-seq	https://www.encodeproject.org/files/ENCFF390JIZ/@@download/ENCFF390JIZ.bigWig
https://www.encodeproject.org/files/ENCFF045CUG/@@download/ENCFF045CUG.bed.gz

 H3K4me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF301RJN/@@download/ENCFF301RJN.bigWig
https://www.encodeproject.org/files/ENCFF277AOQ/@@download/ENCFF277AOQ.bed.gz

 H3K9me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF562HZQ/@@download/ENCFF562HZQ.bigWig
https://www.encodeproject.org/files/ENCFF918VFL/@@download/ENCFF918VFL.bed.gz

 CTCF ChIP-seq	https://www.encodeproject.org/files/ENCFF147GRN/@@download/ENCFF147GRN.bigWig
https://www.encodeproject.org/files/ENCFF821AQO/@@download/ENCFF821AQO.bed.gz

 POLR2A ChIP-seq	https://www.encodeproject.org/files/ENCFF942TZX/@@download/ENCFF942TZX.bigWig
https://www.encodeproject.org/files/ENCFF322DAE/@@download/ENCFF322DAE.bed.gz

 RAD21 ChIP-seq	https://www.encodeproject.org/files/ENCFF883FUW/@@download/ENCFF883FUW.bed.gz

 H3K36me3	https://www.encodeproject.org/files/ENCFF504KOV/@@download/ENCFF504KOV.bed.gz

 H3K4me1	https://www.encodeproject.org/files/ENCFF006PXB/@@download/ENCFF006PXB.bed.gz

**K562	:**

Hi-C	https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/dcfcb009-f006-4ab8-a4c7-af72be58c12c/4DNFITUOMFUQ.hic

H3K27ac ChIP-seq	https://www.encodeproject.org/files/ENCFF094XCU/@@download/ENCFF094XCU.bigWig
https://www.encodeproject.org/files/ENCFF038DDS/@@download/ENCFF038DDS.bed.gz

H3K4me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF845BFY/@@download/ENCFF845BFY.bigWig
https://www.encodeproject.org/files/ENCFF909PMV/@@download/ENCFF909PMV.bed.gz

H3K9me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF559MMQ/@@download/ENCFF559MMQ.bigWig
https://www.encodeproject.org/files/ENCFF371GMJ/@@download/ENCFF371GMJ.bed.gz

CTCF ChIP-seq	https://www.encodeproject.org/files/ENCFF268ZPN/@@download/ENCFF268ZPN.bigWig
https://www.encodeproject.org/files/ENCFF221SKA/@@download/ENCFF221SKA.bed.gz

 POLR2A ChIP-seq	https://www.encodeproject.org/files/ENCFF914WIS/@@download/ENCFF914WIS.bigWig
https://www.encodeproject.org/files/ENCFF355MNE/@@download/ENCFF355MNE.bed.gz

 RAD21 ChIP-seq	https://www.encodeproject.org/files/ENCFF930WPG/@@download/ENCFF930WPG.bed.gz

 SMC3 ChIP-seq	https://www.encodeproject.org/files/ENCFF289LLT/@@download/ENCFF289LLT.bed.gz

 EZH2	https://www.encodeproject.org/files/ENCFF804RVA/@@download/ENCFF804RVA.bed.gz

 H3K36me3	https://www.encodeproject.org/files/ENCFF053DAC/@@download/ENCFF053DAC.bed.gz

 H3K4me1	https://www.encodeproject.org/files/ENCFF759NWD/@@download/ENCFF759NWD.bed.gz
 
**IMR90:**

Hi-C	https://data.4dnucleome.org/files-processed/4DNFIH7TH4MF/

 H3K27ac ChIP-seq	https://www.encodeproject.org/files/ENCFF907GKJ/@@download/ENCFF907GKJ.bigWig
https://www.encodeproject.org/files/ENCFF730BVO/@@download/ENCFF730BVO.bed.gz

 H3K4me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF811ZFE/@@download/ENCFF811ZFE.bigWig
https://www.encodeproject.org/files/ENCFF018CAH/@@download/ENCFF018CAH.bed.gz

 H3K9me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF733CJA/@@download/ENCFF733CJA.bigWig
https://www.encodeproject.org/files/ENCFF098XMT/@@download/ENCFF098XMT.bed.gz

 CTCF ChIP-seq	https://www.encodeproject.org/files/ENCFF895EDK/@@download/ENCFF895EDK.bigWig
https://www.encodeproject.org/files/ENCFF307XFM/@@download/ENCFF307XFM.bed.gz

 POLR2A ChIP-seq	https://www.encodeproject.org/files/ENCFF444OMD/@@download/ENCFF444OMD.bigWig
https://www.encodeproject.org/files/ENCFF448ZOJ/@@download/ENCFF448ZOJ.bed.gz

 H3K36me3	https://www.encodeproject.org/files/ENCFF449ADN/@@download/ENCFF449ADN.bed.gz

 H3K4me1	https://www.encodeproject.org/files/ENCFF611UWF/@@download/ENCFF611UWF.bed.gz

**HepG2:**

Hi-C	https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/25104375-a588-46e6-a382-663cee6c332f/4DNFICSTCJQZ.hic

 H3K27ac ChIP-seq	https://www.encodeproject.org/files/ENCFF847MVG/@@download/ENCFF847MVG.bigWig
https://www.encodeproject.org/files/ENCFF886SZT/@@download/ENCFF886SZT.bed.gz

 H3K4me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF359LQU/@@download/ENCFF359LQU.bigWig
https://www.encodeproject.org/files/ENCFF185WOB/@@download/ENCFF185WOB.bed.gz

 H3K9me3 ChIP-seq	https://www.encodeproject.org/files/ENCFF284AIG/@@download/ENCFF284AIG.bigWig
https://www.encodeproject.org/files/ENCFF353MWR/@@download/ENCFF353MWR.bed.gz

 CTCF ChIP-seq	https://www.encodeproject.org/files/ENCFF105VBF/@@download/ENCFF105VBF.bigWig
https://www.encodeproject.org/files/ENCFF612ZUY/@@download/ENCFF612ZUY.bed.gz

 POLR2A ChIP-seq	https://www.encodeproject.org/files/ENCFF425QWO/@@download/ENCFF425QWO.bigWig
https://www.encodeproject.org/files/ENCFF354VWZ/@@download/ENCFF354VWZ.bed.gz

 RAD21 ChIP-seq	https://www.encodeproject.org/files/ENCFF304NRB/@@download/ENCFF304NRB.bed.gz

 SMC3 ChIP-seq	https://www.encodeproject.org/files/ENCFF028OBL/@@download/ENCFF028OBL.bed.gz

 H3K36me3	https://www.encodeproject.org/files/ENCFF489ZNJ/@@download/ENCFF489ZNJ.bed.gz

 H3K4me1	https://www.encodeproject.org/files/ENCFF532AGT/@@download/ENCFF532AGT.bed.gz

# Main steps:
1. Download the data set
   
3. Extract Hi-C data
   
4. Identify 2 types of TAD at different degrees

5. Identify sub-compartment in calder
   
6. Identify AB area with fanc tool

7. Process AB data
    
8. process and generate the result

9. The proportion of different chip-seq peaks anchored by different classes
    
10. Different types of chip-seq signal values

# all experiment data and code
all experiment data and code can be downloaded from http://mged.nmdms.ustb.edu.cn/storage/data/28658183.
## step1.The Hi-C matrix with resolution of 1Mb, 50kb and 25kb was extracted by juicer tool.
The IMR90 data sets for example, hic format file to download from https://data.4dnucleome.org/files-processed/4DNFIH7TH4MF/.
generate tuple format:
bash example/0-generateKR_hic_matrix.sh

convert tuple format into N*N matrxi format:
run scripts of example/1-tuple2matrix.ipynb

## step2.use the method CASPIAN to generate TAD files with different resolutions. CASPIAN code can be accessed from https://gitee.com/ghaiyan/caspian.
or you can run the code 
run scripts of example/2-caspian.ipynb


## step3.generate the A/B compartments using fanc tool.
you can use fanc tool by accessing https://vaquerizaslab.github.io/fanc/api/analyse/compartments.html 

## step4.generate the sub-compartments using CALDER tool.

you can use CALDER tool by accessing https://github.com/CSOgroup/CALDER to get the chrX_sub_compartments.bed using Hi-C file with (start_location, end_location, contact_value).


library(CALDER)

for (chr in 1:22) {

  chr_filename <- paste0("/home/rstudio/calder/gm12878/chr", chr, "_vc_50kb.hic")

  out_dir <- paste0("./gm12878/chr", chr)
  
  CALDER_main(chr_filename, chr=chr, bin_size=50E3, out_dir=out_dir, sub_domains=TRUE, save_intermediate_data=FALSE)

}



## step5: process the epi files.

run analysis_epi.py.

## step6.generate or process the hypergraph and A/B compartments, sub-compartments.
run scripts of example/3-data_process.ipynb

## step7.Hypergraph learning was used to generate sub-compartments, and the similarity was compared with fanc and CALDER A/B compartment and su-compartment.
run scripts of example/4-hypergraph.ipynb

Or you can change the data of code "example/3-dataprocess.py", and run "python example/3-dataprocess.py"

## step8. plot figures
run scripts of Landmark position.ipynb

## step9 evaluete 
run scripts of "example/5-evalute-anchor-factors.py" and "example/6-evalute-factors-density.py"

