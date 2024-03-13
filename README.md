# A Method for chromatin domain partitioning based on hypergraph clustering
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

## step8. plot figures
run scripts of Landmark position.ipynb
