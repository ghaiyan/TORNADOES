# TORNADOES
A Method for chromatin domain partitioning based on hypergraph clustering
# Before processing
Before running code, please use the method CASPIAN to generate TAD files with different resolutions. CASPIAN code can be accessed from https://gitee.com/ghaiyan/caspian.
# step1: process ChIP-seq data and generate A/B compartment using other tools
running the code in the file of data_process.ipynb
# step2: generate the hypergraph and classify TADs into different types: two, three or four
running the code in the file of hypergraph.ipynb using the jupter notebook.
# step3: evaluate
running the code in the file of Landmark position.ipynb
