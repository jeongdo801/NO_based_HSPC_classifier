# HSPC cell-type classifier using single-cell NO profiles
This is the script developed for the provisional paper Grimes and Jeong et al. titled "Cell type-specific consequences of mosaic structural variants in hematopoietic stem and progenitor cells"


## Overview of this workflow
STEP1. Pre-requirement step - preparation of single-cell NO information using scNOVA analysis <br>
scNOVA: [https://github.com/friendsofstrandseq/mosaicatcher-pipeline](https://github.com/jeongdo801/scNOVA)

<br/><br/>
# Setup
1. **input files**
	* git lfs install
	* git clone https://github.com/jeongdo801/scNOVA.git
        * install dependencies (see further below)
2. **models**
	* Add your single-cell bam and index files (input_bam/*.bam)
	* Add key result files from mosaicatcher output in the input_user folder
		* input_user/simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE.txt
		* input_user/stran<br/><br/>
3. **main outcomes**
  * These cell types include: hematopoietic stem cells (HSCs), multipotent progenitors (MPPs), lymphoid-primed multipotent progenitors (LMPPs), common lymphoid progenitors (CLPs), plasmacytoid-dendritic cells (pDC), common myeloid progenitors (CMPs), granulocyte-macrophage progenitors (GMPs), and megakaryocyte–erythroid progenitors (MEPs) 
