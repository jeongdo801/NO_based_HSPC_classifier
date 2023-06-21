# HSPC cell-type classifier using single-cell NO profiles
This is the script developed for the provisional paper Grimes and Jeong et al. titled "Cell type-specific consequences of mosaic structural variants in hematopoietic stem and progenitor cells"



## pre-requirement step
1. preparation of single-cell NO information using scNOVA analysis <br>
scNOVA: [[https://github.com/friendsofstrandseq/mosaicatcher-pipeline](https://github.com/jeongdo801/scNOVA)](https://github.com/jeongdo801/scNOVA)
> Jeong and Grimes *et al.*, 2022 (doi: 10.1038/s41587-022-01551-4)

<br/><br/>
## Setup and run the script
1. **input files**
	* Single-cell level NO table : `result/{SAMPLE}_sort_geneid.txt`
2. **models**
	* scripts/matlab_scMNase_final_model_BM_github.mat
	* scripts/matlab_scMNase_final_model_CB_github.mat
3. **scripts**
	* scripts/scNOVA_celltype_prediction_HSPC_BM_classifier.m
	* scripts/scNOVA_celltype_prediction_HSPC_CB_classifier.m
4. **main outcomes**
  * It classifies single-cells into one of eight HSPC cell-types based on NO at gene body. These cell types include: hematopoietic stem cells (HSCs), multipotent progenitors (MPPs), lymphoid-primed multipotent progenitors (LMPPs), common lymphoid progenitors (CLPs), plasmacytoid-dendritic cells (pDC), common myeloid progenitors (CMPs), granulocyte-macrophage progenitors (GMPs), and megakaryocyte–erythroid progenitors (MEPs) 
