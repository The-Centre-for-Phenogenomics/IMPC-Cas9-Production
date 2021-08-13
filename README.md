# Files and Code for Cas9 Production Manuscript

_Hillary Elrick, Elif Acar_

_August 12, 2021_

## Directory Structure

* annotation - _annotation files_
	* genelists
		* MGI_ENSID_Mapping.xlsx
		* MGI_Entrez_Mapping.xlsx
		* MGI_Genelist.txt
		* protein-coding-mgi_14-02-2021.csv - _(protein-coding genes generated by `generate_protein_coding.sh`)_
	* GEO - _embryonic expression data_
		* annotation/GEO/GPL1261-56135.txt
		* annotation/GEO/GSE11224_series_matrix.txt
	* histone_modification
		* mm10_H3K27ac_annotatedpeaks.tsv
		* mm10_H3K27me3_annotatedpeaks.tsv
	* human_orthologs - _annotations from human ortholog genes_
		* gnomad_v2_1_1_lof_metrics_by_gene.xlsx
		* Human_mouse_orthologues_for_IMPC.xlsx 
		* HumanEssentiality_MouseOrthologs.xlsx
		* omim_to_ensembl.tsv
	* no_founder_reasons - _files generated by centres explaining GLT failures that had founders_
		*  YYYYMMDD_ExperimentsWithFoundersNoGLT_X.xlsx - _where X is the centre code_
	* mi_attempts_list_20201009.xlsx
	* all_genes_viability.csv _viability info for genes from IMPC release 13.0_
	* AnnotationInfo.xlsx - _gene calculations: length, GC contents, CpG sites, CpG percentage_
	* cytoBand.txt - _cytogenic banding information from UCSC_
	* mm10.chrom.sizes - _mm10 chromosome sizes from UCSC_

* data - _data files_
	* Clean-IMPC_Cas9_2020-10-09.csv _data file containing the latest attempt of each edited gene_
	* ST3.csv _data file containing the repeated attempts of edited genes by each production centre_

* python - _python code files_
	* clean_data.py - _cleaning script_
	* add_annotation_info.py - _add supplementary annotation information_

* R - _R code files_
	* Cas9Production_Stats_Analysis.Rnw _statistical analysis of data_
	* Supplementary_GLM_Analysis.Rmd _supplementary analysis_

## Running

### Requirements
* python>=3.8.5
* pandas
* numpy
* R

On command-line, run:
```bash
python3 python/clean_data.py --input=<iMits-Cas9-YYYY-MM-DD>.xlsx --output=<output_directory>
```
Replacing `<iMits-Cas9-YYYY-MM-DD>.xlsx` with your iMits download file and `<output_directory>` with an empty directory (will be created if it doesn't exist)
**N.B. Before running, edit the spreadsheet from iMits to fix the headers and remove the 'Status Name' column from the F1 Colonies Tab**

This will run all the cleaning steps and write the following files into your output_directory:
* clean data csv
* clean data xlsx
* clean data repeats csv
* clean data repleast xlsx
* summary of all steps performed

Then, to add additional annotation information for the supplementary R anaylsis into failures, run:
```bash
python3 python/add_annotation_info.py <previous-step-output>/Clean-<iMits-Cas9-YYYY-MM-DD>.xlsx <OutputFileName>.xlsx
```

Finally, to run:
* Statistical Analysis,  open `Cas9Production_Stats_Analysis.Rnw` in RStudio, set working directory to source pane, install the listed packages, download the required data files, and Knit the document. 
* Supplementary Analysis, open `Supplementary_GLM_Analysis.Rmd` in RStudio, install the listed packages, and Knit the document.

