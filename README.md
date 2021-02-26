# Files and Code for Cas9 Production Manuscript

## Directory Structure
* input - _contains input files_
	* iMits-Cas9_YYYY-MM-DD.xlsx - _Cas9 data exported from iMits on day specified_

* annotation - _annotation files_
	* IMPC-viability_release-13.0.csv - _viability info for genes from IMPC_
	* MGI-ProteinCoding_YYYY-MM-DD.csv - _protein coding genes (as defined by MGI)_
	* human_orthologs/ - _annotations from human ortholog genes_
		* gnomad_v2_1_1_lof_metrics_by_gene.xlsx - _gnomad download_
		* HumanEssentiality_MouseOrthologs.xlsx - _essentiality data_
		* Human_mouse_orthologues_for_IMPC.xlsx - _ortholog mapping_
	* no_founder_reasons/ - _files generated by centres explaining GLT failures that had founders_
		*  YYYYMMDD_ExperimentsWithFoundersNoGLT_X.xlsx - _where X is the centre code_

* python - _python code files_
	* annotate.py - _add annotation information_
	* clean_data.py - _cleaning script_

* R - _R code files_
	* ChIP_SeqAnalysis.Rmd - _supplementary analysis with ChIPSeq data__

## Running

### Requirements
* python>=3.6
* pandas
* R

On command-line, run:
```bash
python3 python/clean_data.py --input=<iMits-Cas9-YYYY-MM-DD>.xlsx --output=<output_directory>
```
Replacing `<iMits-Cas9-YYYY-MM-DD>.xlsx` with your iMits download file and `<output_directory>` with an empty directory (will be created if it doesn't exist)

This will run all the cleaning steps and write the following files into your output_directory:
* clean data csv
* clean data xlsx
* clean data repeats csv
* clean data repeats xlsx
* summary of all steps performed

**N.B. Before running a new export, edit the spreadsheet to remove the 'Status Name' column from the F1 Colonies Tab**
**Additionally, the header columns may need to be manually merged so they're all the same "level"**