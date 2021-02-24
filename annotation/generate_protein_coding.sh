#!/bin/bash

output=$1
if [ -z "$output" ]; then
	output=protein-coding-mgi_$(date +"%d-%m-%Y").csv
fi

wget -q -O - http://informatics.jax.org/downloads/reports/MGI_EntrezGene.rpt | tr -d '"' | perl -F'\t' -lane 'BEGIN{print "MGI_ID,MGI_Symbol"} print (join(",", @F[0,1])) if ((($F[2] eq "O") && ($F[6] eq "Gene")) && ((index($F[14], "protein-coding") != -1) || (index($F[14], "protein_coding") != -1)))' > $output
