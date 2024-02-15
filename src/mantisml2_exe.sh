#!/bin/bash

#SBATCH --job-name=mantisml2_asthma
#SBATCH --output=/scratch/gen1/nnp5/Gene_analysis_tmp/logerror/%x-%j.out
#SBATCH --time=1:0:0
#SBATCH --mem=30gb
#SBATCH --account=gen1
#SBATCH --export=NONE



source /home/n/nnp5/miniconda3/bin/activate /home/n/nnp5/miniconda3/envs/mantis_ml_2
mantisml2 -d "asthma" -o input/mantisml_asthma_out_HPO -g OT

#Rationale: use mantisml2 to retrieve asthma genes and compare them to the 98 genes from V2G.


#workdir: Gene_analysis/

#set up MANTIS-ML2 as here: https://github.com/astrazeneca-cgr-publications/mantis-ml-release-2.0/tree/main

#Genomics England does not have asthma as a term
#HPO has: 'Asthma'
#OT has: ‘Chronic Obstructive Asthma’, ‘asthma’, ’childhood onset asthma’

#1Run for each four terms
#2Compare results
#3Compare results with the 98 genes with an enrichment analysis: are the 98 genes
#4Enriched in mantisml2-derived asthma genes?


#mantisml2 -d "asthma" -o input/mantisml_asthma_out_OT -g OT