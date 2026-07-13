# FUNCTIONAL PROFILING USING PICRUSt2
# PICRUSt2 is used for functional prediction.
# A QIIME2 plugin tutorial is available at:
# https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial
#
# Installation instructions for the QIIME2 plugin are available at:
# https://library.qiime2.org/plugins/q2-picrust2/13/
#
# Alternatively, PICRUSt2 can be installed using Conda or from source,
# following the instructions available at:
# https://github.com/picrust/picrust2/wiki/Installation
#
# Note that Conda must be installed before using either installation method.

# Once PICRUSt2 has been installed, activate the corresponding Conda environment.

# 1. Sequence placement
place_seqs.py -s refseqs.fna -o placed_seqs.tre -p 1 --intermediate placement_working
# The refseqs.fna file should be exported from the phyloseq object generated previously.

# 2. Hidden-state prediction
hsp.py -i 16S -t placed_seqs.tre -o marker_nsti_predicted.tsv.gz -p 1 -n
hsp.py -i EC -t placed_seqs.tre -o EC_predicted.tsv.gz -p 1
hsp.py -i KO -t placed_seqs.tre -o KO_predicted.tsv.gz -p 1

# 3. Metagenome prediction
metagenome_pipeline.py -i study_seqs.biom \
                       -m marker_nsti_predicted.tsv.gz \
                       -f EC_predicted.tsv.gz \
                       -o EC_metagenome_out

metagenome_pipeline.py -i study_seqs.biom \
                       -m marker_nsti_predicted.tsv.gz \
                       -f KO_predicted.tsv.gz \
                       -o KO_metagenome_out

# 4. Pathway abundance inference
pathway_pipeline.py -i EC_metagenome_out/pred_metagenome_unstrat.tsv.gz \
                    -o pathways_out \
                    --intermediate minpath_working \
                    -p 1
