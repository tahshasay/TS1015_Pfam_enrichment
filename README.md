# TS1015_Pfam_enrichment

#These data and scripts are supplement to Say, T. E., & Degnan, S. M. (2019). Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512

#To ascertain the over-represented functional domains, a Pfam enrichment analysis was performed on the sPLS-DA gene list (see above) against a background reference (all genes in the genome) using a previously published R script (Chandran et al., 2009); see additional supplementary data here. 

#This folder contains the files and scripts used to perform the PFAM enrichment on the genes identified by the sPLS-DA.

INPUT_files:
'Aqu2.1_Pfam_annotation_final.txt' - pfam annotation file.

'pfam_description.txt' - pfam descriptions file.

'q03.01iii_TS1015_PC1_transcripts_2017.07.25_AQlist.txt' - list of sPLS-DA genes identified on the 1st component. 

OUTPUT_files:
Pfam enrichment plot is in Figure S3.

Reference

Chandran, D., Tai, Y. C., Hather, G., Dewdney, J., Denoux, C., Burgess, D. G., . . . Wildermuth, M. C. (2009). Temporal global expression data reveal known and novel salicylate-impacted processes and regulators mediating powdery mildew growth and reproduction on Arabidopsis. Plant Physiology, 149(3), 1435-1451. doi:10.1104/pp.108.132985
