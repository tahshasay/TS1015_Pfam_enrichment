## This code identifies Pfam domains that are enriched among a selected set of genes.
## Created 6-24-08 by Gregory Hather.
##
## The inputs are:
## "Aq_All.txt" -- A list of all the genes considered (AGI numbers). In our case, it was all the chromosomal genes represented on the ATH1 Affymetrix microarray.
## "Pfam_annotation_final.txt" -- A table with two columns. The first column is the AGI number, and the second column is the Pfam domain.
## "pfam_description.txt" -- A table with two columns. The first column is the Pfam domain, and the second column is the description of that domain.
## "genes-with-no-conversion.txt" -- A list of AGI numbers which should be excluded because there was AGI2Uniprot conversion.
## 1.txt -- A list of AGI numbers for the subset of gene of interest.
##
## The output is:
## output_file -- a table listing the enriched domains
##
##
## Note:
## "converted-swisspfam-arath.txt", "swisspfam-descriptions.txt", and "genes-with-noconversion.txt" were generated using perl and python scripts.
## The inputs to the scripts were swisspfam (from Pfam 22.0) and AGI2Uniprot.20080418 (from TAIR).
## All the genes in "converted-swisspfam-arath.txt" are in "ATH1all.txt".
## No genes from "converted-swisspfam-arath.txt" are in "genes-with-noconversion.txt".


################################################ to change / 1 up from  output folder
setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/pfam_enrichment");
################################################ to change

# refer to the two new files below - now with Aqu2.1 model naming system
ATH1all = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_rmtext_Aq_all.txt", what = list(""))[[1]]


################################################ to change
# conduct a pfam enrichment using only the "unique" pfam domains within each Aqu2.1 model
#converted_swisspfam = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_Pfam_annotation_final.txt", what = list("",""), sep = "\t")
#
# total abundance of pfam domains 
converted_swisspfam = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/Aqu2.1_Pfam_annotation_final.txt", what = list("",""), sep = "\t")
################################################ to change



# no need to rename below
swisspfam_descriptions = scan("/Users/tahshasay/Documents/Aqu2.1_Pfam/WH_PFAM_enrichment_2016.10.25_renamed/pfam_description.txt", what = list("",""), sep ="\t", quote = "")

# 2016.10.25
# input files- listAq - copied from below location (using new v R and ~group analysis, and lists fo Aq models were extracted from VENNY diagrams (hence incl. shared etc). 
######################
#/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/DESeq2_TS1015_h03.02_~group_2016.10.24_newRv/TS1015_VENNY_lists/TS1015_t1.t9_group_adjp0.05_FC/t1.t9_group_adjp_0.05_posFC_nat-nat.txt
######################



# to change: DEG1248
# to change: uniq
################################################ to change
top_genes = scan("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/TS1015_PC1_transcripts_2017.07.25_AQlist.txt", what = list(""))[[1]]
output_file = "TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04.txt"
################################################ always write as .txt file (file layout)


## split genes by domains
domain_grouping = split(converted_swisspfam[[1]], converted_swisspfam[[2]])

## Do a hypergeometric test for enrichment of each domain among top_genes
## The balls in the urn are all the genes that are in ATH1all and not in genes_with_no_conversion.
## q: Genes with given domain that are in the top set are the white balls drawn from the urn.
## m: Genes with given domain are white balls.
## n: Genes without given domain are black balls.
## k: The balls drawn are the genes in the top set.

## loop over each domain
p_value_vector = c()
enrichment_vector = c()
q_vector = c()
m_vector = c()
for (index in 1:length(domain_grouping)){
domain = names(domain_grouping)[index]
genes_with_domain = domain_grouping[[index]]
q = length(as.vector(na.omit(match(genes_with_domain, top_genes))))
m = length(genes_with_domain)
n = length(ATH1all) - m
k = length(top_genes)
p_value = phyper(q - 1, m, n, k, lower.tail = FALSE)
p_value_vector[index] = p_value
enrichment = q/(m*k/(n+m))
enrichment_vector[index] = enrichment
q_vector[index] = q
m_vector[index] = m
}

## adjust p-values
library(multtest)
adjusted_p = mt.rawp2adjp(p_value_vector, proc = c("BH"))

## write result
p_cut_off = 0.1
significant_indices = as.vector(na.omit(ifelse(adjusted_p$adjp[,2] < p_cut_off, adjusted_p$index, NA)))
significant_enrichments = enrichment_vector[significant_indices]
indices_by_enrichment = significant_indices[sort(significant_enrichments, decreasing = TRUE, index.return = TRUE)[[2]]]
sorted_descriptions = swisspfam_descriptions[[2]][match(names(domain_grouping)[indices_by_enrichment], swisspfam_descriptions[[1]])]
sorted_Pfam_ID = names(domain_grouping)[indices_by_enrichment]
sorted_full_descriptions = paste(sorted_descriptions, " (", sorted_Pfam_ID, ")", sep = "")
sorted_enrichment = enrichment_vector[indices_by_enrichment]
sorted_p_values = adjusted_p$adjp[match(indices_by_enrichment, adjusted_p$index),2]
sorted_q = q_vector[indices_by_enrichment]
sorted_m = m_vector[indices_by_enrichment]
length(top_genes)
length(ATH1all)
sorted_gene_groups = c()
for (index in 1:length(indices_by_enrichment)){
temp_genes = top_genes[as.vector(na.omit(match(domain_grouping[[indices_by_enrichment[index]]], top_genes)))]
sorted_gene_groups[index] = paste(temp_genes, sep = ", ", collapse = ", ")
}
output = cbind(sorted_full_descriptions, round(sorted_enrichment, digits = 1), as.character(format(signif(sorted_p_values, digits = 1), scientific = TRUE)), sorted_q, sorted_m, sorted_gene_groups)

write(t(output), file = output_file, sep = "\t", ncolumns = 6)


writeLines(capture.output(sessionInfo()), "sessionInfo.txt")



 

####################################################################################
# Visualising data: 
# Dot plot - pfam enrichment 
#
# First saw this dot plot here:
#http://www.bioconductor.org/packages/devel/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html
# Statistical analysis and visualization of functional profiles for genes and gene clusters
# Code cited above adapted by T.Say for CEL-seq2 data
#
####################################################################################

# Goal: "make a geom_point plot of the pfam domains to easily compared similarities and differences between gene lists. 

library(ggplot2)

setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/pfam_enrichment")


#-----------------------------------------------
# Importing multiple .csv / .txt files into R
#-----------------------------------------------

# http://stackoverflow.com/questions/11433432/importing-multiple-csv-files-into-r

library(plyr)

# conduct loop to read in all .csv files listed above
# for (i in 1:length(temp)) assign(temp[i], read.csv(temp[i]))
# for (i in 1:length(lists)) assign(lists[i], read.csv(lists[i]))

# Or, without assign, and to demonstrate (1) how the file name can be cleaned up and (2) show how to use list2env, you can try the following:
#list2env(
 # lapply(setNames(temp, make.names(gsub("*.csv$", "", temp))), 
  #       read.csv), envir = .GlobalEnv)


#lists <- list.files("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.04_DESeq2_DEA_Design-group_2016.10.31/VENNY_lists/VEN_binary_data_for_R/annotated_lists_pfam/pfam_abund/ven_adjp0.1_posFC")

# alt
# This assumes that you have those CSVs in a single directory--your current working directory--and that all of them have the lower-case extension .csv.
lists = list.files(pattern="*.txt")

# check
lists

#list2env(
 #lapply(setNames(lists, make.names(gsub("*.csv$", "", lists))), 
        #read.csv, header=TRUE, row.names=1), envir = .GlobalEnv)      
        

# imp - no header present! added below bit.          
list2env(
 lapply(setNames(lists, make.names(gsub("*.txt$", "", lists))), 
        read.table, sep='\t'), envir = .GlobalEnv)      
        # header = TRUE if you want it to keep a space for you - but no header in these files

#----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #-----       


# adjp0.1_t9_posFC_lgt.drk_and_lgt.nat > adjp0.1_t9_posFC_lgt.drk_and_lgt.nat ((22))
# do not include part of file name referred to above. 

head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1)

# add collumn names manually to each data frame
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/colnames.html
colnames(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1) <- c("pfam", "enrichment", "pval", "Count", "noGenome", "genes")
head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1)


write.table(as.data.frame(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1),file="../pfam_enrichment_w_header/TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

# add column with test information
TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1$test <- "cluster1_upt1"

head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1) # check to see new column


#----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #----- #-----       
# to change



head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9)

# add collumn names manually to each data frame
# https://stat.ethz.ch/R-manual/R-devel/library/base/html/colnames.html
colnames(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9) <- c("pfam", "enrichment", "pval", "Count", "noGenome", "genes")
head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9)
TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9

write.table(as.data.frame(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9),file="../pfam_enrichment_w_header/TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

# add column with test information
TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9$test <- "cluster2_upt9"

head(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9) # check to see new column


#-----------------------------------------------------------------------------------------
#------------------------------------------------------------
# append tables together
# ie. each of the enrichment results into the one data matrix
#------------------------------------------------------------ 


# remove any with 0 enrichments
# adjp0.1_t9_posFC_lgt.nat,

# order matters!!
# order will appear on plot
splsda <- rbind(TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t1, TS1015_PC1_transcripts_2017.07.25_splsda_abd_pfam_2017.08.04_t9)


#######################
# to change
# # and manually create folder in wd
#######################
#write.table(as.data.frame(splsda),file="../pfam_enrichment_l00.00_heatmap/splsda_2017.08.04_abd_pfam_clusters.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')
 write.table(as.data.frame(splsda),file="../pfam_enrichment_dotplot/splsda_2017.08.04_abd_pfam_clusters.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

# if using row.names = TRUE do below
# manually F and R )1, )2, )3 with ) = all the same
# Error: duplicate row names are not allowed
# for some reason the "same" domains are kept separate.  
# therefore need to add a "name" to row.names column ie. "pfam" 
# no need to do above when you eliminate row.names by adding a header for this column


# all.trts <- read.table("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.04_DESeq2_DEA_Design-group_2016.10.31/select_lfc_2017.07.13/lfc1/TS1015_t1.t9.all.trts_resOrdered_Sig_adjp0.1_lfc1_posFC_2017.07.13_eds.txt", header=TRUE, sep="\t")

# subset data by pval
sig <- subset(splsda, pval < 0.05)


#######################
# to change
# # and manually create folder in wd
#######################
write.table(as.data.frame(sig),file="../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_pval0.05.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')
 

######################## colour options

#cols <- colorRampPalette(c("darkblue","purple","red"), space = "rgb")
#myPal <- cols(11) 

 	#scale_fill_gradient(values = myPal) # continuous https://stackoverflow.com/questions/31560916/colour-palette-with-ggplot
 	#scale_fill_manual(values = myPal) # discrete +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	#scale_size(range = c(0.1, 1)) +

library(viridis)
#hmcol <- rev(viridis(256))
#hmcol <- magma(256)
#hmcol <- rev(plasma(256))

# for "script" used to creat viridis scale go to:
# https://gist.github.com/hrbrmstr/f31899d067300d621baa



## WORKING
pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Count_tidy.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = Count, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA, 
       title="pfam_enrichment_abd_posFC") +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 1)) +
	theme(axis.text.y = element_text(size = 3), 
		axis.ticks = element_blank(), 
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(size = 3, angle = 90),
		axis.title.y = element_text(size = 6, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 6),
		title = element_text(size = 4),
		legend.text=element_text(size=3), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.5, "lines"), # legend height
        	plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()




## WORKING
pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Count_tidy_pval0.05.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = Count, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA, 
       title="pfam_enrichment_abd_posFC") +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 1)) +
	theme(axis.text.y = element_text(size = 3), 
		axis.ticks = element_blank(), 
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(size = 3, angle = 90),
		axis.title.y = element_text(size = 8, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 8),
		title = element_text(size = 4),
		legend.text=element_text(size=3), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.5, "lines"), # legend height
        	plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()


## BEST!!!!
pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Enrichment_adjp0.05_tidy.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = enrichment, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA, 
       title="pfam_enrichment_abd_posFC") +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 2.8)) + # size of dots
	theme(axis.text.y = element_text(size = 6), # x axis (pfam) domain
		axis.ticks = element_blank(), 
		panel.grid.major.x = element_blank(),
		axis.text.x = element_text(size = 8, angle = 90),
		axis.title.y = element_text(size = 8, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 8),
		title = element_text(size = 4),
		legend.text=element_text(size=5), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.7, "lines"), # legend height
        	plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()


pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Enrichment_adjp0.05_tidy_nobkgrd.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = enrichment, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA, 
       #title="pfam_enrichment_abd_posFC"
       colour = "P-value",
       size = "Enrichment score"
       ) +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 2.8)) + # size of dots
	theme(axis.text.y = element_text(size = 6), # x axis (pfam) domain
		axis.ticks = element_blank(), 
		panel.grid.major.y = element_line(colour = "grey84"),
		panel.grid.major.x = element_blank(),
		#panel.background = element_blank(),
		axis.text.x = element_text(size = 12, angle = 90),
		axis.title.y = element_text(size = 12, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 12),
		title = element_text(size = 8),
		legend.text=element_text(size=5), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.7, "lines"), # legend height
        	panel.background = element_rect(fill = "white", colour = "grey84"),
        	#plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()


### best
pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Enrichment_adjp0.05_tidy_nobkgrd_new.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = enrichment, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA, 
       #title="pfam_enrichment_abd_posFC"
       colour = "P-value",
       size = "Enrichment score"
       ) +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 2.8)) + # size of dots
	theme(axis.text.y = element_text(size = 6), # x axis (pfam) domain # changes size of plot DO NOT CHANGE!
		axis.ticks = element_blank(), 
		panel.grid.major.y = element_blank(),
		panel.grid.major.x = element_blank(),
		#panel.background = element_blank(),
		axis.text.x = element_text(size = 12, angle = 90),
		axis.title.y = element_text(size = 12, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 12),
		title = element_text(size = 4),
		legend.text=element_text(size=5), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.7, "lines"), # legend height
        	panel.background = element_rect(fill = "white", colour = "grey50"),
        	#plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()


pdf("../pfam_enrichment_dotplot/TS1015_splsda_2017.08.04_abd_pfam_clusters_2017.07.13_abd_Enrichment_adjp0.05_tidy_nobkgrd2.pdf", width=3, height=7) 
 ggplot(sig, aes(x = test, y = pfam, size = enrichment, colour = pval)) + 
 	viridis::scale_color_viridis(discrete = FALSE)+ # https://stackoverflow.com/questions/42728144/discretizing-viridis-ggplot-color-scale
 	labs(#x="Ven intersection", 
       y="PFAM Domain", 
       fill=NA,
       colour = "P-value",
       size = "Enrichment score" 
       #title="pfam_enrichment_abd_posFC"
       ) +  
 	geom_point() +
 	#scale_colour_gradient(low = "blue") +
 	#scale_x_discrete(limits = c("drk-nat",
     #"shrd-all","uniq-drk"), expand=c(0.1,0)) +
	scale_size(range = c(0.1, 2.8)) + # size of dots
	theme(axis.text.y = element_text(size = 6), # x axis (pfam) domain
		axis.ticks = element_blank(), 
		#panel.grid.major.y = element_line(colour = "grey84"),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_blank(),
		#panel.background = element_blank(),
		axis.text.x = element_blank(),
		#axis.text.x = element_text(size = 12, angle = 90),
		axis.title.y = element_text(size = 12, angle = 90), #rel(0.7) ie x
		axis.title.x = element_text(size = 12),
		#title = element_text(size = 12),
		legend.text=element_text(size=8), # legend
		legend.justification = "top",
		legend.key.width=unit(0.5, "lines"), # legend width
        	legend.key.height=unit(0.7, "lines"), # legend height
        	panel.background = element_rect(fill = "white", colour = "grey50"),
        	#plot.background=element_blank()
		) 
		#+
		#guides(colour = "none") #rm legend. http://docs.ggplot2.org/current/guides.html
dev.off()


#######################
# to change
#######################
pdf("TS1015_t9_DEGtrt_all_adjp0.1_posFC_abd_ven_2017.07.13_abd_other.pdf") 
 ggplot(sig, aes(x = test, y = pfam, size = enrichment, color = test)) +
	geom_point() +
	scale_size(range = c(0.1, 2))
dev.off()




###########################################
# to change
#
#‘sink’ diverts R output to a connection.
sink("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_h03.06_DESeq2_DEA_Design-group_t9_2017.05.15_nolfc/VEN/pfam_enrichment_l00.00_heatmap/session_info_sink_TS1015_DEGt1.t9_all_adjp0.1_posFC_2017.07.13.txt")
#sessionInfo()
#sink()
###########################################


#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")







