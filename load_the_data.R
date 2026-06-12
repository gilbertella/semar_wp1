## -------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list=ls())
library(readr)
library(dplyr)
library(tidyr)
library(phyloseq)


## -------------------------------------------------------------------------------------------------------------------------------------------------------
pats_file <- 'data/isolate_clinical_data.txt'
pat_data <- read_delim(pats_file,'\t')

# Convert percentages to proportions
tmp_data <- pat_data %>% mutate(R_AP=R_AP/100,R_FQ=R_FQ/100,R_3GC=R_3GC/100,R_AG=R_AG/100)

hospital <- tmp_data %>% filter(Source == "Hospital")
wwtp <- tmp_data %>% filter(Source == "Wastewater")

# Change names of isolate and lcinical data
hospital_resistance <- hospital %>% dplyr::select("Country","R_AP","R_FQ","R_3GC","R_AG","MAR")
colnames(hospital_resistance) <- c("Country","cl_R_AP","cl_R_FQ","cl_R_3GC","cl_R_AG","cl_R_MAR")
wwtp_resistance <- wwtp %>% dplyr::select("Country","R_AP","R_FQ","R_3GC","R_AG","MAR")
colnames(wwtp_resistance) <- c("Country","wwtp_R_AP","wwtp_R_FQ","wwtp_R_3GC","wwtp_R_AG","wwtp_R_MAR")

isolate_merged <- merge(wwtp_resistance,hospital_resistance, by = "Country")


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Paths to data files
meta_data_file <- 'data/metadata.txt'
resfinder_count_file <- 'data/semar_sample_lineage_arg_counts_90.txt'


# Load the meta data and resfinder data
meta_data <- read_delim(meta_data_file,'\t')
arg_data <- read_delim(resfinder_count_file,'\t')


# Convert the arg read counts to fragments
n_columns = dim(arg_data)[2]
n_samples <- n_columns -1
arg_counts_fragments <- ceiling(arg_data[,2:n_columns]/2)
first_column <- colnames(arg_counts_fragments)[1]
arg_data_fragments <- arg_counts_fragments %>% mutate(Gene=arg_data$Gene, .before=first_column)


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Get the counts together with lineage
lineage_to_use <- c("Class","cl_90","Family", "Gene", "ARG")
arg_data_lineage <- arg_data_fragments %>% separate(Gene,lineage_to_use,';')
rownames(arg_data_lineage) <- arg_data_lineage$ARG

# Remove the disinfectant resistance genes
arg_data_lineage <- arg_data_lineage %>% filter(!Class == "Disinfectant")

# Get only the lineages
arg_lineage <- arg_data_lineage %>% dplyr::select(all_of(lineage_to_use))


# Get only the counts
arg_counts <- arg_data_lineage %>% dplyr::select(-all_of(lineage_to_use))
rownames(arg_counts) <- arg_data_lineage$ARG


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Normalize the args by arg length (if we want to compare different genes against each other)
ARG_lengths <- read_delim("data/gene_length.txt",'\t')

# join data frame with read lengths
arg_data_tmp <- inner_join(arg_data_lineage,ARG_lengths,by='ARG')
rownames(arg_data_tmp) <- arg_data_tmp$ARG


#extract only the counts together with the arg length
arg_counts_tmp <- arg_data_tmp %>% dplyr::select(-Class,-cl_90,-Family,-Gene,-ARG)


# Normalize each cell by the corresponding arg length
n <- ncol(arg_counts_tmp)
arg_stop <- n-1
length_normalized_args <- sweep(arg_counts_tmp[,1:arg_stop], 1,
      arg_counts_tmp$length, FUN="/")
rownames(length_normalized_args) <- arg_data_lineage$ARG



## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Get total arg counts length normalized and add to meta_data frame
res_count <- apply(length_normalized_args,2,sum)

row.names(meta_data) <- meta_data$Sample
res_frame <- data.frame(res_count)
res_frame <- res_frame %>% mutate(Sample = rownames(res_frame))
new_meta_data <- merge(meta_data, res_frame, by="Sample")
rownames(new_meta_data) <- new_meta_data$Sample


# Get total arg counts and add to meta_data frame
res_count_raw <- apply(arg_counts,2,sum)
new_meta_data <- merge(new_meta_data, as.data.frame(res_count_raw), by=0, all=TRUE)
rownames(new_meta_data) <- new_meta_data$Sample
new_meta_data <- new_meta_data %>% dplyr::select(-Row.names)

# Normalize total arg counts
info_frame <- new_meta_data %>% mutate(res_count_frag_norm=res_count/Mfragments)
info_frame <- info_frame %>% mutate(res_count_bac_norm=res_count/bacterial_16s)
info_frame <- info_frame %>% mutate(intI1_count_bac_norm=intI1/bacterial_16s)

# Merge the metagenomic data and the isolate/clinical data
info_frame <- merge(info_frame,isolate_merged,by="Country")


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Create phyloseq objects
frag_norm_arg <- length_normalized_args
arg_taxa <- as.matrix(arg_lineage)
rownames(arg_taxa) <- arg_lineage$ARG
colnames(arg_taxa)


my_meta_data <- sample_data(info_frame)
sample_names(my_meta_data) <- my_meta_data$Sample

# Create one phyloseq normalized by length and total million fragments
ps_args_frag_norm <- phyloseq(otu_table(frag_norm_arg, taxa_are_rows=TRUE), tax_table(arg_taxa),my_meta_data)
otu_table(ps_args_frag_norm) <- otu_table(sweep(otu_table(ps_args_frag_norm), 2, sample_data(ps_args_frag_norm)$Mfragments, FUN="/"),taxa_are_row=TRUE)
# Convert to FPKM
otu_table(ps_args_frag_norm) <- otu_table(ps_args_frag_norm)*1000


# Create one phyloseq normalized by length and bacterial 16s
ps_args_bac_norm <- phyloseq(otu_table(frag_norm_arg, taxa_are_rows=TRUE), tax_table(arg_taxa),my_meta_data)
otu_table(ps_args_bac_norm) <- otu_table(sweep(otu_table(ps_args_bac_norm), 2, sample_data(ps_args_bac_norm)$bacterial_16s, FUN="/"),taxa_are_row=TRUE)


# Create one with just the counts
ps_args_count <- phyloseq(otu_table(arg_counts, taxa_are_rows=TRUE), tax_table(arg_taxa),my_meta_data)


## -------------------------------------------------------------------------------------------------------------------------------------------------------
semar_info_frame <- info_frame
variables_to_be_saved <- c("ps_args_count","ps_args_bac_norm","ps_args_frag_norm","semar_info_frame")
rm(list=setdiff(ls(), variables_to_be_saved))


## -------------------------------------------------------------------------------------------------------------------------------------------------------
# Define file paths
## -------------------------------------------------------------------------------------------------------------------------------------------------------
variables_to_be_saved <- c("ps_args_count","ps_args_bac_norm","ps_args_frag_norm","semar_info_frame")
rm(list=setdiff(ls(), variables_to_be_saved))

