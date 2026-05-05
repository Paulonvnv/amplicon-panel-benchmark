#where_is_the_script
# Load libraries----
source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_required_libraries.R')
source('~/Documents/Github/Plasmodium_WGS_analysis/functions_libraries/rGenome_functions.R')
source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_functions.R')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/Rcpp_functions.cpp')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/hmmloglikelihood.cpp')


if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome_allelesupdated')){
  
  # Load Genome data-----
  
  WGS_all_samples_rGenome = read_rGenome('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome', format = 'tsv', sep = '\t')
  
  sum(colnames(WGS_all_samples_rGenome@gt) != WGS_all_samples_rGenome@metadata$Sample_id)
  colnames(WGS_all_samples_rGenome@gt) = WGS_all_samples_rGenome@metadata$Sample_id
  
  # Filter out samples----
  
  dim(WGS_all_samples_rGenome@gt)
  
  WGS_all_samples_rGenome@metadata %>%
    summarise(nSamples = n(), 
              .by = c(site_of_collection_world_region, site_of_collection_snl0)) %>%
    arrange(site_of_collection_world_region, site_of_collection_snl0)
  
  WGS_all_samples_rGenome = 
    filter_samples(WGS_all_samples_rGenome,
                   (WGS_all_samples_rGenome@metadata$PvGTSeq_SampAmpRate >= .75 &
                      WGS_all_samples_rGenome@metadata$PvAmpliSeq_SampAmpRate >= .75 &
                      WGS_all_samples_rGenome@metadata$rhAmpSeq_SampAmpRate >= .75 &
                      WGS_all_samples_rGenome@metadata$WGS_SampAmpRate >= .75))
  
  
  WGS_all_samples_rGenome = filter_samples(WGS_all_samples_rGenome,
                                           !(
                                             WGS_all_samples_rGenome@metadata$batch %in% c('Duraisingh', 'Ancient') |
                                               WGS_all_samples_rGenome@metadata$site_of_collection_snl0 %in% c('P. simium', 'Unknown')
                                           ))
  
  WGS_all_samples_rGenome@metadata %>%
    summarise(nSamples = n(), 
              .by = c(site_of_collection_world_region, site_of_collection_snl0)) %>%
    arrange(site_of_collection_world_region, site_of_collection_snl0)
  
  dim(WGS_all_samples_rGenome@gt)
  
  # Calculate amplification rate per Loci-----
  WGS_all_samples_rGenome@loci_table$Locus_Ampl_Rate = LocusAmplRate(WGS_all_samples_rGenome, update = F)
  
  WGS_all_samples_rGenome@loci_table %>%
    ggplot(aes(x = Locus_Ampl_Rate)) + 
    geom_histogram()
  
  ## Remove loci with less than 75% of genome coverage----
  WGS_all_samples_rGenome = filter_loci(WGS_all_samples_rGenome, 
                                        v = WGS_all_samples_rGenome@loci_table$Locus_Ampl_Rate >= .75)
  
  dim(WGS_all_samples_rGenome@gt)
  
  # Update allele counts----
  
  names(WGS_all_samples_rGenome@loci_table)
  
  WGS_all_samples_rGenome = update_allele_lables(WGS_all_samples_rGenome, n = 500)
  
  names(WGS_all_samples_rGenome@loci_table)
  
  write_rGenome(WGS_all_samples_rGenome, 
                format = 'tsv',
                name = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome_allelesupdated',
                sep = '\t')
}else{
  
  WGS_all_samples_rGenome = read_rGenome(
    file = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome_allelesupdated',
    format = 'tsv',
    sep = '\t')
}

sum(colnames(WGS_all_samples_rGenome@gt) != WGS_all_samples_rGenome@metadata$Sample_id)
colnames(WGS_all_samples_rGenome@gt) = WGS_all_samples_rGenome@metadata$Sample_id

dim(WGS_all_samples_rGenome@gt)
dim(WGS_all_samples_rGenome@loci_table)

sum(rownames(WGS_all_samples_rGenome@gt) != rownames(WGS_all_samples_rGenome@loci_table)) 

rownames(WGS_all_samples_rGenome@gt)[1:10]

colnames(WGS_all_samples_rGenome@loci_table)

## Remove monomorphic sites----

sum(WGS_all_samples_rGenome@loci_table$Cardinality == 1)

dim(WGS_all_samples_rGenome@gt)

WGS_all_samples_rGenome = filter_loci(WGS_all_samples_rGenome, 
                                               v = WGS_all_samples_rGenome@loci_table$Cardinality > 1)

dim(WGS_all_samples_rGenome@gt)

# Calculate minor allele frequency
major_freq = sapply(WGS_all_samples_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})

WGS_all_samples_rGenome@loci_table$major_freq = major_freq

gene_description = get_gene_description_rGenome(WGS_all_samples_rGenome@loci_table, gff = '~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PlasmoDB-67_PvivaxP01.gff')

WGS_all_samples_rGenome@loci_table = cbind(WGS_all_samples_rGenome@loci_table,
                                           gene_description)

length(unique(WGS_all_samples_rGenome@loci_table$gene_id))


# 
# write_rGenome(WGS_all_samples_rGenome, format = 'tsv', name = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_broad_Pv4_WGS_rGenome', sep = '\t')
# 
# WGS_all_samples_rGenome_bachup = read_rGenome(format = 'tsv', file = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_broad_Pv4_WGS_rGenome', sep = '\t')
# 
# WGS_all_samples_rGenome = WGS_all_samples_rGenome_bachup
## Keep sites with MAF > 0.001


sum(major_freq <= 1 - 5/ncol(WGS_all_samples_rGenome@gt))


WGS_all_samples_rGenome_01 = filter_loci(WGS_all_samples_rGenome, 
                                               v = WGS_all_samples_rGenome@loci_table$major_freq <= 1 - 5/ncol(WGS_all_samples_rGenome@gt))



dim(WGS_all_samples_rGenome_01@gt)
length(unique(WGS_all_samples_rGenome_01@loci_table$gene_id))

# Identify the type of polymophism of the variant site----
WGS_all_samples_rGenome_01@loci_table$type_of_polymorphism = get_type_of_polymorphism(WGS_all_samples_rGenome_01)

unique(WGS_all_samples_rGenome_01@loci_table$type_of_polymorphism)

## Remove Homopolymers and STRs----
WGS_all_samples_rGenome_01SNPsfINDELs = 
  filter_loci(WGS_all_samples_rGenome_01,
              v = !(WGS_all_samples_rGenome_01@loci_table$type_of_polymorphism %in% 
                      c('INDEL:Homopolymer', 'INDEL:Dinucleotide_STR')))

dim(WGS_all_samples_rGenome_01SNPsfINDELs@gt)

# Calculate again sample amplification rate----
WGS_all_samples_rGenome_01SNPsfINDELs = SampleAmplRate(WGS_all_samples_rGenome_01SNPsfINDELs)

WGS_all_samples_rGenome_01SNPsfINDELs@metadata %>%
  ggplot(aes(x = SampleAmplRate)) +
  geom_histogram()

# WGS_all_samples_rGenome_01SNPsfINDELs@metadata %>%
#   summarise(nSamples75 = sum(SampleAmplRate >= .75),
#             nSamples50 = sum(SampleAmplRate >= .5),
#             .by = site_of_collection_snl0)


# Remove samples with less than 75% og coverage----
# WGS_all_samples_rGenome_01SNPsfINDELs_75 = 
#   filter_samples(WGS_all_samples_rGenome_01SNPsfINDELs, 
#                  v = WGS_all_samples_rGenome_01SNPsfINDELs@metadata$SampleAmplRate >= .75) # This must be 0.75
# 
# nrow(WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata)

WGS_all_samples_rGenome_01SNPsfINDELs_75 = WGS_all_samples_rGenome_01SNPsfINDELs

dim(WGS_all_samples_rGenome_01SNPsfINDELs_75@gt)

# update allele frequencies ----

WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Allele_Counts = NULL
WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Alleles = NULL
WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Cardinality = NULL
WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$major_freq = NULL

names(WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table)

WGS_all_samples_rGenome_01SNPsfINDELs_75 = update_allele_lables(WGS_all_samples_rGenome_01SNPsfINDELs_75)

names(WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table)
dim(WGS_all_samples_rGenome_01SNPsfINDELs_75@gt)

## Remove monomorphic sites----

sum(WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Cardinality == 1)

# WGS_all_samples_rGenome_01SNPsfINDELs = filter_loci(WGS_all_samples_rGenome_01SNPsfINDELs, 
#                                                v = WGS_all_samples_rGenome_01SNPsfINDELs@loci_table$Cardinality > 1)
# 
# nrow(WGS_all_samples_rGenome_01SNPsfINDELs@gt)

# Calculate minor allele frequency
major_freq = sapply(WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})

WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$major_freq = major_freq

WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table %>%
  ggplot(aes(x = major_freq)) + 
  geom_histogram()



sum(WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$major_freq <= 1 - 5/ncol(WGS_all_samples_rGenome_01SNPsfINDELs_75@gt))

# WGS_all_samples_rGenome_01SNPsfINDELs_75 = filter_loci(WGS_all_samples_rGenome_01SNPsfINDELs_75, 
#                                                   v = WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$major_freq <= 1 - 5/ncol(WGS_all_samples_rGenome_01SNPsfINDELs_75@gt))

dim(WGS_all_samples_rGenome_01SNPsfINDELs_75@gt)

# Keep only biallelic SNPS----

WGS_all_samples_rGenome_01SNPsbiallelic_75 = 
  filter_loci(WGS_all_samples_rGenome_01SNPsfINDELs_75,
              v = (WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$type_of_polymorphism == 'SNP' &
                     WGS_all_samples_rGenome_01SNPsfINDELs_75@loci_table$Cardinality == 2))


dim(WGS_all_samples_rGenome_01SNPsbiallelic_75@loci_table)



# Calculate observed heterozygousity and within host divergence----
WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$ObsHet = get_ObsHet(WGS_all_samples_rGenome_01SNPsbiallelic_75, by = 'sample')
WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$Fws = get_Fws_rGenome(obj = WGS_all_samples_rGenome_01SNPsbiallelic_75)

WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata %>%
  ggplot(aes(y = Fws, x = ObsHet)) +
  geom_point(alpha = .3) + 
  geom_hline(yintercept = .975)+
  geom_vline(xintercept = .005)

# Define monoclonal and polyclonal infections----
Monoclonals = WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata %>%
  filter(Fws >= .975, ObsHet <= 0.005) %>%
  select(Sample_id) %>% unlist()

Polyclonals = WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata %>%
  filter(!(Sample_id %in% Monoclonals)) %>%
  select(Sample_id) %>% unlist()

WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata %<>%
  mutate(Clonality = case_when(
    Sample_id %in% Monoclonals ~ 'Monoclonal',
    .default = 'Polyclonal'
  ))


nrow(WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata)
ncol(WGS_all_samples_rGenome_01SNPsbiallelic_75@gt)


WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata$ObsHet = 
  WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$ObsHet

WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata$Fws = 
  WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$Fws

WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata$Clonality = 
  WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$Clonality


# Generate table of all pairwise comparisons----

all_pairs_full = as.data.frame(t(combn(WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata$Sample_id, 2)))

names(all_pairs_full) = c('Yi', 'Yj')

all_pairs_full = left_join(
  all_pairs_full,
  WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(all_pairs_full) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


all_pairs_full = left_join(
  all_pairs_full,
  WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(all_pairs_full) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')


all_pairs_by_country = all_pairs_full %>% 
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0)


# Select monoclonal infections----


Monoclonals = 
  WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata %>%
  filter(Clonality == 'Monoclonal') %>%
  select(Sample_id) %>%
  unlist()

monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsfINDELs_75 = 
  filter_samples(WGS_all_samples_rGenome_01SNPsfINDELs_75,
                 v = WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata$Clonality == 'Monoclonal')

monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75 = 
  filter_samples(WGS_all_samples_rGenome_01SNPsbiallelic_75,
                 v = WGS_all_samples_rGenome_01SNPsbiallelic_75@metadata$Clonality == 'Monoclonal')


write.table(monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$Sample_id,
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/monoclonals_Pv_all_samples.bed',
            sep = '\t',
            quote = F,
            row.names = F,
            col.names = F
            )

# monoclonals_PvBroad_SelectedPos2_rGenome = 
#   SampleAmplRate(obj = monoclonals_PvBroad_SelectedPos2_rGenome)
# 
# monoclonals_PvBroad_SelectedPos2_rGenome = 
#   filter_samples(monoclonals_PvBroad_SelectedPos2_rGenome,
#                  v = monoclonals_PvBroad_SelectedPos2_rGenome@metadata$SampleAmplRate >= .5)


# Monoclonal pairs by country----

Monoclonal_pairs = as.data.frame(t(combn(Monoclonals, 2)))

names(Monoclonal_pairs) = c('Yi', 'Yj')

Monoclonal_pairs = left_join(
  Monoclonal_pairs,
  WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


Monoclonal_pairs = left_join(
  Monoclonal_pairs,
  WGS_all_samples_rGenome_01SNPsfINDELs_75@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')

# Run pairwise IBS on monoclonal infections----

monoclonals_gt_01SNPsbiallelic_75 = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@gt))

monoclonals_gt_01SNPsbiallelic_75 = t(monoclonals_gt_01SNPsbiallelic_75)


monoclonals_gt_01SNPsfINDELs_75 = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsfINDELs_75@gt))

monoclonals_gt_01SNPsfINDELs_75 = t(monoclonals_gt_01SNPsfINDELs_75)

monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci = create_loci(loci_table = monoclonals_gt_01SNPsbiallelic_75,
                                                                     metadata = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata)

dim(monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@loci_table)


monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci = create_loci(loci_table = monoclonals_gt_01SNPsfINDELs_75,
                                                    metadata = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsfINDELs_75@metadata)

dim(monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci@loci_table)


monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci@metadata %>%
  summarise(nSamples = n(), 
            .by = c(site_of_collection_world_region, site_of_collection_snl0)) %>%
  arrange(site_of_collection_world_region, site_of_collection_snl0)

sum(rownames(monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci@loci_table) != monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci@metadata$Sample_id)


if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsbiallelic_75_dist.csv')){
  
  Monoclonal_pairs_01SNPsbiallelic_75_dist = NULL
  for(w in 1:500){
    start_time = Sys.time()
    Monoclonal_pairs_01SNPsbiallelic_75_dist = 
      rbind(Monoclonal_pairs_01SNPsbiallelic_75_dist,
            pairwise_euclidean(
              obj = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci, parallel = T, pairs = Monoclonal_pairs,
              w = w, n = 500
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_pairs_01SNPsbiallelic_75_dist %<>%
    mutate(rhat = 1 - euDist)
  
  Monoclonal_pairs_01SNPsbiallelic_75_dist = 
    left_join(Monoclonal_pairs_01SNPsbiallelic_75_dist,
              Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #sum(Monoclonal_pairs_01SNPsbiallelic_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_pairs_01SNPsbiallelic_75_dist, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsbiallelic_75_dist.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_pairs_01SNPsbiallelic_75_dist = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsbiallelic_75_dist.csv')
  
}


Monoclonal_pairs_01SNPsbiallelic_75_dist %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_pairs_01SNPsbiallelic_75_dist %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')




if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsfINDELs_75_dist.csv')){
  
  Monoclonal_pairs_01SNPsfINDELs_75_dist = NULL
  for(w in 1:500){
    start_time = Sys.time()
    Monoclonal_pairs_01SNPsfINDELs_75_dist = 
      rbind(Monoclonal_pairs_01SNPsfINDELs_75_dist,
            pairwise_euclidean(
              obj = monoclonals_PvBroad_SelectedPos2_01SNPsfINDELs_75_loci, parallel = T, pairs = Monoclonal_pairs,
              w = w, n = 500
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_pairs_01SNPsfINDELs_75_dist %<>%
    mutate(rhat = 1 - euDist)
  
  Monoclonal_pairs_01SNPsfINDELs_75_dist = 
    left_join(Monoclonal_pairs_01SNPsfINDELs_75_dist,
              Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #names(Monoclonal_pairs_01SNPsfINDELs_75_dist) = c(names(Monoclonal_pairs_01SNPsfINDELs_75_dist)[-ncol(Monoclonal_pairs_01SNPsfINDELs_75_dist)], 'site_of_collection_snl0')
  
  
  #sum(Monoclonal_pairs_01SNPsfINDELs_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_pairs_01SNPsfINDELs_75_dist, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsfINDELs_75_dist.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_pairs_01SNPsfINDELs_75_dist = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_01SNPsfINDELs_75_dist.csv')
  
}

Monoclonal_pairs_01SNPsfINDELs_75_dist$site_of_collection_snl0 = NULL

names(Monoclonal_pairs_01SNPsfINDELs_75_dist)

unique(Monoclonal_pairs_01SNPsfINDELs_75_dist$site_of_collection_snl0)

Monoclonal_pairs_01SNPsfINDELs_75_dist %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_pairs_01SNPsfINDELs_75_dist %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')


plot(x = Monoclonal_pairs_01SNPsfINDELs_75_dist$rhat, y = Monoclonal_pairs_01SNPsbiallelic_75_dist$rhat)


# Identify clonal groups in WGS----

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_clusters_wgs.tsv')){
  thres = 1
  
  temp_cluster = plot_network(pairwise_relatedness = Monoclonal_pairs_01SNPsbiallelic_75_dist,  
                              threshold = thres,
                              metadata = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata,
                              sample_id = 'Sample_id',
                              levels = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata %>%
                                select(site_of_collection_world_region) %>% unlist() %>% unique(),
                              group_by = 'site_of_collection_world_region',
                              colors = brewer.pal(11, 'Set3'),
                              vertex.size = 4,
                              method = 'fruchtermanreingold')
  
  
  summary_clusters = temp_cluster$clusters %>% 
    summarise(Total = n(),
              nsampClus = sum(grepl('Cluster', Cluster)),
              nSampSing = n() - sum(grepl('Cluster', Cluster)),
              nclusters = length(unique(.$Cluster[grepl('Cluster', Cluster)]))
    )
  
  summary_clusters['Threshold'] = thres
  
  temp_summary_clusters = summary_clusters
  
  while(temp_summary_clusters['nclusters'] > 1# & thres > 0.60
  ){
    thres = thres - 0.01
    
    temp_cluster = plot_network(pairwise_relatedness = Monoclonal_pairs_01SNPsbiallelic_75_dist,  
                                threshold = thres,
                                metadata = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata,
                                sample_id = 'Sample_id',
                                levels = monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata %>%
                                  select(site_of_collection_world_region) %>% unlist() %>% unique(),
                                group_by = 'site_of_collection_world_region',
                                colors = brewer.pal(11, 'Set3'),
                                vertex.size = 4,
                                method = 'fruchtermanreingold')
    
    temp_summary_clusters = temp_cluster$clusters %>% 
      summarise(Total = n(),
                nsampClus = sum(grepl('Cluster', Cluster)),
                nSampSing = n() - sum(grepl('Cluster', Cluster)),
                nclusters = length(unique(.$Cluster[grepl('Cluster', Cluster)])))
    
    temp_summary_clusters['Threshold'] = thres
    
    summary_clusters = rbind(summary_clusters,
                             temp_summary_clusters)
    
    print(thres)
    
  }
  
  write.table(summary_clusters, '~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_clusters_wgs.tsv', quote = F, row.names = F, sep = '\t')
  
}else{
  summary_clusters = read.table('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_clusters_wgs.tsv', sep = '\t', header = T)  
}


summary_clusters$deriv_nsampClus = summary_clusters$nsampClus
summary_clusters$deriv_nClus = summary_clusters$nclusters


summary_clusters$nclusters
summary_clusters$nSampSing
summary_clusters$nsampClus

summary_clusters[-1,]$deriv_nsampClus = summary_clusters[-1,]$nsampClus - summary_clusters[-nrow(summary_clusters),]$nsampClus
summary_clusters[-1,]$deriv_nClus = summary_clusters[-1,]$nclusters - summary_clusters[-nrow(summary_clusters),]$nclusters

summary_clusters[-nrow(summary_clusters),]$deriv_nClus = summary_clusters[-nrow(summary_clusters),]$nclusters - summary_clusters[-1,]$nclusters


summary_clusters %>% 
  # mutate(deriv_nsampClus = case_when(
  #   deriv_nsampClus < 100 ~ deriv_nsampClus,
  #   .default = NA
  # ))%>%
  pivot_longer(cols = c(nclusters,
                        deriv_nClus,
                        nsampClus,
                        nSampSing,
                        deriv_nsampClus
  ), 
  names_to = 'Metric',
  values_to = 'Value') %>%
  ggplot(aes(x = Threshold, y = Value)) + 
  geom_line() + 
  facet_grid(Metric~., scales = 'free')

summary_clusters %>% 
  ggplot(aes(x = Threshold, y = nclusters)) + 
  geom_line(linewidth = 2)


summary_clusters %>% 
  ggplot(aes(x = Threshold, y = deriv_nClus)) + 
  geom_line(linewidth = 2) +
  scale_x_continuous(breaks = seq(.87, 1, .01)) + 
  theme(axis.text.x = element_text(angle = 90))



full_set_network_098  = plot_ggnetwork(pairwise_relatedness = Monoclonal_pairs_01SNPsbiallelic_75_dist,
                                       threshold = .97, 
                                       metadata = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata,
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_world_region', vertex.size = 4)

full_set_network_098$plot_network


cluster_wgsIBS098 = plot_network(pairwise_relatedness = Monoclonal_pairs_01SNPsbiallelic_75_dist,  
                                 threshold = .99,
                                 metadata = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata,
                                 sample_id = 'Sample_id',
                                 levels = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata %>%
                                   select(site_of_collection_world_region) %>% unlist() %>% unique(),
                                 group_by = 'site_of_collection_world_region',
                                 colors = brewer.pal(11, 'Set3'),
                                 vertex.size = 4,
                                 method = 'fruchtermanreingold')

selected_samples_ids = cluster_wgsIBS098$clusters %>%
  filter(!grepl('Cluster', Cluster)) %>%
  select(Sample_id) %>% unlist()

cluster_list = cluster_wgsIBS098$clusters %>%
  filter(grepl('Cluster', Cluster)) %>%
  select(Cluster) %>% unlist() %>% unique()

for(Cluster_of_interest in cluster_list){
  samples_in_cluster = cluster_wgsIBS098$clusters %>%
    filter(Cluster == Cluster_of_interest) %>%
    select(Sample_id) %>% unlist()
  
  set.seed(1000)
  selected_samples_ids = c(selected_samples_ids,
                           sample(samples_in_cluster, 3, replace = F))
}

length(selected_samples_ids)

write.table(selected_samples_ids, '~/Documents/Github/PvGTSeq_paper_draft/Outputs/selected_sample_ids.bed', sep = '\n', quote = F, row.names = F, col.names = F)


# WGS PCoA----

evectors_wgs_world_regions_withoutfilter = GRM_evectors(dist_table = Monoclonal_pairs_01SNPsbiallelic_75_dist, 
                                                        k = nrow(monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata %>%
                                                                   filter(Sample_id %in% selected_samples_ids,
                                                                          possible_origin != 'Imported'
                                                                   )
                                                        ),
                                                        metadata = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata%>%
                                                          filter(Sample_id %in% selected_samples_ids,
                                                                 possible_origin != 'Imported'),
                                                        Pop = 'site_of_collection_world_region',
                                                        fastGRM = F,
                                                        method = 'princomp'
)

evectors_wgs_world_regions_withoutfilter$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = site_of_collection_world_region
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('A) World regions using\n',
                      ncol(monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@loci_table),
                      ' SNPs from WGS'
  ),
  x = paste0('PCoA 1 (', round(evectors_wgs_world_regions_withoutfilter$contrib[1], 2), '%)'),
  y = paste0('PCoA 2 (', round(evectors_wgs_world_regions_withoutfilter$contrib[2], 2), '%)'),
  color = 'Regions',
  shape = 'Regions'
  ) +
  theme(#legend.position = 'inside', legend.position.inside = c(.8, .8),
    legend.position = 'bottom',
    #legend.justification = 'left',
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    title = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 4,
                              theme = theme(legend.title.position = 'left'
                              )))


AF_outliers = evectors_wgs_world_regions_withoutfilter$eigenvector %>%
  filter(PC2 < .2 & site_of_collection_world_region == 'AF') %>% select(Sample_id) %>% unlist()

Asia_outliers = evectors_wgs_world_regions_withoutfilter$eigenvector %>%
  filter(PC2 > .1 & site_of_collection_world_region %in% c('WSEA', 'MSEA', 'ESEA', 'OCE')) %>% select(Sample_id) %>% unlist()

LAC_outliers = evectors_wgs_world_regions_withoutfilter$eigenvector %>%
  filter(PC1 < .5 & site_of_collection_world_region %in% c('LAC')) %>% select(Sample_id) %>% unlist()





evectors_wgs_world_regions = GRM_evectors(dist_table = Monoclonal_pairs_01SNPsbiallelic_75_dist, 
                                          k = nrow(monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata %>%
                                                     filter(Sample_id %in% selected_samples_ids,
                                                            !(Sample_id %in% AF_outliers),
                                                            !(Sample_id %in% Asia_outliers),
                                                            !(Sample_id %in% LAC_outliers),
                                                            possible_origin != 'Imported'
                                                     )
                                          ),
                                          metadata = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata%>%
                                            filter(Sample_id %in% selected_samples_ids,
                                                   !(Sample_id %in% AF_outliers),
                                                   !(Sample_id %in% Asia_outliers),
                                                   !(Sample_id %in% LAC_outliers),
                                                   possible_origin != 'Imported'),
                                          Pop = 'site_of_collection_world_region',
                                          fastGRM = F,
                                          method = 'princomp')


PCoA_no_clonal_wgs_world_regions = evectors_wgs_world_regions$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = PC2, 
             color = site_of_collection_world_region
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('A) World regions using\n',
                      ncol(monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@loci_table),
                      ' biallelic SNPs from WGS'
  ),
    x = paste0('PCoA 1 (', round(evectors_wgs_world_regions$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_wgs_world_regions$contrib[2], 2), '%)'),
    color = 'Regions',
    shape = 'Regions'
  ) +
  theme(legend.position = 'inside', legend.position.inside = c(.8, .8),
    #legend.position = 'bottom',
    legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12))# +
  # guides(color = guide_legend(ncol = 4,
  #                             theme = theme(legend.title.position = 'left'
  #                             )))

PCoA_no_clonal_wgs_world_regions

# LAC PCoA----

monoclonals_LAC_rGenome_01SNPsbiallelic_75 = filter_samples(monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75,
                                                            monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$site_of_collection_world_region == 'LAC')


sum(colnames(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt) != monoclonals_LAC_rGenome_01SNPsbiallelic_75@metadata$Sample_id)

monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Allele_Counts = NULL
monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Alleles = NULL
monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Cardinality = NULL
monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$major_freq = NULL

names(monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table)

monoclonals_LAC_rGenome_01SNPsbiallelic_75 = update_allele_lables(monoclonals_LAC_rGenome_01SNPsbiallelic_75)

names(monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table)


## Remove monomorphic sites----

sum(monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Cardinality == 1)

monoclonals_LAC_rGenome_01SNPsbiallelic_75 = filter_loci(monoclonals_LAC_rGenome_01SNPsbiallelic_75,
                                                         v = monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Cardinality > 1)

nrow(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt)

# Calculate minor allele frequency
major_freq = sapply(monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})

monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$major_freq = major_freq

monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table %>%
  ggplot(aes(x = major_freq)) + 
  geom_histogram()

dim(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt)

sum(monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$major_freq <= 1 - 5/ncol(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt))

monoclonals_LAC_rGenome_01SNPsbiallelic_75 = filter_loci(monoclonals_LAC_rGenome_01SNPsbiallelic_75, 
                                                         v = monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table$major_freq <= 1 - 5/ncol(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt))


monoclonals_LAC_rGenome_01SNPsbiallelic_75@metadata %>%
  #filter(Sample_id %in% selected_samples_ids) %>%
  summarise(nSamples = n(),
            .by = c(site_of_collection_snl0))

dim(monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt)

monoclonals_LAC_rGenome_01SNPsbiallelic_75_non_clonal =
  filter_samples(monoclonals_LAC_rGenome_01SNPsbiallelic_75,
                 monoclonals_LAC_rGenome_01SNPsbiallelic_75@metadata$Sample_id %in% selected_samples_ids)


dim(monoclonals_LAC_rGenome_01SNPsbiallelic_75_non_clonal@gt)

monoclonals_LAC_rGenome_01SNPsbiallelic_75@metadata %>%
  filter(Sample_id %in% selected_samples_ids) %>%
  summarise(nSamples = n(),
            .by = c(site_of_collection_snl0))

monoclonals_LAC_01SNPsbiallelic_75_gt = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt))

monoclonals_LAC_01SNPsbiallelic_75_gt = gsub(':\\d+', '', monoclonals_LAC_rGenome_01SNPsbiallelic_75@gt)

monoclonals_LAC_01SNPsbiallelic_75_gt = t(monoclonals_LAC_01SNPsbiallelic_75_gt)

monoclonals_LAC_01SNPsbiallelic_75_loci = create_loci(loci_table = monoclonals_LAC_01SNPsbiallelic_75_gt, 
                                                      metadata = monoclonals_LAC_rGenome_01SNPsbiallelic_75@metadata)


allele_freqs = get_allele_freq(monoclonals_LAC_01SNPsbiallelic_75_loci)

monoclonals_LAC_01SNPsbiallelic_75_loci@freq_table = allele_freqs

nrow(allele_freqs)
ncol(monoclonals_LAC_01SNPsbiallelic_75_loci@loci_table)


monoclonals_LAC_01SNPsbiallelic_75_loci@markers = 
  monoclonals_LAC_rGenome_01SNPsbiallelic_75@loci_table %>% S4Vectors::rename('CHROM' = 'chromosome',
                                                                   'POS' = 'pos')



Monoclonal_LAC_pairs = as.data.frame(t(combn(monoclonals_LAC_01SNPsbiallelic_75_loci@metadata$Sample_id, 2)))

names(Monoclonal_LAC_pairs) = c('Yi', 'Yj')

Monoclonal_LAC_pairs = left_join(
  Monoclonal_LAC_pairs,
  monoclonals_LAC_01SNPsbiallelic_75_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(Monoclonal_LAC_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


Monoclonal_LAC_pairs = left_join(
  Monoclonal_LAC_pairs,
  monoclonals_LAC_01SNPsbiallelic_75_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(Monoclonal_LAC_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')

nrow(Monoclonal_LAC_pairs)/5000

sum(rownames(monoclonals_LAC_01SNPsbiallelic_75_loci@loci_table) != monoclonals_LAC_01SNPsbiallelic_75_loci@metadata$Sample_id)



if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:5000){
    start_time = Sys.time()
    Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_LAC_01SNPsbiallelic_75_loci, parallel = T, pairs = Monoclonal_LAC_pairs, max_k = 20,
              w = w, n = 5000
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %<>%
    mutate(euDist = 1 - rhat)
  
  Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD = 
    left_join(Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,
              Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #sum(Monoclonal_LAC_pairs_01SNPsbiallelic_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD.csv')
  
}


max(Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD$rhat)

Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')


Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == 'Panama', Yj_site_of_collection_snl0 == 'Panama')

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_wgs_LAC.tsv')){
  thres = 1
  
  temp_cluster = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,  
                              threshold = thres,
                              metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata,
                              sample_id = 'Sample_id',
                              levels = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                select(site_of_collection_snl0) %>% unlist() %>% unique(),
                              group_by = 'site_of_collection_snl0',
                              colors = brewer.pal(11, 'Set3'),
                              vertex.size = 4,
                              method = 'fruchtermanreingold')
  
  
  summary_IBDclusters_wgs_LAC = temp_cluster$clusters %>% 
    summarise(Total = n(),
              nsampClus = sum(grepl('Cluster', Cluster)),
              nSampSing = n() - sum(grepl('Cluster', Cluster)),
              nclusters = length(unique(.$Cluster[grepl('Cluster', Cluster)]))
    )
  
  summary_IBDclusters_wgs_LAC['Threshold'] = thres
  
  temp_summary_clusters = summary_IBDclusters_wgs_LAC
  
  while(temp_summary_clusters['nclusters'] > 1 | thres == 1
  ){
    thres = thres - 0.01
    
    temp_cluster = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,  
                                threshold = thres,
                                metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata,
                                sample_id = 'Sample_id',
                                levels = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                  select(site_of_collection_snl0) %>% unlist() %>% unique(),
                                group_by = 'site_of_collection_snl0',
                                colors = brewer.pal(11, 'Set3'),
                                vertex.size = 4,
                                method = 'fruchtermanreingold')
    
    temp_summary_clusters = temp_cluster$clusters %>% 
      summarise(Total = n(),
                nsampClus = sum(grepl('Cluster', Cluster)),
                nSampSing = n() - sum(grepl('Cluster', Cluster)),
                nclusters = length(unique(.$Cluster[grepl('Cluster', Cluster)])))
    
    temp_summary_clusters['Threshold'] = thres
    
    summary_IBDclusters_wgs_LAC = rbind(summary_IBDclusters_wgs_LAC,
                                        temp_summary_clusters)
    
    print(thres)
    
  }
  
  write.table(summary_IBDclusters_wgs_LAC, '~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_wgs_LAC.tsv', quote = F, row.names = F, sep = '\t')
  
}else{
  summary_IBDclusters_wgs_LAC = read.table('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_wgs_LAC.tsv', sep = '\t', header = T)  
}



summary_IBDclusters_wgs_LAC$deriv_nsampClus = summary_IBDclusters_wgs_LAC$nsampClus
summary_IBDclusters_wgs_LAC$deriv_nClus = summary_IBDclusters_wgs_LAC$nclusters


summary_IBDclusters_wgs_LAC$nclusters
summary_IBDclusters_wgs_LAC$nSampSing
summary_IBDclusters_wgs_LAC$nsampClus
#summary_clusters$deriv = NA

summary_IBDclusters_wgs_LAC[-1,]$deriv_nsampClus = summary_IBDclusters_wgs_LAC[-1,]$nsampClus - summary_IBDclusters_wgs_LAC[-nrow(summary_IBDclusters_wgs_LAC),]$nsampClus
summary_IBDclusters_wgs_LAC[-1,]$deriv_nClus = summary_IBDclusters_wgs_LAC[-1,]$nclusters - summary_IBDclusters_wgs_LAC[-nrow(summary_IBDclusters_wgs_LAC),]$nclusters

summary_IBDclusters_wgs_LAC[-nrow(summary_IBDclusters_wgs_LAC),]$deriv_nClus = summary_IBDclusters_wgs_LAC[-nrow(summary_IBDclusters_wgs_LAC),]$nclusters - summary_IBDclusters_wgs_LAC[-1,]$nclusters


summary_IBDclusters_wgs_LAC %>% 
  mutate(deriv_nsampClus = case_when(
    deriv_nsampClus < 100 ~ deriv_nsampClus,
    .default = NA
  ))%>%
  pivot_longer(cols = c(nclusters,
                        deriv_nClus,
                        nsampClus,
                        nSampSing,
                        deriv_nsampClus
  ), 
  names_to = 'Metric',
  values_to = 'Value') %>%
  ggplot(aes(x = Threshold, y = Value)) + 
  geom_line() + 
  facet_grid(Metric~., scales = 'free')

summary_IBDclusters_wgs_LAC %>% 
  ggplot(aes(x = Threshold, y = nclusters)) + 
  geom_line(linewidth = 2)


summary_IBDclusters_wgs_LAC %>% 
  ggplot(aes(x = Threshold, y = deriv_nClus)) + 
  geom_line(linewidth = 2) +
  scale_x_continuous(breaks = seq(0, 1, .02)) + 
  theme(axis.text.x = element_text(angle = 90))


full_set_network_095  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,
                                       threshold = .95, 
                                       metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                         filter(!(Sample_id %in% LAC_outliers),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_snl0', vertex.size = 4)

full_set_network_095$plot_network

cluster_wgsIBD095 = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,  
                                 threshold = .95,
                                 metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata,
                                 sample_id = 'Sample_id',
                                 levels = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                   select(site_of_collection_snl0) %>% unlist() %>% unique(),
                                 group_by = 'site_of_collection_snl0',
                                 colors = brewer.pal(11, 'Set3'),
                                 vertex.size = 4,
                                 method = 'fruchtermanreingold')


selected_LAC_samples_ids = cluster_wgsIBD095$clusters %>%
  filter(!grepl('Cluster', Cluster)) %>%
  select(Sample_id) %>% unlist()

LAC_cluster_list = cluster_wgsIBD095$clusters %>%
  filter(grepl('Cluster', Cluster)) %>%
  select(Cluster) %>% unlist() %>% unique()

for(Cluster_of_interest in LAC_cluster_list){
  samples_in_cluster = cluster_wgsIBD095$clusters %>%
    filter(Cluster == Cluster_of_interest) %>%
    select(Sample_id) %>% unlist()
  
  set.seed(1000)
  selected_LAC_samples_ids = c(selected_LAC_samples_ids,
                               sample(samples_in_cluster, 3, replace = F))
}


full_set_network_048  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,
                                       threshold = 0.48, 
                                       metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                         filter(!(Sample_id %in% LAC_outliers),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_snl0', vertex.size = 4)

full_set_network_048$plot_network

full_set_network_026  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD,
                                       threshold = 0.26, 
                                       metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                         filter(!(Sample_id %in% LAC_outliers),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_snl0', vertex.size = 4)

full_set_network_026$plot_network





IBD_evectors_wgs_LAC_withoutfilters = GRM_evectors(dist_table = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD, 
                                    k = nrow(monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                               filter(Sample_id %in% selected_samples_ids,
                                                      !(Sample_id %in% LAC_outliers),
                                                      !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana')))
                                    ),
                                    metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata%>%
                                      filter(Sample_id %in% selected_samples_ids,
                                             !(Sample_id %in% LAC_outliers),
                                             !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                    Pop = 'site_of_collection_snl0',
                                    fastGRM = F,
                                    cor = F,
                                    method = 'princomp') # fastSVD 


IBD_evectors_wgs_LAC_withoutfilters$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = -PC2, 
             color = site_of_collection_snl0
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('A) LAC countries using\n',
                      ncol(monoclonals_LAC_01SNPsbiallelic_75_loci@loci_table),
                      ' SNPs from WGS'
  ),
    x = paste0('PCoA 1 (', round(IBD_evectors_wgs_LAC_withoutfilters$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(IBD_evectors_wgs_LAC_withoutfilters$contrib[2], 2), '%)'),
    color = 'Study sites',
    shape = 'Study sites'
  ) +
  theme(#legend.position = 'inside', legend.position.inside = c(.7, .5),
    legend.position = 'bottom',    
    #legend.title = element_blank(),
        legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 4,
                              theme = theme(legend.title.position = 'left'
                              )))


Brazil_outliers = IBD_evectors_wgs_LAC_withoutfilters$eigenvector %>%
  filter(PC2 > 1) %>% select(Sample_id) %>% unlist()



IBD_evectors_wgs_LAC = GRM_evectors(dist_table = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD, 
                                    k = nrow(monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                               filter(Sample_id %in% selected_samples_ids,
                                                      !(Sample_id %in% LAC_outliers),
                                                      #!(Sample_id %in% Brazil_outliers),
                                                      !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana')))
                                    ),
                                    metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata%>%
                                      filter(Sample_id %in% selected_samples_ids,
                                             !(Sample_id %in% LAC_outliers),
                                             #!(Sample_id %in% Brazil_outliers),
                                             !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                    Pop = 'site_of_collection_snl0',
                                    fastGRM = F,
                                    cor = F,
                                    method = 'princomp') # fastSVD 


PCoA_no_clonal_wgs_LAC = IBD_evectors_wgs_LAC$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = -PC2, 
             color = site_of_collection_snl0
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('A) LAC countries using\n',
                      ncol(monoclonals_LAC_01SNPsbiallelic_75_loci@loci_table),
                      ' SNPs from WGS'
  ),
  x = paste0('PCoA 1 (', round(IBD_evectors_wgs_LAC$contrib[1], 2), '%)'),
  y = paste0('PCoA 2 (', round(IBD_evectors_wgs_LAC$contrib[2], 2), '%)'),
  color = 'Country',
  shape = 'Country'
  ) +
  theme(legend.position = 'inside', legend.position.inside = c(.4, .05),
    #legend.position = 'bottom',    
    legend.title = element_blank(),
    legend.justification = 'left',
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    title = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 3,
                              theme = theme(legend.title.position = 'left'
                              )))

PCoA_no_clonal_wgs_LAC





ggdraw() + 
  draw_plot(PCoA_no_clonal_wgs_world_regions + 
              theme(legend.position = 'inside', 
                    legend.title = element_blank(),
                    legend.position.inside = c(.7, .7)
                    ),
            x = 0,
            width = .5,
            y = .5, 
            height = .5) +
  draw_plot(PCoA_no_clonal_PvGTSeq_world_regions + theme(legend.position = 'none'),
            x = 0.5,
            width = .5,
            y = .5, 
            height = .5) +
  draw_plot(PCoA_no_clonal_PvAmpliSeq_world_regions + theme(legend.position = 'none'),
            x = 0,
            width = .5,
            y = .0, 
            height = .5) +
  draw_plot(PCoA_no_clonal_rhAmpSeq_world_regions + theme(legend.position = 'none'),
            x = 0.5,
            width = .5,
            y = 0, 
            height = .5)



ggdraw() + 
  draw_plot(PCoA_no_clonal_wgs_LAC + 
              scale_y_continuous(limits = c(-.75, 1))+
              theme(legend.position = 'inside', 
                    legend.title = element_blank(),
                    legend.position.inside = c(.15, .9)
              ),
            x = 0,
            width = .5,
            y = .5, 
            height = .5) +
  draw_plot(PCoA_no_clonal_PvGTSeq_LAC + theme(legend.position = 'none'),
            x = 0.5,
            width = .5,
            y = .5, 
            height = .5) +
  draw_plot(PCoA_no_clonal_PvAmpliSeq_LAC + theme(legend.position = 'none'),
            x = 0,
            width = .5,
            y = .0, 
            height = .5) +
  draw_plot(PCoA_no_clonal_rhAmpSeq_LAC + theme(legend.position = 'none'),
            x = 0.5,
            width = .5,
            y = 0, 
            height = .5)





LAC_IBD_sample_ids = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata%>%
  filter(Sample_id %in% selected_samples_ids,
         !(Sample_id %in% LAC_outliers),
         #!(Sample_id %in% Brazil_outliers),
         !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))) %>%
  select(Sample_id) %>% unlist()




IBD_values_LAC = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %>%
  filter(Yi %in% LAC_IBD_sample_ids, Yj %in% LAC_IBD_sample_ids) %>%
  S4Vectors::rename('rhat' = 'IBDwgs')


IBD_values_LAC = left_join(IBD_values_LAC,
                           Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %>%
                             filter(Yi %in% LAC_IBD_sample_ids, Yj %in% LAC_IBD_sample_ids) %>%
                             S4Vectors::rename('rhat' = 'PvGTSeq') %>%
                             select(Yi, Yj, PvGTSeq),
                           by = join_by(Yi, Yj)
                           )


IBD_values_LAC = left_join(IBD_values_LAC,
                           Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %>%
                             filter(Yi %in% LAC_IBD_sample_ids, Yj %in% LAC_IBD_sample_ids) %>%
                             S4Vectors::rename('rhat' = 'PvAmpliSeq') %>%
                             select(Yi, Yj, PvAmpliSeq),
                           by = join_by(Yi, Yj)
)


IBD_values_LAC = left_join(IBD_values_LAC,
                           Monoclonal_LAC_pairs_rhAmpSeq_01SNPsbiallelic_75_IBD %>%
                             filter(Yi %in% LAC_IBD_sample_ids, Yj %in% LAC_IBD_sample_ids) %>%
                             S4Vectors::rename('rhat' = 'rhAmpSeq') %>%
                             select(Yi, Yj, rhAmpSeq),
                           by = join_by(Yi, Yj)
)


IBD_corr_plot = IBD_values_LAC %>%
  pivot_longer(cols = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'),
               names_to = 'Method',
               values_to = 'IBD') %>%
  mutate(Method = factor(Method, levels = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'))) %>%
  ggplot(aes(x = IBDwgs, y = IBD, color = Method))+
  geom_point(alpha = .1) +
  geom_line(data = data.frame(IBDwgs = rep(seq(0, 1, .01), times = 3),
                              IBD = rep(seq(0, 1, .01), times = 3),
                              Method = factor(rep(c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'), each = 101), 
                                              levels = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'))
                              ),
            aes(x = IBDwgs, y = IBD), color = 'firebrick3') +
  scale_color_manual(values = c('gray10', 'firebrick3', 'dodgerblue3')) + 
  scale_x_continuous(breaks = seq(0, 1, .25), labels = c('0', '0.25', '0.5', '0.75', '1')) + 
  facet_grid(.~ Method) +
  theme_minimal() + 
  theme(legend.position = 'none',
        axis.text = element_text(size = 11)) +
  labs(x = 'Pairwise IBD by WGS', y = 'Pairwise IBD by amplicon panels', title = 'A) Pairwise IBD measured using WGS vs amplicon Panels')



# RMSE_IBD_plot = IBD_values_LAC %>%
#   pivot_longer(cols = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'),
#                names_to = 'Method',
#                values_to = 'IBD') %>%
#   mutate(Method = factor(Method, levels = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq')),
#          Error = (IBD - IBDwgs)^2,
#          IBD_cats = case_when(
#            IBDwgs < .05 ~ '0.00 - 0.049',
#            IBDwgs >= .05 & IBDwgs < .1 ~ '0.05 - 0.099',
#            IBDwgs >= .1 & IBDwgs < .15 ~ '0.10 - 0.149',
#            IBDwgs >= .15 & IBDwgs < .2 ~ '0.15 - 0.199',
#            IBDwgs >= .2 & IBDwgs < .25 ~ '0.20 - 0.249',
#            IBDwgs >= .25 & IBDwgs < .3 ~ '0.25 - 0.299',
#            IBDwgs >= .3 & IBDwgs < .35 ~ '0.30 - 0.349',
#            IBDwgs >= .35 & IBDwgs < .4 ~ '0.35 - 0.399',
#            IBDwgs >= .4 & IBDwgs < .45 ~ '0.40 - 0.449',
#            IBDwgs >= .45 & IBDwgs < .5 ~ '0.45 - 0.499',
#            IBDwgs >= .5 & IBDwgs < .55 ~ '0.50 - 0.549',
#            IBDwgs >= .55 & IBDwgs < .6 ~ '0.55 - 0.599',
#            IBDwgs >= .6 & IBDwgs < .65 ~ '0.60 - 0.649',
#            IBDwgs >= .65 & IBDwgs < .7 ~ '0.65 - 0.699',
#            IBDwgs >= .7 & IBDwgs < .75 ~ '0.70 - 0.749',
#            IBDwgs >= .75 & IBDwgs < .8 ~ '0.75 - 0.799',
#            IBDwgs >= .8 & IBDwgs < .85 ~ '0.80 - 0.849',
#            IBDwgs >= .85 & IBDwgs < .9 ~ '0.85 - 0.899',
#            IBDwgs >= .9 & IBDwgs < .95 ~ '0.90 - 0.949',
#            IBDwgs >= .95 ~ '0.95 - 1.000'
#          )
#          ) %>%
#   summarise(RMSE = sqrt(mean(Error)), 
#             lower_limit = sqrt(n()*mean(Error)/qchisq(p = .05, df = n() - 1, lower.tail = TRUE)), 
#             upper_limit = sqrt(n()*mean(Error)/qchisq(p = .05, df = n() - 1, lower.tail = FALSE)), 
#             .by = c(Method, IBD_cats)) %>%
#   #summarise(RMSE = sqrt(mean(Error)), se = sqrt(sd(Error)), .by = c(Method, IBD_cats)) %>%
#   #mutate(lower_limit = ifelse(RMSE - se >= 0, RMSE - se, 0), upper_limit  = RMSE + se) %>%
#   ggplot(aes(y = RMSE, x = IBD_cats, group = Method)) +
#   geom_point()+
#   geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit))+
#   geom_line()+
#   scale_y_continuous(limits = c(0, .7), breaks = seq(0, 1, .1)) + 
#   facet_grid(.~ Method)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 90))+
#   labs(x = 'Pairwise IBD by WGS',
#        title = 'B) Root mean squared error (RMSE) of pairwised IBD with respect to WGS')




RMSE_IBD_plot = IBD_values_LAC %>%
  pivot_longer(cols = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'),
               names_to = 'Method',
               values_to = 'IBD') %>%
  mutate(Method = factor(Method, levels = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq')),
         Error = (IBD - IBDwgs)^2,
         IBD_cats = case_when(
           IBDwgs < .05 ~ '0.00 - 0.049',
           IBDwgs >= .05 & IBDwgs < .1 ~ '0.05 - 0.099',
           IBDwgs >= .1 & IBDwgs < .15 ~ '0.10 - 0.149',
           IBDwgs >= .15 & IBDwgs < .2 ~ '0.15 - 0.199',
           IBDwgs >= .2 & IBDwgs < .25 ~ '0.20 - 0.249',
           IBDwgs >= .25 & IBDwgs < .3 ~ '0.25 - 0.299',
           IBDwgs >= .3 & IBDwgs < .35 ~ '0.30 - 0.349',
           IBDwgs >= .35 & IBDwgs < .4 ~ '0.35 - 0.399',
           IBDwgs >= .4 & IBDwgs < .45 ~ '0.40 - 0.449',
           IBDwgs >= .45 & IBDwgs < .5 ~ '0.45 - 0.499',
           IBDwgs >= .5 & IBDwgs < .55 ~ '0.50 - 0.549',
           IBDwgs >= .55 & IBDwgs < .6 ~ '0.55 - 0.599',
           IBDwgs >= .6 & IBDwgs < .65 ~ '0.60 - 0.649',
           IBDwgs >= .65 & IBDwgs < .7 ~ '0.65 - 0.699',
           IBDwgs >= .7 & IBDwgs < .75 ~ '0.70 - 0.749',
           IBDwgs >= .75 & IBDwgs < .8 ~ '0.75 - 0.799',
           IBDwgs >= .8 & IBDwgs < .85 ~ '0.80 - 0.849',
           IBDwgs >= .85 & IBDwgs < .9 ~ '0.85 - 0.899',
           IBDwgs >= .9 & IBDwgs < .95 ~ '0.90 - 0.949',
           IBDwgs >= .95 ~ '0.95 - 1.000'
         )
  ) %>%
  summarise(RMSE = sqrt(mean(Error)), 
            lower_limit = sqrt(n()*mean(Error)/qchisq(p = .05, df = n() - 1, lower.tail = TRUE)), 
            upper_limit = sqrt(n()*mean(Error)/qchisq(p = .05, df = n() - 1, lower.tail = FALSE)), 
            .by = c(Method, IBD_cats)) %>%
  #summarise(RMSE = sqrt(mean(Error)), se = sqrt(sd(Error)), .by = c(Method, IBD_cats)) %>%
  #mutate(lower_limit = ifelse(RMSE - se >= 0, RMSE - se, 0), upper_limit  = RMSE + se) %>%
  ggplot(aes(y = RMSE, x = IBD_cats, color = Method, group = Method)) +
  geom_point()+
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = .1)+
  geom_line()+
  geom_hline(yintercept = .1, linetype = 2) +
  scale_y_continuous(limits = c(0, .5), breaks = seq(0, 1, .1)) + 
  scale_color_manual(values = c('gray10', 'firebrick3', 'dodgerblue3')) + 
  #facet_grid(.~ Method)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = 'inside', legend.position.inside = c(.85, .8), legend.title = element_blank(),
        axis.text = element_text(size = 11)
        )+
  labs(x = 'Pairwise IBD by WGS',
       title = 'B) Root mean squared error (RMSE) of pairwised IBD with respect to WGS')





CV_IBD_plot = IBD_values_LAC %>%
  pivot_longer(cols = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'),
               names_to = 'Method',
               values_to = 'IBD') %>%
  mutate(Method = factor(Method, levels = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq')),
         #Error = sqrt((IBD - IBDwgs)^2),
         Error = IBD,
         IBD_cats = case_when(
           IBDwgs < .05 ~ '0.00 - 0.049',
           IBDwgs >= .05 & IBDwgs < .1 ~ '0.05 - 0.099',
           IBDwgs >= .1 & IBDwgs < .15 ~ '0.10 - 0.149',
           IBDwgs >= .15 & IBDwgs < .2 ~ '0.15 - 0.199',
           IBDwgs >= .2 & IBDwgs < .25 ~ '0.20 - 0.249',
           IBDwgs >= .25 & IBDwgs < .3 ~ '0.25 - 0.299',
           IBDwgs >= .3 & IBDwgs < .35 ~ '0.30 - 0.349',
           IBDwgs >= .35 & IBDwgs < .4 ~ '0.35 - 0.399',
           IBDwgs >= .4 & IBDwgs < .45 ~ '0.40 - 0.449',
           IBDwgs >= .45 & IBDwgs < .5 ~ '0.45 - 0.499',
           IBDwgs >= .5 & IBDwgs < .55 ~ '0.50 - 0.549',
           IBDwgs >= .55 & IBDwgs < .6 ~ '0.55 - 0.599',
           IBDwgs >= .6 & IBDwgs < .65 ~ '0.60 - 0.649',
           IBDwgs >= .65 & IBDwgs < .7 ~ '0.65 - 0.699',
           IBDwgs >= .7 & IBDwgs < .75 ~ '0.70 - 0.749',
           IBDwgs >= .75 & IBDwgs < .8 ~ '0.75 - 0.799',
           IBDwgs >= .8 & IBDwgs < .85 ~ '0.80 - 0.849',
           IBDwgs >= .85 & IBDwgs < .9 ~ '0.85 - 0.899',
           IBDwgs >= .9 & IBDwgs < .95 ~ '0.90 - 0.949',
           IBDwgs >= .95 ~ '0.95 - 1.000'
         )
  ) %>%
  summarise(CV = sd(Error, na.rm = T)/mean(Error, na.rm = T), 
            lower_limit = sd(Error, na.rm = T)/mean(Error, na.rm = T) + qt(p = .05, df = n() - 1, lower.tail = TRUE) * (sd(Error, na.rm = T)/mean(Error, na.rm = T)/sqrt(2*n())), 
            upper_limit = sd(Error, na.rm = T)/mean(Error, na.rm = T) + qt(p = .05, df = n() - 1, lower.tail = FALSE) * (sd(Error, na.rm = T)/mean(Error, na.rm = T)/sqrt(2*n())), 
            .by = c(Method, IBD_cats)) %>%
  #summarise(RMSE = sqrt(mean(Error)), se = sqrt(sd(Error)), .by = c(Method, IBD_cats)) %>%
  #mutate(lower_limit = ifelse(RMSE - se >= 0, RMSE - se, 0), upper_limit  = RMSE + se) %>%
  ggplot(aes(y = CV, x = IBD_cats, color = Method, group = Method)) +
  geom_point()+
  geom_errorbar(aes(ymin = lower_limit, ymax = upper_limit), width = .1)+
  geom_line() +
  geom_hline(yintercept = .1, linetype = 2) +
  scale_y_continuous(limits = c(0, .2), breaks = seq(0, 1, .05)) + 
  scale_color_manual(values = c('gray10', 'firebrick3', 'dodgerblue3')) + 
  #facet_grid(.~ Method)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = 'inside', legend.position.inside = c(.85, .8), legend.title = element_blank(),
        axis.text = element_text(size = 11)
  ) +
  labs(x = 'Pairwise IBD by WGS',
       title = 'C) Coefficient of variation (CV) of pairwised IBD with respect to WGS')

qt(p = .05, df = 3000 - 1, lower.tail = TRUE)
qt(p = .05, df = 3000 - 1, lower.tail = FALSE)


ggdraw() +
  draw_plot(IBD_corr_plot,
            x = 0, 
            width = 1,
            y = .66,
            height = .34) +
  draw_plot(RMSE_IBD_plot,
            x = 0, 
            width = 1,
            y = .33,
            height = .33) +
  draw_plot(CV_IBD_plot,
            x = 0, 
            width = 1,
            y = 0,
            height = .33)


monoclonals_PvGTSeq_LAC_loci@metadata%>%
  filter(Sample_id %in% selected_samples_ids, site_of_collection_snl0 == 'Brazil') %>% View()




Error_test_data = IBD_values_LAC %>%
  pivot_longer(cols = c('PvGTSeq', 'PvAmpliSeq', 'rhAmpSeq'),
               names_to = 'Method',
               values_to = 'IBD') %>%
  mutate(Error = (IBD - IBDwgs)^2,
         IBD_cats = case_when(
           IBDwgs < .05 ~ '0.00 - 0.049',
           IBDwgs >= .05 & IBDwgs < .1 ~ '0.05 - 0.099',
           IBDwgs >= .1 & IBDwgs < .15 ~ '0.10 - 0.149',
           IBDwgs >= .15 & IBDwgs < .2 ~ '0.15 - 0.199',
           IBDwgs >= .2 & IBDwgs < .25 ~ '0.20 - 0.249',
           IBDwgs >= .25 & IBDwgs < .3 ~ '0.25 - 0.299',
           IBDwgs >= .3 & IBDwgs < .35 ~ '0.30 - 0.349',
           IBDwgs >= .35 & IBDwgs < .4 ~ '0.35 - 0.399',
           IBDwgs >= .4 & IBDwgs < .45 ~ '0.40 - 0.449',
           IBDwgs >= .45 & IBDwgs < .5 ~ '0.45 - 0.499',
           IBDwgs >= .5 & IBDwgs < .55 ~ '0.50 - 0.549',
           IBDwgs >= .55 & IBDwgs < .6 ~ '0.55 - 0.599',
           IBDwgs >= .6 & IBDwgs < .65 ~ '0.60 - 0.649',
           IBDwgs >= .65 & IBDwgs < .7 ~ '0.65 - 0.699',
           IBDwgs >= .7 & IBDwgs < .75 ~ '0.70 - 0.749',
           IBDwgs >= .75 & IBDwgs < .8 ~ '0.75 - 0.799',
           IBDwgs >= .8 & IBDwgs < .85 ~ '0.80 - 0.849',
           IBDwgs >= .85 & IBDwgs < .9 ~ '0.85 - 0.899',
           IBDwgs >= .9 & IBDwgs < .95 ~ '0.90 - 0.949',
           IBDwgs >= .95 ~ '0.95 - 1.000'
         )
  )


Error_test_data %>%
  summarise(Mean = sqrt(mean(Error)), .by = c(Method, IBD_cats)) %>% View()


Error_test_result = NULL


categories = c('0.00 - 0.049',
               '0.05 - 0.099',
               '0.10 - 0.149',
               '0.15 - 0.199',
               '0.20 - 0.249',
               '0.25 - 0.299',
               '0.30 - 0.349',
               '0.35 - 0.399',
               '0.40 - 0.449',
               '0.45 - 0.499',
               '0.50 - 0.549',
               '0.55 - 0.599',
               '0.60 - 0.649',
               '0.65 - 0.699',
               '0.70 - 0.749',
               '0.75 - 0.799',
               '0.80 - 0.849',
               '0.85 - 0.899',
               '0.90 - 0.949',
               '0.95 - 1.000')



for(category in categories){
  
  for(method in c('PvAmpliSeq', 'rhAmpSeq')){
    
    PvGTSeq_error = Error_test_data %>% filter(Method == 'PvGTSeq', IBD_cats == category) %>% select(Error) %>% unlist()

    Method2_error = Error_test_data %>% filter(Method == method, IBD_cats == category) %>% select(Error) %>% unlist()
    
    t_test_result = t.test(x = PvGTSeq_error, y = Method2_error, alternative = 'less', paired = T)
    
    Error_test_result = rbind(Error_test_result,
                              data.frame(comparison = paste0('PvGTSeq - ', method),
                                         IBD_category = category,
                                         difference = t_test_result$estimate,
                                         p_value = t_test_result$p.value))
  }
  
}












IBD_evectors_wgs_LAC = GRM_evectors(dist_table = Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD, 
                                    k = nrow(monoclonals_LAC_01SNPsbiallelic_75_loci@metadata %>%
                                               filter(Sample_id %in% selected_samples_ids,
                                                      !(Sample_id %in% c('ERR2309697', 'ERR2678957', 'ERR2678961')),
                                                      !(Sample_id %in% c('ERR2351947', 'ERR2351948', 'ERR2678958', 'ERR2678959', 'ERR2678960')),
                                                      !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua')))
                                    ),
                                    metadata = monoclonals_LAC_01SNPsbiallelic_75_loci@metadata%>%
                                      filter(Sample_id %in% selected_samples_ids,
                                             !(Sample_id %in% c('ERR2309697', 'ERR2678957', 'ERR2678961')),
                                             !(Sample_id %in% c('ERR2351947', 'ERR2351948', 'ERR2678958', 'ERR2678959', 'ERR2678960')),
                                             !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua'))),
                                    Pop = 'site_of_collection_snl0',
                                    fastGRM = F,
                                    cor = F,
                                    method = 'princomp') # fastSVD 



Paired12 = brewer.pal(11, 'Paired')

PCoA_no_clonal_wgs_LAC = IBD_evectors_wgs_LAC$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = -PC2, 
             color = factor(site_of_collection_snl0, levels = c('Brazil', 'Venezuela', 'Guyana', 'FrenchGuiana', 'Trinidad',
                                                                'Peru', 'Colombia', 'Panama', 'Honduras', 'Guatemala', 'Mexico', 'USA'))
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_manual(values = c(Paired12, 'gray10'))+
  theme_minimal()+
  labs(title = paste0('A) LAC countries using\n',
                      ncol(monoclonals_LAC_01SNPsbiallelic_75_loci@loci_table),
                      ' SNPs from WGS'
  ),
  x = paste0('PCoA 1 (', round(IBD_evectors_wgs_LAC$contrib[1], 2), '%)'),
  y = paste0('PCoA 2 (', round(IBD_evectors_wgs_LAC$contrib[2], 2), '%)'),
  color = 'Study sites',
  shape = 'Study sites'
  ) +
  theme(#legend.position = 'inside', legend.position.inside = c(.7, .5),
    legend.position = 'bottom',    
    #legend.title = element_blank(),
    legend.justification = 'left',
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    title = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 4,
                              theme = theme(legend.title.position = 'left'
                              )))

PCoA_no_clonal_wgs_LAC

Monoclonal_LAC_pairs_01SNPsbiallelic_75_IBD %>% filter(Yi_site_of_collection_snl0 == 'USA' | Yj_site_of_collection_snl0 == 'USA') %>% filter(rhat > .5) %>% arrange(rhat)
