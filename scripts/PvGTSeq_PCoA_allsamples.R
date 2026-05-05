
source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_required_libraries.R')
source('~/Documents/Github/Plasmodium_WGS_analysis/functions_libraries/rGenome_functions.R')
source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_functions.R')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/Rcpp_functions.cpp')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/hmmloglikelihood.cpp')

# PvGTSeq ----

PvGTSeq_markers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_markersTable.csv')

PvGTSeq_markers %<>% filter(use !='DRS')

# Extract PvGTSeq coordinates from rGenome object----

PvGTSeq_coordinates = NULL
for(amplicon in 1:nrow(PvGTSeq_markers)){
  amplicon_chromosome = gsub('v1$', 'v2', PvGTSeq_markers[amplicon,][['chromosome']])
  amplicon_start = PvGTSeq_markers[amplicon,][['start']]
  amplicon_end = PvGTSeq_markers[amplicon,][['end']]
  
  PvGTSeq_coordinates = c(PvGTSeq_coordinates,
                             paste(amplicon_chromosome,
                                   amplicon_start:amplicon_end, sep = '_'))
  
}

# load or generate rGenome object----

PvGTSeq_rGenome = read_rGenome('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_all_samples_PvGTSeq_rGenome', 
                               format = 'tsv',
                               sep = '\t')


sum(colnames(PvGTSeq_rGenome@gt) != PvGTSeq_rGenome@metadata$Sample_id)
colnames(PvGTSeq_rGenome@gt) = PvGTSeq_rGenome@metadata$Sample_id

# Verify that there are no duplicated samples per individual----

sum(duplicated(PvGTSeq_rGenome@metadata$metadata_id))

# Filter out samples----

PvGTSeq_rGenome@metadata %>%
  summarise(nSamples = n(), 
            .by = c(batch, site_of_collection_world_region, site_of_collection_snl0, country)) %>%
  arrange(batch, site_of_collection_world_region, site_of_collection_snl0)
  


PvGTSeq_rGenome = 
  filter_samples(PvGTSeq_rGenome,
                 (PvGTSeq_rGenome@metadata$PvGTSeq_SampAmpRate >= .75 &
                    PvGTSeq_rGenome@metadata$PvAmpliSeq_SampAmpRate >= .75 &
                    PvGTSeq_rGenome@metadata$rhAmpSeq_SampAmpRate >= .75 &
                    PvGTSeq_rGenome@metadata$WGS_SampAmpRate >= .75))


PvGTSeq_rGenome = filter_samples(PvGTSeq_rGenome,
                                         !(
                                           PvGTSeq_rGenome@metadata$batch %in% c('Duraisingh', 'Ancient') |
                                             PvGTSeq_rGenome@metadata$site_of_collection_snl0 %in% c('P. simium', 'Unknown')
                                         ))

PvGTSeq_rGenome@metadata %>%
  summarise(nSamples = n(), 
            .by = c(site_of_collection_world_region, site_of_collection_snl0)) %>%
  arrange(site_of_collection_world_region, site_of_collection_snl0)

dim(PvGTSeq_rGenome@gt)



# Keep only positions for geographic differentiation---

dim(PvGTSeq_rGenome@gt)

PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, v = rownames(PvGTSeq_rGenome@gt) %in% PvGTSeq_coordinates)

dim(PvGTSeq_rGenome@gt)


# Check number of amplicons covered----

PvGTSeq_rGenome@loci_table$amplicon = sapply(1:nrow(PvGTSeq_rGenome@loci_table), function(i){
  
  variant_site_chromosome = gsub('v2$', 'v1', PvGTSeq_rGenome@loci_table[i, ][['CHROM']])
  variant_site_position = PvGTSeq_rGenome@loci_table[i, ][['POS']]
  
  PvGTSeq_markers[PvGTSeq_markers$chromosome == variant_site_chromosome &
                    PvGTSeq_markers$start <= variant_site_position &
                    PvGTSeq_markers$end >= variant_site_position,
  ][['amplicon']]
  
})

print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome@loci_table$amplicon))))


# Update allele frequencies and remove monomorphic sites----

PvGTSeq_rGenome = update_allele_lables(PvGTSeq_rGenome, n = 10)

## Remove monomorphic sites----
dim(PvGTSeq_rGenome@gt)

PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, 
                              v = PvGTSeq_rGenome@loci_table$Cardinality > 1)

dim(PvGTSeq_rGenome@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome@loci_table$amplicon))))

# Calculate amplification rate per Loci-----
PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate = LocusAmplRate(PvGTSeq_rGenome, update = F)

min(PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate)

PvGTSeq_rGenome@loci_table %>%
  ggplot(aes(x = Locus_Ampl_Rate))+
  geom_histogram()


# # ## Remove loci with less than 50% of genome coverage----
# PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome,
#                                  v = PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate >= .75)
# 
# dim(PvGTSeq_rGenome@gt)
# print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome@loci_table$amplicon))))
# names(PvGTSeq_rGenome@metadata)

# Identify the type of polymophism of the variant site----
PvGTSeq_rGenome@loci_table$type_of_polymorphism = get_type_of_polymorphism(PvGTSeq_rGenome)

unique(PvGTSeq_rGenome@loci_table$type_of_polymorphism)

## Remove Homopolymers and STRs----
PvGTSeq_rGenome = 
  filter_loci(PvGTSeq_rGenome,
              v = !(PvGTSeq_rGenome@loci_table$type_of_polymorphism %in% 
                      c('INDEL:Homopolymer', 'INDEL:Dinucleotide_STR')))

dim(PvGTSeq_rGenome@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome@loci_table$amplicon))))

# Calculate minor allele frequency----

PvGTSeq_major_freq = sapply(PvGTSeq_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})

max(PvGTSeq_major_freq)


dim(PvGTSeq_rGenome@gt)

PvGTSeq_rGenome@loci_table$major_freq = PvGTSeq_major_freq

sum(PvGTSeq_major_freq <= 1 - 5/ncol(PvGTSeq_rGenome@gt))

PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, 
                                 v = PvGTSeq_rGenome@loci_table$major_freq <= 1 - 5/ncol(PvGTSeq_rGenome@gt))


dim(PvGTSeq_rGenome@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome@loci_table$amplicon))))

# Keep monoclonal samples with 75% of wgs coverage----

monoclonals_Pv_all_samples = unlist(read.table('~/Documents/Github/PvGTSeq_paper_draft/Outputs/monoclonals_Pv_all_samples.bed',
                                               header = F,
                                               sep = '\t'
                                               ))

names(monoclonals_Pv_all_samples) = NULL

sum(!(monoclonals_Pv_all_samples %in% PvGTSeq_rGenome@metadata$Sample_id))

PvGTSeq_rGenome@metadata %<>% mutate(
  Clonality = case_when(
    Sample_id %in% monoclonals_Pv_all_samples ~ 'Monoclonal',
    .default = 'Polyclonal'
  )
)


PvGTSeq_rGenome_biallelic_SNPs = filter_loci(PvGTSeq_rGenome,
                                             (PvGTSeq_rGenome@loci_table$type_of_polymorphism == 'SNP' &
                                                PvGTSeq_rGenome@loci_table$Cardinality == 2
                                              ))

dim(PvGTSeq_rGenome_biallelic_SNPs@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(PvGTSeq_rGenome_biallelic_SNPs@loci_table$amplicon))))



PvGTSeq_rGenome@metadata$Fws = get_Fws_rGenome(PvGTSeq_rGenome_biallelic_SNPs)
PvGTSeq_rGenome@metadata$ObsHet = get_ObsHet(PvGTSeq_rGenome_biallelic_SNPs, by = 'sample')

PvGTSeq_rGenome@metadata %>%
  ggplot(aes(x = ObsHet, y = Fws, color = Clonality)) + 
  geom_point() +
  geom_vline(xintercept = 0.0125)+
  geom_hline(yintercept = .95)+
  scale_y_continuous(limits = c(0, 1))


# PvGTSeq_rGenome@metadata %>%
#   summarise(nMonoclonals = sum(Fws >= .95 & ObsHet <= 0.0125),
#             nPolyclonals = sum(Fws < .95 | ObsHet > 0.0125),
#             nSamples = n(),
#             .by = site_of_collection_snl0)
# 
# 
# 
# PvGTSeq_rGenome@metadata %<>%
#   mutate(Clonality = case_when(
#     Fws >= .95 & ObsHet <= 0.0125 ~ 'Monoclonal',
#     Fws < .95 | ObsHet > 0.0125 ~ 'Polyclonal'))


monoclonals_PvGTSeq_rGenome = 
  filter_samples(PvGTSeq_rGenome,
                 PvGTSeq_rGenome@metadata$Clonality == 'Monoclonal')

dim(monoclonals_PvGTSeq_rGenome@gt)


# Update allele frequencies and remove monomorphic sites----

monoclonals_PvGTSeq_rGenome = update_allele_lables(monoclonals_PvGTSeq_rGenome, n = 1)

## Remove monomorphic sites----

sum(monoclonals_PvGTSeq_rGenome@loci_table$Cardinality == 1)

monoclonals_PvGTSeq_rGenome = filter_loci(monoclonals_PvGTSeq_rGenome, 
                                             v = monoclonals_PvGTSeq_rGenome@loci_table$Cardinality > 1)

dim(monoclonals_PvGTSeq_rGenome@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(monoclonals_PvGTSeq_rGenome@loci_table$amplicon))))

# Calculate minor allele frequency
PvGTSeq_major_freq = sapply(monoclonals_PvGTSeq_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvGTSeq_major_freq)

monoclonals_PvGTSeq_rGenome@loci_table$major_freq = PvGTSeq_major_freq

sum(monoclonals_PvGTSeq_rGenome@loci_table$major_freq <= 1 - 5/ncol(monoclonals_PvGTSeq_rGenome@gt))

monoclonals_PvGTSeq_rGenome = filter_loci(monoclonals_PvGTSeq_rGenome, 
                                             v = monoclonals_PvGTSeq_rGenome@loci_table$major_freq <= 1 - 5/ncol(monoclonals_PvGTSeq_rGenome@gt))

dim(monoclonals_PvGTSeq_rGenome@gt)
print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(monoclonals_PvGTSeq_rGenome@loci_table$amplicon))))

monoclonals_PvGTSeq_rGenome@metadata$sample_ampl_rate = SampleAmplRate(monoclonals_PvGTSeq_rGenome, update = F)

monoclonals_PvGTSeq_rGenome@metadata %>% 
  ggplot(aes(x = sample_ampl_rate)) +
  geom_histogram()


# Create loci object----

monoclonals_PvGTSeq_gt = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_PvGTSeq_rGenome@gt))

monoclonals_PvGTSeq_gt = t(monoclonals_PvGTSeq_gt)

monoclonals_PvGTSeq_loci = create_loci(loci_table = monoclonals_PvGTSeq_gt,
                                          metadata = monoclonals_PvGTSeq_rGenome@metadata)

dim(monoclonals_PvGTSeq_loci@loci_table)


PvGTSeq_allele_freqs = get_allele_freq(monoclonals_PvGTSeq_loci)

monoclonals_PvGTSeq_loci@freq_table = PvGTSeq_allele_freqs

nrow(PvGTSeq_allele_freqs)
dim(monoclonals_PvGTSeq_loci@loci_table)


monoclonals_PvGTSeq_loci@markers = 
  monoclonals_PvGTSeq_rGenome@loci_table %>% S4Vectors::rename('CHROM' = 'chromosome',
                                                       'POS' = 'pos')


rownames(monoclonals_PvGTSeq_loci@loci_table) = monoclonals_PvGTSeq_loci@metadata$Sample_id

# Monoclonal pairs ----

PvGTSeq_Monoclonals = monoclonals_PvGTSeq_loci@metadata$Sample_id

PvGTSeq_Monoclonal_pairs = as.data.frame(t(combn(PvGTSeq_Monoclonals, 2)))

names(PvGTSeq_Monoclonal_pairs) = c('Yi', 'Yj')

PvGTSeq_Monoclonal_pairs = left_join(
  PvGTSeq_Monoclonal_pairs,
  monoclonals_PvGTSeq_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(PvGTSeq_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


PvGTSeq_Monoclonal_pairs = left_join(
  PvGTSeq_Monoclonal_pairs,
  monoclonals_PvGTSeq_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(PvGTSeq_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')

nrow(PvGTSeq_Monoclonal_pairs)/5000

# Calculate IBD for world regions---

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:10000){
    start_time = Sys.time()
    Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvGTSeq_loci, parallel = T, pairs = PvGTSeq_Monoclonal_pairs, max_k = 20,
              w = w, n = 10000
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %<>%
    mutate(euDist = 1 - rhat)
  
  Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
    left_join(Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
              PvGTSeq_Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  write.csv(Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')
  
}



names(Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD)


Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')



evectors_PvGTSeq_world_regions = GRM_evectors(dist_table = Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD, 
                                                 k = nrow(monoclonals_PvGTSeq_loci@metadata %>%
                                                            filter(Sample_id %in% selected_samples_ids,
                                                                   !(Sample_id %in% AF_outliers),
                                                                   !(Sample_id %in% Asia_outliers),
                                                                   !(Sample_id %in% LAC_outliers),
                                                                   possible_origin != 'Imported')
                                                 ),
                                                 metadata = monoclonals_PvGTSeq_loci@metadata%>%
                                                   filter(Sample_id %in% selected_samples_ids,
                                                          !(Sample_id %in% AF_outliers),
                                                          !(Sample_id %in% Asia_outliers),
                                                          !(Sample_id %in% LAC_outliers),
                                                          possible_origin != 'Imported'),
                                                 Pop = 'site_of_collection_world_region',
                                                 method = 'princomp'
)




PCoA_no_clonal_PvGTSeq_world_regions = evectors_PvGTSeq_world_regions$eigenvector %>%
  ggplot(aes(x = -PC1, 
             y = PC2, 
             color = site_of_collection_world_region
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('B) World regions using ', 
                      nrow(monoclonals_PvGTSeq_loci@markers), 
                      ' variant sites\n across ',
                      length(unique(monoclonals_PvGTSeq_loci@markers$amplicon)),
                      '/', nrow(PvGTSeq_markers),
                      ' PvGTSeq amplicons'),
    x = paste0('PCoA 1 (', round(evectors_PvGTSeq_world_regions$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvGTSeq_world_regions$contrib[2], 2), '%)'),
    color = 'Study sites',
    shape = 'Study sites'
  ) +
  theme(legend.position = 'none',
        #legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12))# +
  # guides(color = guide_legend(ncol = 4,
  #                             theme = theme(legend.title.position = 'left'
  #                             )))

PCoA_no_clonal_PvGTSeq_world_regions




# Calculate IBD for LAC---

## Create loci object----

dim(monoclonals_PvGTSeq_rGenome@gt) 

monoclonals_PvGTSeq_LAC_rGenome = 
  filter_samples(monoclonals_PvGTSeq_rGenome,
                 monoclonals_PvGTSeq_rGenome@metadata$site_of_collection_world_region == 'LAC'
                 )

monoclonals_PvGTSeq_LAC_rGenome@loci_table$Cardinality = NULL
monoclonals_PvGTSeq_LAC_rGenome@loci_table$Alleles = NULL
monoclonals_PvGTSeq_LAC_rGenome@loci_table$Allele_Counts = NULL
monoclonals_PvGTSeq_LAC_rGenome@loci_table$major_freq = NULL

monoclonals_PvGTSeq_LAC_rGenome = update_allele_lables(monoclonals_PvGTSeq_LAC_rGenome)

sum(monoclonals_PvGTSeq_LAC_rGenome@loci_table$Cardinality >= 2)


monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
            monoclonals_PvGTSeq_LAC_rGenome@loci_table$Cardinality >= 2)

dim(monoclonals_PvGTSeq_LAC_rGenome@gt)

print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(monoclonals_PvGTSeq_LAC_rGenome@loci_table$amplicon))))

monoclonals_PvGTSeq_LAC_rGenome@metadata  %>%
  summarise(nSamples = n(), 
            .by = c(site_of_collection_world_region, site_of_collection_snl0)) %>%
  arrange(site_of_collection_world_region, site_of_collection_snl0)



# Calculate minor allele frequency
PvGTSeq_major_freq = sapply(monoclonals_PvGTSeq_LAC_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvGTSeq_major_freq)
1 - 5/ncol(monoclonals_PvGTSeq_LAC_rGenome@gt)

sum(PvGTSeq_major_freq <= 1 - 5/ncol(monoclonals_PvGTSeq_LAC_rGenome@gt))

monoclonals_PvGTSeq_LAC_rGenome@loci_table$major_freq = PvGTSeq_major_freq

monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
                                             v = monoclonals_PvGTSeq_LAC_rGenome@loci_table$major_freq <= 1 - 5/ncol(monoclonals_PvGTSeq_LAC_rGenome@gt))

dim(monoclonals_PvGTSeq_LAC_rGenome@gt)

print(paste0('the number of amplicons with partial coverage in the data set is: ', length(unique(monoclonals_PvGTSeq_LAC_rGenome@loci_table$amplicon))))

# sum(monoclonals_PvGTSeq_LAC_rGenome@loci_table$type_of_polymorphism != 'SNP')
# 
# monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
#                                               monoclonals_PvGTSeq_LAC_rGenome@loci_table$type_of_polymorphism == 'SNP')

dim(monoclonals_PvGTSeq_LAC_rGenome@gt)

#monoclonals_PvGTSeq_gt = gsub(':\\d+', '', monoclonals_PvGTSeq_LAC_rGenome@gt)
monoclonals_PvGTSeq_gt = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_PvGTSeq_LAC_rGenome@gt))

monoclonals_PvGTSeq_gt = t(monoclonals_PvGTSeq_gt)

monoclonals_PvGTSeq_LAC_loci = create_loci(loci_table = monoclonals_PvGTSeq_gt,
                                       metadata = monoclonals_PvGTSeq_LAC_rGenome@metadata)

dim(monoclonals_PvGTSeq_LAC_loci@loci_table)


# monoclonals_PvGTSeq_LAC_loci@metadata %<>%
#   mutate(Strata = case_when(
#     site_of_collection_snl0 %in% c('Colombia') ~ 'Colombia',
#     site_of_collection_snl0 %in% c('El Salvador', 'Honduras', 'Mexico', 'Nicaragua', 'Panama') ~ 'Meso America',
#     site_of_collection_snl0 %in% c('Brazil', 'Peru') ~ 'Amazon',
#     site_of_collection_snl0 %in% c('Guyana', 'Venezuela') ~ 'Guainia Shield'
#   ))
# 
# monoclonals_PvGTSeq_LAC_loci_no_clonal = filter_samples(monoclonals_PvGTSeq_LAC_loci,
#                                                         monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id %in% selected_LAC_samples_ids
#                                                         )

PvGTSeq_allele_freqs = get_allele_freq(monoclonals_PvGTSeq_LAC_loci#, by = 'Strata'
                                       )

monoclonals_PvGTSeq_LAC_loci@freq_table = PvGTSeq_allele_freqs

dim(monoclonals_PvGTSeq_LAC_loci@loci_table)


monoclonals_PvGTSeq_LAC_loci@markers = 
  monoclonals_PvGTSeq_LAC_rGenome@loci_table %>% S4Vectors::rename('CHROM' = 'chromosome',
                                                    'POS' = 'pos')


rownames(monoclonals_PvGTSeq_LAC_loci@loci_table) = monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id


# Monoclonal pairs ----

PvGTSeq_Monoclonals_LAC = monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id

PvGTSeq_LAC_Monoclonal_pairs = as.data.frame(t(combn(PvGTSeq_Monoclonals_LAC, 2)))

names(PvGTSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj')

PvGTSeq_LAC_Monoclonal_pairs = left_join(
  PvGTSeq_LAC_Monoclonal_pairs,
  monoclonals_PvGTSeq_LAC_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(PvGTSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


PvGTSeq_LAC_Monoclonal_pairs = left_join(
  PvGTSeq_LAC_Monoclonal_pairs,
  monoclonals_PvGTSeq_LAC_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(PvGTSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')


nrow(PvGTSeq_LAC_Monoclonal_pairs)/5000

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:5000){
    start_time = Sys.time()
    Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvGTSeq_LAC_loci, parallel = T, pairs = PvGTSeq_LAC_Monoclonal_pairs, max_k = 20,
              #freq_table = monoclonals_PvGTSeq_LAC_loci@freq_table, by = 'Strata',
              w = w, n = 5000
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %<>%
    mutate(euDist = 1 - rhat)
  
  Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
    left_join(Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
              PvGTSeq_LAC_Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #sum(Monoclonal_LAC_pairs_01SNPsbiallelic_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')
  
}



names(Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD)


Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')


dim(monoclonals_PvGTSeq_LAC_loci@loci_table)


evectors_PvGTSeq_LAC = GRM_evectors(dist_table = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD, 
                                       k = nrow(monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                                  filter(Sample_id %in% selected_samples_ids,
                                                         !(Sample_id %in% LAC_outliers),
                                                         #!(Sample_id %in% Brazil_outliers),
                                                         !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana')))
                                       ),
                                       metadata = monoclonals_PvGTSeq_LAC_loci@metadata%>%
                                         filter(Sample_id %in% selected_samples_ids,
                                                !(Sample_id %in% LAC_outliers),
                                                #!(Sample_id %in% Brazil_outliers),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                       Pop = 'site_of_collection_snl0',
                                       method = 'princomp')



PCoA_no_clonal_PvGTSeq_LAC = evectors_PvGTSeq_LAC$eigenvector %>%
  ggplot(aes(x = PC1, 
             y = -PC2, 
             color = site_of_collection_snl0
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('B) LAC countries using ', 
                      nrow(monoclonals_PvGTSeq_LAC_loci@markers), 
                      ' variant sites\nacross ',
                      length(unique(monoclonals_PvGTSeq_LAC_loci@markers$amplicon)),
                      '/', nrow(PvGTSeq_markers),
                      ' PvGTSeq amplicons'),
    x = paste0('PCoA 1 (', round(evectors_PvGTSeq_LAC$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvGTSeq_LAC$contrib[2], 2), '%)'),
    color = 'Study sites',
    shape = 'Study sites'
  ) +
  theme(legend.position = 'none',
        #legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12)) #+
  # guides(color = guide_legend(ncol = 4,
  #                             theme = theme(legend.title.position = 'left'
  #                             )))

PCoA_no_clonal_PvGTSeq_LAC



if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_PvGTSeq_LAC.tsv')){
  thres = 1
  
  temp_cluster = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,  
                              threshold = thres,
                              metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                              sample_id = 'Sample_id',
                              levels = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                select(site_of_collection_snl0) %>% unlist() %>% unique(),
                              group_by = 'site_of_collection_snl0',
                              colors = brewer.pal(11, 'Set3'),
                              vertex.size = 4,
                              method = 'fruchtermanreingold')
  
  
  summary_IBDclusters_PvGTSeq_LAC = temp_cluster$clusters %>% 
    summarise(Total = n(),
              nsampClus = sum(grepl('Cluster', Cluster)),
              nSampSing = n() - sum(grepl('Cluster', Cluster)),
              nclusters = length(unique(.$Cluster[grepl('Cluster', Cluster)]))
    )
  
  summary_IBDclusters_PvGTSeq_LAC['Threshold'] = thres
  
  temp_summary_clusters = summary_IBDclusters_PvGTSeq_LAC
  
  while(temp_summary_clusters['nclusters'] > 1# & thres > 0.60
  ){
    thres = thres - 0.01
    
    temp_cluster = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,  
                                threshold = thres,
                                metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                                sample_id = 'Sample_id',
                                levels = monoclonals_PvGTSeq_LAC_loci@metadata %>%
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
    
    summary_IBDclusters_PvGTSeq_LAC = rbind(summary_IBDclusters_PvGTSeq_LAC,
                                        temp_summary_clusters)
    
    print(thres)
    
  }
  
  write.table(summary_IBDclusters_PvGTSeq_LAC, '~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_PvGTSeq_LAC.tsv', quote = F, row.names = F, sep = '\t')
  
}else{
  summary_IBDclusters_PvGTSeq_LAC = read.table('~/Documents/Github/PvGTSeq_paper_draft/Outputs/summary_IBDclusters_PvGTSeq_LAC.tsv', sep = '\t', header = T)  
}


summary_IBDclusters_PvGTSeq_LAC$deriv_nsampClus = summary_IBDclusters_PvGTSeq_LAC$nsampClus
summary_IBDclusters_PvGTSeq_LAC$deriv_nClus = summary_IBDclusters_PvGTSeq_LAC$nclusters


summary_IBDclusters_PvGTSeq_LAC$nclusters
summary_IBDclusters_PvGTSeq_LAC$nSampSing
summary_IBDclusters_PvGTSeq_LAC$nsampClus
#summary_clusters$deriv = NA

summary_IBDclusters_PvGTSeq_LAC[-1,]$deriv_nsampClus = summary_IBDclusters_PvGTSeq_LAC[-1,]$nsampClus - summary_IBDclusters_PvGTSeq_LAC[-nrow(summary_IBDclusters_PvGTSeq_LAC),]$nsampClus
summary_IBDclusters_PvGTSeq_LAC[-1,]$deriv_nClus = summary_IBDclusters_PvGTSeq_LAC[-1,]$nclusters - summary_IBDclusters_PvGTSeq_LAC[-nrow(summary_IBDclusters_PvGTSeq_LAC),]$nclusters

summary_IBDclusters_PvGTSeq_LAC[-nrow(summary_IBDclusters_PvGTSeq_LAC),]$deriv_nClus = summary_IBDclusters_PvGTSeq_LAC[-nrow(summary_IBDclusters_PvGTSeq_LAC),]$nclusters - summary_IBDclusters_PvGTSeq_LAC[-1,]$nclusters


summary_IBDclusters_PvGTSeq_LAC %>% 
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

summary_IBDclusters_PvGTSeq_LAC %>% 
  ggplot(aes(x = Threshold, y = nclusters)) + 
  geom_line(linewidth = 2)


summary_IBDclusters_PvGTSeq_LAC %>% 
  ggplot(aes(x = Threshold, y = deriv_nClus)) + 
  geom_line(linewidth = 2) +
  scale_x_continuous(breaks = seq(0, 1, .02)) + 
  theme(axis.text.x = element_text(angle = 90))



PvGTSeq_LAC_network_100  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 1, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                            filter(!(Sample_id %in% LAC_outliers),
                                                   !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_100$plot_network

PvGTSeq_LAC_network_059  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 0.59, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                            filter(!(Sample_id %in% LAC_outliers),
                                                   !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_059$plot_network

PvGTSeq_LAC_network_050 = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                         threshold = 0.5, 
                                         metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                           filter(!(Sample_id %in% LAC_outliers),
                                                  !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                         sample_id = 'Sample_id',
                                         color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_050$plot_network


PvGTSeq_LAC_network_045 = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 0.45, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                           filter(!(Sample_id %in% LAC_outliers),
                                                  !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_045$plot_network


PvGTSeq_LAC_network_034  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                       threshold = 0.34, 
                                       metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                         filter(!(Sample_id %in% LAC_outliers),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_034$plot_network

PvGTSeq_LAC_network_033  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 0.42, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                            filter(!(Sample_id %in% LAC_outliers),
                                                   !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua', 'USA', 'Trinidad', 'Guatemala', 'FrenchGuiana'))),
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_033$plot_network



PvGTSeq_clusters = plot_network(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,  
                            threshold = .98,
                            metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                            sample_id = 'Sample_id',
                            levels = monoclonals_PvGTSeq_LAC_loci@metadata %>%
                              select(site_of_collection_snl0) %>% unlist() %>% unique(),
                            group_by = 'site_of_collection_snl0',
                            colors = brewer.pal(11, 'Set3'),
                            vertex.size = 4,
                            method = 'fruchtermanreingold')
















