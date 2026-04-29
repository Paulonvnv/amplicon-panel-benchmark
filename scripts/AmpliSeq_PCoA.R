
# PvAmpliSeq ----

PvAmpliSeq_markers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvAmpliSeq_markers_table.csv')

PvAmpliSeq_markers %<>% filter(!(gene_name %in% c('ABCE1',
                                                  'pvcrt_o', 
                                                  'pvdhfr', 
                                                  'pvdhps', 
                                                  'pvdmt2', 
                                                  'pvk13', 
                                                  'pvmdr1', 
                                                  'pvmdr2', 
                                                  'pvmrp1', 
                                                  'pvmrp2', 
                                                  'pvp13k')))



# Extract PvAmpliSeq coordinates from rGenome object----

amplicon = 1

PvAmpliSeq_coordinates = NULL
for(amplicon in 1:nrow(PvAmpliSeq_markers)){
  amplicon_chromosome = PvAmpliSeq_markers[amplicon,][['chromosome']]
  amplicon_start = PvAmpliSeq_markers[amplicon,][['start']]
  amplicon_end = PvAmpliSeq_markers[amplicon,][['end']]
  
  PvAmpliSeq_coordinates = c(PvAmpliSeq_coordinates,
                          paste(amplicon_chromosome,
                                amplicon_start:amplicon_end, sep = '_'))
  
}

PvAmpliSeq_rGenome = read_rGenome('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_broad_Pv4_AmpliSeq_rGenome/', 
                                format = 'tsv',
                                sep = '\t')


sum(duplicated(PvAmpliSeq_rGenome@metadata$Sample_id))
colnames(PvAmpliSeq_rGenome@gt) = PvAmpliSeq_rGenome@metadata$Sample_id

sum(rownames(PvAmpliSeq_rGenome@metadata) != colnames(PvAmpliSeq_rGenome@gt))

PvAmpliSeq_rGenome@metadata = left_join(PvAmpliSeq_rGenome@metadata,
                                        All_broad_Pv4_PvGTSeq_metadata,
                                        by = 'Sample_id')

PvAmpliSeq_rGenome = filter_samples(PvAmpliSeq_rGenome,
                                    !(duplicated(PvAmpliSeq_rGenome@metadata$PvGTSeq_id)))

sum(duplicated(PvAmpliSeq_rGenome@metadata$PvGTSeq_id))


PvAmpliSeq_rGenome@metadata %<>% mutate(
  Old_Sample_id = Sample_id,
  Sample_id = PvGTSeq_id
)

colnames(PvAmpliSeq_rGenome@gt) = PvAmpliSeq_rGenome@metadata$Sample_id
rownames(PvAmpliSeq_rGenome@metadata) = PvAmpliSeq_rGenome@metadata$Sample_id


names(PvAmpliSeq_rGenome@metadata)


PvAmpliSeq_rGenome@metadata %<>% mutate(
  site_of_collection_world_region = 
    case_when(grepl('Latin', possible_site_of_infection_world_region) ~ 'LAC',
              site_of_collection_snl0  == 'Thailand' ~ 'WSEA',
              site_of_collection_snl0 == 'China' ~ 'EAS',
              site_of_collection_snl0 == 'Bangladesh' ~ 'WSEA',
              site_of_collection_snl0 == 'Bhutan' ~ 'WSEA',
              site_of_collection_snl0 == 'Madagascar' ~ 'AF',
              site_of_collection_snl0 == 'Mauritania' ~ 'AF',
              site_of_collection_snl0 == 'Myanmar' ~ 'WSEA',
              site_of_collection_snl0 == 'North Korea' ~ 'EAS',
              site_of_collection_snl0 == 'Sudan' ~ 'AF',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Panama' ~ 'LAC',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Venezuela' ~ 'LAC',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Guyana' ~ 'LAC',
              possible_site_of_infection_world_region == 'West Africa' ~ 'AF',
              possible_site_of_infection_world_region == 'Asia' ~ 'WAS',
              .default = possible_site_of_infection_world_region))

dim(PvAmpliSeq_rGenome@gt)

PvAmpliSeq_rGenome = filter_loci(PvAmpliSeq_rGenome, v = rownames(PvAmpliSeq_rGenome@gt) %in% PvAmpliSeq_coordinates)

dim(PvAmpliSeq_rGenome@gt)

# Update allele frequencies and remove monomorphic sites----

PvAmpliSeq_rGenome = update_allele_lables(PvAmpliSeq_rGenome, n = 10)

## Remove monomorphic sites----
PvAmpliSeq_rGenome = filter_loci(PvAmpliSeq_rGenome, 
                                 v = PvAmpliSeq_rGenome@loci_table$Cardinality > 1)

dim(PvAmpliSeq_rGenome@gt)

# Calculate amplification rate per Loci-----
PvAmpliSeq_rGenome@loci_table$Locus_Ampl_Rate = LocusAmplRate(PvAmpliSeq_rGenome, update = F)

min(PvAmpliSeq_rGenome@loci_table$Locus_Ampl_Rate)

## Remove loci with less than 50% of genome coverage----
PvAmpliSeq_rGenome = filter_loci(PvAmpliSeq_rGenome, 
                                 v = PvAmpliSeq_rGenome@loci_table$Locus_Ampl_Rate >= .75)

dim(PvAmpliSeq_rGenome@gt)

names(PvAmpliSeq_rGenome@metadata)

# Identify the type of polymophism of the variant site----
PvAmpliSeq_rGenome@loci_table$type_of_polymorphism = get_type_of_polymorphism(PvAmpliSeq_rGenome)

unique(PvAmpliSeq_rGenome@loci_table$type_of_polymorphism)

## Remove Homopolymers and STRs----
PvAmpliSeq_rGenome = 
  filter_loci(PvAmpliSeq_rGenome,
              v = !(PvAmpliSeq_rGenome@loci_table$type_of_polymorphism %in% 
                      c('INDEL:Homopolymer', 'INDEL:Dinucleotide_STR')))

dim(PvAmpliSeq_rGenome@gt)

# Keep monoclonal samples with 75% of wgs covirage----

length(monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$Sample_id)


monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$Sample_id

monoclonals_PvAmpliSeq_rGenome = 
  filter_samples(PvAmpliSeq_rGenome,
  PvAmpliSeq_rGenome@metadata$Sample_id %in% monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$Sample_id)

dim(monoclonals_PvAmpliSeq_rGenome@gt)


# Update allele frequencies and remove monomorphic sites----

monoclonals_PvAmpliSeq_rGenome = update_allele_lables(monoclonals_PvAmpliSeq_rGenome, n = 1)

## Remove monomorphic sites----

sum(monoclonals_PvAmpliSeq_rGenome@loci_table$Cardinality == 1)

monoclonals_PvAmpliSeq_rGenome = filter_loci(monoclonals_PvAmpliSeq_rGenome, 
                                 v = monoclonals_PvAmpliSeq_rGenome@loci_table$Cardinality > 1)

dim(monoclonals_PvAmpliSeq_rGenome@gt)


# Calculate minor allele frequency
PvAmpliSeq_major_freq = sapply(monoclonals_PvAmpliSeq_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvAmpliSeq_major_freq)

monoclonals_PvAmpliSeq_rGenome@loci_table$major_freq = PvAmpliSeq_major_freq

monoclonals_PvAmpliSeq_rGenome = filter_loci(monoclonals_PvAmpliSeq_rGenome, 
                                 v = monoclonals_PvAmpliSeq_rGenome@loci_table$major_freq <= 1 - 5/996)

dim(monoclonals_PvAmpliSeq_rGenome@gt)

monoclonals_PvAmpliSeq_rGenome@metadata$sample_ampl_rate = SampleAmplRate(monoclonals_PvAmpliSeq_rGenome, update = F)

monoclonals_PvAmpliSeq_rGenome@metadata %>% 
  ggplot(aes(x = sample_ampl_rate)) +
  geom_histogram()

# Create loci object----

monoclonals_PvAmpliSeq_gt = gsub('_\\d+', '', gsub(':\\d+', '', monoclonals_PvAmpliSeq_rGenome@gt))

monoclonals_PvAmpliSeq_gt = t(monoclonals_PvAmpliSeq_gt)

monoclonals_PvAmpliSeq_loci = create_loci(loci_table = monoclonals_PvAmpliSeq_gt,
                                                                metadata = monoclonals_PvAmpliSeq_rGenome@metadata)

dim(monoclonals_PvAmpliSeq_loci@loci_table)


PvAmpliSeq_allele_freqs = get_allele_freq(monoclonals_PvAmpliSeq_loci)

monoclonals_PvAmpliSeq_loci@freq_table = PvAmpliSeq_allele_freqs

nrow(PvAmpliSeq_allele_freqs)
dim(monoclonals_PvAmpliSeq_loci@loci_table)


monoclonals_PvAmpliSeq_loci@markers = 
  monoclonals_PvAmpliSeq_rGenome@loci_table %>% rename('CHROM' = 'chromosome',
                                                       'POS' = 'pos')


rownames(monoclonals_PvAmpliSeq_loci@loci_table) = monoclonals_PvAmpliSeq_loci@metadata$Sample_id

rownames(monoclonals_PvAmpliSeq_loci@loci_table)

# Monoclonal pairs by country----

PvAmpliSeq_Monoclonal_pairs = as.data.frame(t(combn(Monoclonals, 2)))

names(PvAmpliSeq_Monoclonal_pairs) = c('Yi', 'Yj')

PvAmpliSeq_Monoclonal_pairs = left_join(
  PvAmpliSeq_Monoclonal_pairs,
  monoclonals_PvAmpliSeq_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(PvAmpliSeq_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


PvAmpliSeq_Monoclonal_pairs = left_join(
  PvAmpliSeq_Monoclonal_pairs,
  monoclonals_PvAmpliSeq_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(PvAmpliSeq_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')


# Calculate IBD ---

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:1000){
    start_time = Sys.time()
    Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvAmpliSeq_loci, parallel = T, pairs = Monoclonal_pairs, max_k = 20,
              w = w, n = 1000
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %<>%
    mutate(euDist = 1 - rhat)
  
  Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = 
    left_join(Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD,
              Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #sum(Monoclonal_LAC_pairs_01SNPsbiallelic_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv')
  
}



names(Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD)


Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')



evectors_PvAmpliSeq_world_regions = GRM_evectors(dist_table = Monoclonal_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD, 
                                                 k = nrow(monoclonals_PvAmpliSeq_loci@metadata %>%
                                                            filter(Sample_id %in% selected_samples_ids, 
                                                                   site_of_collection_world_region != 'Pv4')
                                                 ),
                                                 metadata = monoclonals_PvAmpliSeq_loci@metadata%>%
                                                   filter(Sample_id %in% selected_samples_ids, site_of_collection_world_region != 'Pv4'),
                                                 Pop = 'site_of_collection_world_region',
                                                 method = 'princomp'
)

monoclonals_PvAmpliSeq_loci@markers$amplicon = sapply(1:nrow(monoclonals_PvAmpliSeq_loci@markers), function(i){
  
  variant_site_chromosome = monoclonals_PvAmpliSeq_loci@markers[i, ][['chromosome']]
  variant_site_position = monoclonals_PvAmpliSeq_loci@markers[i, ][['pos']]
  
  PvAmpliSeq_markers[PvAmpliSeq_markers$chromosome == variant_site_chromosome &
                       PvAmpliSeq_markers$start <= variant_site_position &
                       PvAmpliSeq_markers$end >= variant_site_position,
  ][['amplicon']]
  
})

nrow(PvAmpliSeq_markers)

length(unique(monoclonals_PvAmpliSeq_loci@markers$amplicon))

PCoA_no_clonal_PvAmpliSeq_world_regions = evectors_PvAmpliSeq_world_regions$eigenvector %>%
  ggplot(aes(x = -PC1, 
             y = PC2, 
             color = site_of_collection_world_region
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('C) World regions using ', 
                      nrow(monoclonals_PvAmpliSeq_loci@markers %>% filter(type_of_polymorphism == 'SNP')), 
                      ' SNPs from\n',
                      length(unique(monoclonals_PvAmpliSeq_loci@markers$amplicon)),
                      ' partial amplicons from PvAmpliSeq'),
    x = paste0('PCoA 1 (', round(evectors_PvAmpliSeq_world_regions$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvAmpliSeq_world_regions$contrib[2], 2), '%)'),
    color = 'Study sites',
    shape = 'Study sites'
  ) +
  theme(legend.position = 'bottom',
        legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12)) +
  guides(color = guide_legend(ncol = 4,
                              theme = theme(legend.title.position = 'left'
                              )))

PCoA_no_clonal_PvAmpliSeq_world_regions



# Calculate IBD for LAC---

## Create loci object----

dim(monoclonals_PvAmpliSeq_rGenome@gt) 

monoclonals_PvAmpliSeq_LAC_rGenome = 
  filter_samples(monoclonals_PvAmpliSeq_rGenome,
                 monoclonals_PvAmpliSeq_rGenome@metadata$site_of_collection_world_region == 'LAC'
  )

monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Cardinality = NULL
monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Alleles = NULL
monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Allele_Counts = NULL
monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$major_freq = NULL

monoclonals_PvAmpliSeq_LAC_rGenome = update_allele_lables(monoclonals_PvAmpliSeq_LAC_rGenome)

sum(monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Cardinality == 2)


monoclonals_PvAmpliSeq_LAC_rGenome = filter_loci(monoclonals_PvAmpliSeq_LAC_rGenome,
                                                 monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Cardinality >= 2)

dim(monoclonals_PvAmpliSeq_LAC_rGenome@gt)


# Calculate minor allele frequency
PvAmpliSeq_major_freq = sapply(monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvAmpliSeq_major_freq)
1 - 5/488

sum(PvAmpliSeq_major_freq <= 1 - 5/488)

monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$major_freq = PvAmpliSeq_major_freq

monoclonals_PvAmpliSeq_LAC_rGenome = filter_loci(monoclonals_PvAmpliSeq_LAC_rGenome,
                                              v = monoclonals_PvAmpliSeq_LAC_rGenome@loci_table$major_freq <= 1 - 5/488)

dim(monoclonals_PvAmpliSeq_LAC_rGenome@gt)

# sum(monoclonals_PvGTSeq_LAC_rGenome@loci_table$type_of_polymorphism != 'SNP')
# 
# monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
#                                               monoclonals_PvGTSeq_LAC_rGenome@loci_table$type_of_polymorphism == 'SNP')

dim(monoclonals_PvAmpliSeq_LAC_rGenome@gt)
View(monoclonals_PvAmpliSeq_LAC_rGenome@gt)

#monoclonals_PvGTSeq_gt = gsub(':\\d+', '', monoclonals_PvGTSeq_LAC_rGenome@gt)
monoclonals_PvAmpliSeq_gt = gsub('/\\d+', '', gsub(':\\d+', '', monoclonals_PvAmpliSeq_LAC_rGenome@gt))

monoclonals_PvAmpliSeq_gt = t(monoclonals_PvAmpliSeq_gt)

monoclonals_PvAmpliSeq_LAC_loci = create_loci(loci_table = monoclonals_PvAmpliSeq_gt,
                                           metadata = monoclonals_PvAmpliSeq_LAC_rGenome@metadata)

dim(monoclonals_PvAmpliSeq_LAC_loci@loci_table)


# monoclonals_PvGTSeq_LAC_loci@metadata %<>%
#   mutate(Strata = case_when(
#     site_of_collection_snl0 %in% c('Colombia') ~ 'Colombia',
#     site_of_collection_snl0 %in% c('El Salvador', 'Honduras', 'Mexico', 'Nicaragua', 'Panama') ~ 'Meso America',
#     site_of_collection_snl0 %in% c('Brazil', 'Peru') ~ 'Amazon',
#     site_of_collection_snl0 %in% c('Guyana', 'Venezuela') ~ 'Guainia Shield'
#   ))

monoclonals_PvGTSeq_LAC_loci_no_clonal = filter_samples(monoclonals_PvAmpliSeq_LAC_loci,
                                                        monoclonals_PvAmpliSeq_LAC_loci@metadata$Sample_id %in% selected_LAC_samples_ids
)

PvAmpliSeq_allele_freqs = get_allele_freq(monoclonals_PvAmpliSeq_LAC_loci#, by = 'Strata'
)

monoclonals_PvAmpliSeq_LAC_loci@freq_table = PvAmpliSeq_allele_freqs


dim(monoclonals_PvAmpliSeq_LAC_loci@loci_table)


monoclonals_PvAmpliSeq_LAC_loci@markers = 
  monoclonals_PvAmpliSeq_LAC_rGenome@loci_table %>% rename('CHROM' = 'chromosome',
                                                        'POS' = 'pos')


rownames(monoclonals_PvAmpliSeq_LAC_loci@loci_table) = monoclonals_PvAmpliSeq_LAC_loci@metadata$Sample_id

rownames(monoclonals_PvAmpliSeq_LAC_loci@loci_table)

# Monoclonal pairs ----

Monoclonals = monoclonals_PvAmpliSeq_LAC_loci@metadata$Sample_id

PvAmpliSeq_LAC_Monoclonal_pairs = as.data.frame(t(combn(Monoclonals, 2)))

names(PvAmpliSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj')

PvAmpliSeq_LAC_Monoclonal_pairs = left_join(
  PvAmpliSeq_LAC_Monoclonal_pairs,
  monoclonals_PvAmpliSeq_LAC_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yi' == 'Sample_id'))

names(PvAmpliSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0')


PvAmpliSeq_LAC_Monoclonal_pairs = left_join(
  PvAmpliSeq_LAC_Monoclonal_pairs,
  monoclonals_PvAmpliSeq_LAC_loci@metadata[,c('Sample_id', 'site_of_collection_snl0')],
  by = join_by('Yj' == 'Sample_id'))

names(PvAmpliSeq_LAC_Monoclonal_pairs) = c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:500){
    start_time = Sys.time()
    Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvAmpliSeq_LAC_loci, parallel = T, pairs = PvAmpliSeq_LAC_Monoclonal_pairs, max_k = 20,
              w = w, n = 500
            ))
    end_time = Sys.time()
    print(paste0('Window ', w, ' done in:'))
    print(end_time - start_time)
  }
  
  Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %<>%
    mutate(euDist = 1 - rhat)
  
  Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = 
    left_join(Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD,
              PvAmpliSeq_LAC_Monoclonal_pairs[, c('Yi', 'Yj', 'Yi_site_of_collection_snl0', 'Yj_site_of_collection_snl0')], by = c('Yi', 'Yj'))
  
  #sum(Monoclonal_LAC_pairs_01SNPsbiallelic_75_dist$site_of_collection_snl0 %in% c(removed_countries, 'Pv4'))
  
  write.csv(Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD, 
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv',
            quote = F,
            row.names = F)
  
}else{
  
  Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD.csv')
  
}



names(Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD)


Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %>%
  ggplot(aes(x = rhat)) +
  geom_histogram()+
  geom_vline(xintercept = .99)


Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD %>%
  filter(Yi_site_of_collection_snl0 == Yj_site_of_collection_snl0) %>%
  ggplot(aes(x = rhat, color = Yi_site_of_collection_snl0, fill = Yi_site_of_collection_snl0)) +
  geom_density()+
  geom_vline(xintercept = .99)+
  facet_grid(Yi_site_of_collection_snl0 ~., scales = 'free_y')



evectors_PvAmpliSeq_LAC = GRM_evectors(dist_table = Monoclonal_LAC_pairs_PvAmpliSeq_01SNPsbiallelic_75_IBD, 
                                                 k = nrow(monoclonals_PvAmpliSeq_LAC_loci@metadata %>%
                                                            filter(Sample_id %in% selected_LAC_samples_ids,
                                                                   !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua')))
                                                 ),
                                                 metadata = monoclonals_PvAmpliSeq_LAC_loci@metadata%>%
                                                   filter(Sample_id %in% selected_LAC_samples_ids, 
                                                          !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua'))),
                                                 Pop = 'site_of_collection_snl0',
                                                 method = 'princomp'
)

monoclonals_PvAmpliSeq_LAC_loci@markers$amplicon = sapply(1:nrow(monoclonals_PvAmpliSeq_LAC_loci@markers), function(i){
  
  variant_site_chromosome = monoclonals_PvAmpliSeq_LAC_loci@markers[i, ][['chromosome']]
  variant_site_position = monoclonals_PvAmpliSeq_LAC_loci@markers[i, ][['pos']]
  
  PvAmpliSeq_markers[PvAmpliSeq_markers$chromosome == variant_site_chromosome &
                       PvAmpliSeq_markers$start <= variant_site_position &
                       PvAmpliSeq_markers$end >= variant_site_position,
  ][['amplicon']]
  
})

nrow(PvAmpliSeq_markers)

length(unique(monoclonals_PvAmpliSeq_LAC_loci@markers$amplicon))


PCoA_no_clonal_PvAmpliSeq_LAC = evectors_PvAmpliSeq_LAC$eigenvector %>%
  ggplot(aes(x = -PC1, 
             y = -PC2, 
             color = site_of_collection_snl0
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('C) LAC countries using ', 
                      nrow(monoclonals_PvAmpliSeq_LAC_loci@markers %>% filter(type_of_polymorphism == 'SNP')), 
                      ' SNPs from\n',
                      length(unique(monoclonals_PvAmpliSeq_LAC_loci@markers$amplicon)),
                      ' partial amplicons from PvAmpliSeq'),
    x = paste0('PCoA 1 (', round(evectors_PvAmpliSeq_LAC$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvAmpliSeq_LAC$contrib[2], 2), '%)'),
    color = 'Study sites',
    shape = 'Study sites'
  ) +
  theme(legend.position = 'none',
        legend.justification = 'left',
        text = element_text(size = 12),
        axis.text = element_text(size = 12),
        title = element_text(size = 12)) #+
  # guides(color = guide_legend(ncol = 4,
  #                             theme = theme(legend.title.position = 'left'
  #                             )))

PCoA_no_clonal_PvAmpliSeq_LAC





PvAmpliSeq_ampseq = 
  read_ampseq(file = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_broad_Pv4_AmpliSeq_ampseq/', 
              format = 'tsv',
              sep = '\t')




PvAmpliSeq_ampseq@metadata %<>% mutate(
  site_of_collection_world_region = 
    case_when(grepl('Latin', possible_site_of_infection_world_region) ~ 'LAC',
              site_of_collection_snl0  == 'Thailand' ~ 'WSEA',
              site_of_collection_snl0 == 'China' ~ 'EAS',
              site_of_collection_snl0 == 'Bangladesh' ~ 'WSEA',
              site_of_collection_snl0 == 'Bhutan' ~ 'WSEA',
              site_of_collection_snl0 == 'Madagascar' ~ 'AF',
              site_of_collection_snl0 == 'Mauritania' ~ 'AF',
              site_of_collection_snl0 == 'Myanmar' ~ 'WSEA',
              site_of_collection_snl0 == 'North Korea' ~ 'EAS',
              site_of_collection_snl0 == 'Sudan' ~ 'AF',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Panama' ~ 'LAC',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Venezuela' ~ 'LAC',
              possible_site_of_infection_world_region == '' & site_of_collection_snl0 == 'Guyana' ~ 'LAC',
              possible_site_of_infection_world_region == 'West Africa' ~ 'AF',
              possible_site_of_infection_world_region == 'Asia' ~ 'WAS',
              .default = possible_site_of_infection_world_region))


PvAmpliSeq_ampseq = filter_samples(PvAmpliSeq_ampseq, PvAmpliSeq_ampseq@metadata$site_of_collection_world_region != 'Pv4')

PvAmpliSeq_ampseq@markers$pos = round(PvAmpliSeq_ampseq@markers$end - PvAmpliSeq_ampseq@markers$start)

PvAmpliSeq_ampseq = filter_loci(obj = PvAmpliSeq_ampseq, v = PvAmpliSeq_ampseq@markers$amplicon %in% PvAmpliSeq_markers$amplicon, skip_errors = T)

PvAmpliSeq_coverage = get_ReadDepth_coverage(PvAmpliSeq_ampseq, variable = 'site_of_collection_world_region')

PvAmpliSeq_coverage$plot_read_depth_heatmap

PvAmpliSeq_locusAmpRate = locus_amplification_rate(PvAmpliSeq_ampseq, 
                                                threshold = .75,
                                                strata = 'site_of_collection_world_region', update_loci = F)


PvAmpliSeq_locusAmpRate$all_loci_performance_plot +
  labs(title = 'B) PvAmpliSeq',
       y = '# of Amplicons',
       x = 'Propottion of amplified samples')




