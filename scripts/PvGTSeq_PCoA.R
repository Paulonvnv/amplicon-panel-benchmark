

# PvGTSeq ----

PvGTSeq_markers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_markersTable.csv')

PvGTSeq_markers %<>% filter(use !='DRS')

# Extract PvGTSeq coordinates from rGenome object----

amplicon = 1

PvGTSeq_coordinates = NULL
for(amplicon in 1:nrow(PvGTSeq_markers)){
  amplicon_chromosome = PvGTSeq_markers[amplicon,][['chromosome']]
  amplicon_start = PvGTSeq_markers[amplicon,][['start']]
  amplicon_end = PvGTSeq_markers[amplicon,][['end']]
  
  PvGTSeq_coordinates = c(PvGTSeq_coordinates,
                             paste(amplicon_chromosome,
                                   amplicon_start:amplicon_end, sep = '_'))
  
}

PvGTSeq_rGenome = read_rGenome('~/Documents/Github/PvDemography_SAmer/PvWGS/all_broad_WGS_PvGTSeq_rGenome/', 
                                format = 'tsv',
                                sep = '\t')


sum(duplicated(PvGTSeq_rGenome@metadata$Sample_id))
colnames(PvGTSeq_rGenome@gt) = PvGTSeq_rGenome@metadata$Sample_id

sum(rownames(PvGTSeq_rGenome@metadata) != colnames(PvGTSeq_rGenome@gt))


All_broad_Pv4_PvGTSeq_metadata = read.table('~/Documents/Github/DataManagment_NeafseyLab/Metadata/All_Broad_PvSamples/All_broad_Pv4_PvWGS_metadata.tsv', sep = '\t', header = T)


PvGTSeq_rGenome@metadata = left_join(PvGTSeq_rGenome@metadata,
                                        All_broad_Pv4_PvGTSeq_metadata,
                                        by = 'Sample_id')

sum(duplicated(PvGTSeq_rGenome@metadata$PvGTSeq_id))

PvGTSeq_rGenome = filter_samples(PvGTSeq_rGenome,
                                    !(duplicated(PvGTSeq_rGenome@metadata$PvGTSeq_id)))

sum(duplicated(PvGTSeq_rGenome@metadata$PvGTSeq_id))

monoclonals = monoclonals_PvBroad_SelectedPos2_01SNPsbiallelic_75_loci@metadata$Sample_id

sum(!(monoclonals %in% PvGTSeq_rGenome@metadata$PvGTSeq_id))

monoclonals[!(monoclonals %in% PvGTSeq_rGenome@metadata$Sample_id)]

monoclonals[!(monoclonals %in% PvGTSeq_rGenome@metadata$PvGTSeq_id)]

PvGTSeq_rGenome@metadata %<>% mutate(
  Old_Sample_id = Sample_id,
  Sample_id = PvGTSeq_id
)

colnames(PvGTSeq_rGenome@gt) = PvGTSeq_rGenome@metadata$Sample_id
rownames(PvGTSeq_rGenome@metadata) = PvGTSeq_rGenome@metadata$Sample_id


names(PvGTSeq_rGenome@metadata)


PvGTSeq_rGenome@metadata %<>% mutate(
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

dim(PvGTSeq_rGenome@gt)

PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, v = rownames(PvGTSeq_rGenome@gt) %in% PvGTSeq_coordinates)

dim(PvGTSeq_rGenome@gt)

# Update allele frequencies and remove monomorphic sites----

PvGTSeq_rGenome = update_allele_lables(PvGTSeq_rGenome, n = 10)

## Remove monomorphic sites----
PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, 
                                 v = PvGTSeq_rGenome@loci_table$Cardinality > 1)

dim(PvGTSeq_rGenome@gt)

# Calculate amplification rate per Loci-----
PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate = LocusAmplRate(PvGTSeq_rGenome, update = F)

min(PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate)

## Remove loci with less than 50% of genome coverage----
PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, 
                                 v = PvGTSeq_rGenome@loci_table$Locus_Ampl_Rate >= .75)

dim(PvGTSeq_rGenome@gt)

names(PvGTSeq_rGenome@metadata)

# Identify the type of polymophism of the variant site----
PvGTSeq_rGenome@loci_table$type_of_polymorphism = get_type_of_polymorphism(PvGTSeq_rGenome)

unique(PvGTSeq_rGenome@loci_table$type_of_polymorphism)

## Remove Homopolymers and STRs----
PvGTSeq_rGenome = 
  filter_loci(PvGTSeq_rGenome,
              v = !(PvGTSeq_rGenome@loci_table$type_of_polymorphism %in% 
                      c('INDEL:Homopolymer', 'INDEL:Dinucleotide_STR')))

dim(PvGTSeq_rGenome@gt)

# Keep monoclonal samples with 75% of wgs covirage----

length(monoclonals_PvBroad_SelectedPos2_rGenome_01SNPsbiallelic_75@metadata$Sample_id)


monoclonals_PvGTSeq_rGenome = 
  filter_samples(PvGTSeq_rGenome,
                 PvGTSeq_rGenome@metadata$Sample_id %in% monoclonals)

dim(monoclonals_PvGTSeq_rGenome@gt)


# Update allele frequencies and remove monomorphic sites----

monoclonals_PvGTSeq_rGenome = update_allele_lables(monoclonals_PvGTSeq_rGenome, n = 1)

## Remove monomorphic sites----

sum(monoclonals_PvGTSeq_rGenome@loci_table$Cardinality == 1)

monoclonals_PvGTSeq_rGenome = filter_loci(monoclonals_PvGTSeq_rGenome, 
                                             v = monoclonals_PvGTSeq_rGenome@loci_table$Cardinality > 1)

dim(monoclonals_PvGTSeq_rGenome@gt)


# Calculate minor allele frequency
PvGTSeq_major_freq = sapply(monoclonals_PvGTSeq_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvGTSeq_major_freq)

monoclonals_PvGTSeq_rGenome@loci_table$major_freq = PvGTSeq_major_freq

sum(monoclonals_PvGTSeq_rGenome@loci_table$major_freq <= 1 - 5/1023)

monoclonals_PvGTSeq_rGenome = filter_loci(monoclonals_PvGTSeq_rGenome, 
                                             v = monoclonals_PvGTSeq_rGenome@loci_table$major_freq <= 1 - 5/1023)

dim(monoclonals_PvGTSeq_rGenome@gt)

monoclonals_PvGTSeq_rGenome@metadata$sample_ampl_rate = SampleAmplRate(monoclonals_PvGTSeq_rGenome, update = F)

monoclonals_PvGTSeq_rGenome@metadata %>% 
  ggplot(aes(x = sample_ampl_rate)) +
  geom_histogram()

monoclonals_PvGTSeq_rGenome = filter_samples(monoclonals_PvGTSeq_rGenome,
                                             monoclonals_PvGTSeq_rGenome@metadata$sample_ampl_rate >= .75)

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
  monoclonals_PvGTSeq_rGenome@loci_table %>% rename('CHROM' = 'chromosome',
                                                       'POS' = 'pos')


rownames(monoclonals_PvGTSeq_loci@loci_table) = monoclonals_PvGTSeq_loci@metadata$Sample_id

rownames(monoclonals_PvGTSeq_loci@loci_table)

# Monoclonal pairs ----

Monoclonals = monoclonals_PvGTSeq_loci@metadata$Sample_id

PvGTSeq_Monoclonal_pairs = as.data.frame(t(combn(Monoclonals, 2)))

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



# Calculate IBD for world regions---

if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:1000){
    start_time = Sys.time()
    Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvGTSeq_loci, parallel = T, pairs = PvGTSeq_Monoclonal_pairs, max_k = 20,
              w = w, n = 1000
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
                                                                   site_of_collection_world_region != 'Pv4')
                                                 ),
                                                 metadata = monoclonals_PvGTSeq_loci@metadata%>%
                                                   filter(Sample_id %in% selected_samples_ids, site_of_collection_world_region != 'Pv4'),
                                                 Pop = 'site_of_collection_world_region',
                                                 method = 'princomp'
)

monoclonals_PvGTSeq_loci@markers$amplicon = sapply(1:nrow(monoclonals_PvGTSeq_loci@markers), function(i){
  
  variant_site_chromosome = monoclonals_PvGTSeq_loci@markers[i, ][['chromosome']]
  variant_site_position = monoclonals_PvGTSeq_loci@markers[i, ][['pos']]
  
  PvGTSeq_markers[PvGTSeq_markers$chromosome == variant_site_chromosome &
                     PvGTSeq_markers$start <= variant_site_position &
                     PvGTSeq_markers$end >= variant_site_position,
  ][['amplicon']]
  
})

nrow(PvGTSeq_markers)

length(unique(monoclonals_PvGTSeq_loci@markers$amplicon))

monoclonals_PvGTSeq_loci@markers %>% View

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
                      nrow(monoclonals_PvGTSeq_loci@markers %>% filter(type_of_polymorphism == 'SNP')), 
                      ' SNPs from\n',
                      length(unique(monoclonals_PvGTSeq_loci@markers$amplicon)),
                      ' partial amplicons from PvGTSeq'),
    x = paste0('PCoA 1 (', round(evectors_PvGTSeq_world_regions$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvGTSeq_world_regions$contrib[2], 2), '%)'),
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

PCoA_no_clonal_PvGTSeq_world_regions


evectors_PvGTSeq_world_regions$eigenvector %>% 
  filter(PC2 > 1 & site_of_collection_world_region == 'LAC') %>%
  select(Sample_id)


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

sum(monoclonals_PvGTSeq_LAC_rGenome@loci_table$Cardinality == 2)


monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
            monoclonals_PvGTSeq_LAC_rGenome@loci_table$Cardinality >= 2)

dim(monoclonals_PvGTSeq_LAC_rGenome@gt)


# Calculate minor allele frequency
PvGTSeq_major_freq = sapply(monoclonals_PvGTSeq_LAC_rGenome@loci_table$Allele_Counts, function(allele_counts){
  counts = as.numeric(unlist(str_split(gsub('\\d+:', '', allele_counts), ',')))
  freqs = counts/sum(counts)
  major_freq = max(freqs)
})



max(PvGTSeq_major_freq)
1 - 5/482

sum(PvGTSeq_major_freq <= 1 - 5/482)

monoclonals_PvGTSeq_LAC_rGenome@loci_table$major_freq = PvGTSeq_major_freq

monoclonals_PvGTSeq_LAC_rGenome = filter_loci(monoclonals_PvGTSeq_LAC_rGenome,
                                             v = monoclonals_PvGTSeq_LAC_rGenome@loci_table$major_freq <= 1 - 5/482)

dim(monoclonals_PvGTSeq_LAC_rGenome@gt)

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

monoclonals_PvGTSeq_LAC_loci_no_clonal = filter_samples(monoclonals_PvGTSeq_LAC_loci,
                                                        monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id %in% selected_LAC_samples_ids
                                                        )

PvGTSeq_allele_freqs = get_allele_freq(monoclonals_PvGTSeq_LAC_loci#, by = 'Strata'
                                       )

monoclonals_PvGTSeq_LAC_loci@freq_table = PvGTSeq_allele_freqs

ncol(PvGTSeq_allele_freqs$Amazon)
dim(monoclonals_PvGTSeq_LAC_loci@loci_table)


monoclonals_PvGTSeq_LAC_loci@markers = 
  monoclonals_PvGTSeq_LAC_rGenome@loci_table %>% rename('CHROM' = 'chromosome',
                                                    'POS' = 'pos')


rownames(monoclonals_PvGTSeq_LAC_loci@loci_table) = monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id

rownames(monoclonals_PvGTSeq_LAC_loci@loci_table)

# Monoclonal pairs ----

Monoclonals = monoclonals_PvGTSeq_LAC_loci@metadata$Sample_id

PvGTSeq_LAC_Monoclonal_pairs = as.data.frame(t(combn(Monoclonals, 2)))

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



if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD.csv')){
  
  Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = NULL
  for(w in 1:500){
    start_time = Sys.time()
    Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD = 
      rbind(Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
            pairwise_hmmIBD(
              obj = monoclonals_PvGTSeq_LAC_loci, parallel = T, pairs = PvGTSeq_LAC_Monoclonal_pairs, max_k = 20,
              #freq_table = monoclonals_PvGTSeq_LAC_loci@freq_table, by = 'Strata',
              w = w, n = 500
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

evectors_PvGTSeq_LAC$eigenvector %>%
  filter(PC2 > 1.5) %>%
  select(Sample_id)

evectors_PvGTSeq_LAC = GRM_evectors(dist_table = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD, 
                                       k = nrow(monoclonals_PvGTSeq_LAC_loci@metadata %>%
                                                  filter(Sample_id %in% selected_LAC_samples_ids,
                                                         #!(Sample_id %in% c('ERR2309697', 'ERR2678957', 'ERR2678961')),
                                                         #!(Sample_id %in% c('ERR2351947', 'ERR2351948', 'ERR2678958', 'ERR2678959', 'ERR2678960')),
                                                         !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua')))
                                       ),
                                       metadata = monoclonals_PvGTSeq_LAC_loci@metadata%>%
                                         filter(Sample_id %in% selected_LAC_samples_ids, 
                                                #!(Sample_id %in% c('ERR2309697', 'ERR2678957', 'ERR2678961')),
                                                #!(Sample_id %in% c('ERR2351947', 'ERR2351948', 'ERR2678958', 'ERR2678959', 'ERR2678960')),
                                                !(site_of_collection_snl0 %in% c('El Salvador', 'Nicaragua'))),
                                       Pop = 'site_of_collection_snl0',
                                       method = 'princomp'
)

monoclonals_PvGTSeq_LAC_loci@markers$amplicon = sapply(1:nrow(monoclonals_PvGTSeq_LAC_loci@markers), function(i){
  
  variant_site_chromosome = monoclonals_PvGTSeq_LAC_loci@markers[i, ][['chromosome']]
  variant_site_position = monoclonals_PvGTSeq_LAC_loci@markers[i, ][['pos']]
  
  PvGTSeq_markers[PvGTSeq_markers$chromosome == variant_site_chromosome &
                    PvGTSeq_markers$start <= variant_site_position &
                    PvGTSeq_markers$end >= variant_site_position,
  ][['amplicon']]
  
})

length(unique(monoclonals_PvGTSeq_LAC_loci@markers$amplicon))

PCoA_no_clonal_PvGTSeq_LAC = evectors_PvGTSeq_LAC$eigenvector %>%
  ggplot(aes(x = -PC1, 
             y = -PC2, 
             color = site_of_collection_snl0
  ))+
  geom_point(size = 2, alpha = .85)+
  scale_color_brewer(
    palette = 'Paired')+
  theme_minimal()+
  labs(title = paste0('B) LAC countries using ', 
                      nrow(monoclonals_PvGTSeq_LAC_loci@markers %>% filter(type_of_polymorphism == 'SNP')), 
                      ' SNPs from\n',
                      length(unique(monoclonals_PvGTSeq_LAC_loci@markers$amplicon)),
                      ' partial amplicons from PvGTSeq'),
    x = paste0('PCoA 1 (', round(evectors_PvGTSeq_LAC$contrib[1], 2), '%)'),
    y = paste0('PCoA 2 (', round(evectors_PvGTSeq_LAC$contrib[2], 2), '%)'),
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

PCoA_no_clonal_PvGTSeq_LAC




PvGTSeq_ampseq = 
  read_ampseq(file = '~/Documents/Github/PvDemography_SAmer/PvWGS/all_broad_WGS_PvGTSeq_ampseq', 
              format = 'tsv',
              sep = '\t'
              )


PvGTSeq_ampseq@metadata %<>% mutate(
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


PvGTSeq_ampseq = filter_samples(PvGTSeq_ampseq, PvGTSeq_ampseq@metadata$site_of_collection_world_region != 'Pv4')

PvGTSeq_ampseq = filter_loci(PvGTSeq_ampseq, PvGTSeq_ampseq@markers$amplicon %in% PvGTSeq_markers$amplicon)

PvGTSeq_coverage = get_ReadDepth_coverage(PvGTSeq_ampseq, variable = 'site_of_collection_world_region')

PvGTSeq_coverage$plot_read_depth_heatmap

PvGTSeq_locusAmpRate = locus_amplification_rate(PvGTSeq_ampseq, 
                         threshold = .75,
                         strata = 'site_of_collection_world_region', update_loci = F)


PvGTSeq_locusAmpRate$all_loci_performance_plot +
  labs(title = 'A) PvGTSeq',
       y = '# of Amplicons',
       x = 'Propottion of amplified samples')






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



PvGTSeq_LAC_network_048  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 0.47, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_048$plot_network

PvGTSeq_LAC_network_042 = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                          threshold = 0.43, 
                                          metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                                          sample_id = 'Sample_id',
                                          color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_042$plot_network


PvGTSeq_LAC_network_036  = plot_ggnetwork(pairwise_relatedness = Monoclonal_LAC_pairs_PvGTSeq_01SNPsbiallelic_75_IBD,
                                       threshold = 0.37, 
                                       metadata = monoclonals_PvGTSeq_LAC_loci@metadata,
                                       sample_id = 'Sample_id',
                                       color_by = 'site_of_collection_snl0', vertex.size = 4)

PvGTSeq_LAC_network_036$plot_network
