
wd = 'C:/Users/jerry/Documentos/GitHub/amplicon-panel-benchmark/src/'
fd = 'C:/Users/jerry/Documentos/GitHub/amplicon-panel-benchmark/src'
rd = 'C:/Users/jerry/Documentos/GitHub/amplicon-panel-benchmark/docs/reference/Pviv_P01'

# Get functions and libraries----
source(file.path(fd, 'amplseq_required_libraries.R'))
source(file.path(fd, 'functions/amplseq_functions.R'))
source(file.path(fd, 'functions/rGenome_functions.R'))

sourceCpp(file.path(fd, 'functions/Rcpp_functions.cpp'))
sourceCpp(file.path(fd, 'functions/hmmloglikelihood.cpp'))

install.packages('Rtools')

# PvGTSeq ----

PvGTSeq_markers = read.csv(files.path(rd,'PvGTSeq249_markersTable.csv'))

PvGTSeq_markers %<>% filter(use !='DRS')

## Extract PvGTSeq coordinates from rGenome object----

PvGTSeq_coordinates = NULL
for(amplicon in 1:nrow(PvGTSeq_markers)){
  amplicon_chromosome = gsub('v1$', 'v2', PvGTSeq_markers[amplicon,][['chromosome']])
  amplicon_start = PvGTSeq_markers[amplicon,][['start']]
  amplicon_end = PvGTSeq_markers[amplicon,][['end']]
  
  PvGTSeq_coordinates = c(PvGTSeq_coordinates,
                          paste(amplicon_chromosome,
                                amplicon_start:amplicon_end, sep = '_'))
  
}

## load vcf and generate rGenome object----
PvGTSeq_vcf = load_vcf('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/Pv_all_samples_PvGTSeq_SelectedFromPositions.recode.vcf')

PvGTSeq_rGenome = vcf2rGenome(vcf = PvGTSeq_vcf, n = 100, threshold = 5)

## Keep only positions for geographic differentiation---

dim(PvGTSeq_rGenome@gt)

PvGTSeq_rGenome = filter_loci(PvGTSeq_rGenome, v = rownames(PvGTSeq_rGenome@gt) %in% PvGTSeq_coordinates)

dim(PvGTSeq_rGenome@gt)

## Calculate the fraction of amplified loci per sample----

PvGTSeq_rGenome@metadata$PvGTSeq_SampAmpRate = SampleAmplRate(PvGTSeq_rGenome, update = F)

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



## Extract PvAmpliSeq coordinates from rGenome object----

PvAmpliSeq_coordinates = NULL
for(amplicon in 1:nrow(PvAmpliSeq_markers)){
  amplicon_chromosome = gsub('v1$', 'v2', PvAmpliSeq_markers[amplicon,][['chromosome']])
  amplicon_start = PvAmpliSeq_markers[amplicon,][['start']]
  amplicon_end = PvAmpliSeq_markers[amplicon,][['end']]
  
  PvAmpliSeq_coordinates = c(PvAmpliSeq_coordinates,
                             paste(amplicon_chromosome,
                                   amplicon_start:amplicon_end, sep = '_'))
  
}


## load vcf and generate rGenome object----

PvAmpliSeq_vcf = load_vcf('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/Pv_all_samples_AmpliSeq_SelectedFromPositions.recode.vcf')

PvAmpliSeq_rGenome = vcf2rGenome(vcf = PvAmpliSeq_vcf, n = 100, threshold = 5)

## Keep only positions for geographic differentiation---

dim(PvAmpliSeq_rGenome@gt)

PvAmpliSeq_rGenome = filter_loci(PvAmpliSeq_rGenome, v = rownames(PvAmpliSeq_rGenome@gt) %in% PvAmpliSeq_coordinates)

dim(PvAmpliSeq_rGenome@gt)

## Calculate the fraction of amplified loci per sample----

PvAmpliSeq_rGenome@metadata$PvAmpliSeq_SampAmpRate = SampleAmplRate(PvAmpliSeq_rGenome, update = F)


# rhAmpSeq ----

rhAmpSeq_markers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/rhAmpSeq_v2_markers_table.csv')

rhAmpSeq_markers %<>% filter(!grepl('(MIT|MDR1|DHPS)', amplicon_name))

length(rhAmpSeq_markers$amplicon_name)

## Extract rhAmpSeq coordinates from rGenome object----

rhAmpSeq_coordinates = NULL
for(amplicon in 1:nrow(rhAmpSeq_markers)){
  amplicon_chromosome = rhAmpSeq_markers[amplicon,][['chromosome']]
  amplicon_start = rhAmpSeq_markers[amplicon,][['start']]
  amplicon_end = rhAmpSeq_markers[amplicon,][['end']]
  
  rhAmpSeq_coordinates = c(rhAmpSeq_coordinates,
                           paste(amplicon_chromosome,
                                 amplicon_start:amplicon_end, sep = '_'))
  
}

## load vcf and generate rGenome object----

rhAmpSeq_vcf = load_vcf('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/Pv_all_samples_rhAmpSeq_SelectedFromPositions.recode.vcf')

rhAmpSeq_rGenome = vcf2rGenome(vcf = rhAmpSeq_vcf, n = 100, threshold = 5)

## Keep only positions for geographic differentiation---

dim(rhAmpSeq_rGenome@gt)

rhAmpSeq_rGenome = filter_loci(rhAmpSeq_rGenome, v = rownames(rhAmpSeq_rGenome@gt) %in% rhAmpSeq_coordinates)

dim(rhAmpSeq_rGenome@gt)

## Calculate the fraction of amplified loci per sample----

rhAmpSeq_rGenome@metadata$rhAmpSeq_SampAmpRate = SampleAmplRate(rhAmpSeq_rGenome, update = F)


# WGS ----

# Load Genome data-----
load('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome.RData')

# Fix missing data----
WGS_all_samples_rGenome@gt[grepl('NA',WGS_all_samples_rGenome@gt)] = NA

## Calculate the fraction of amplified loci per sample----

WGS_all_samples_rGenome@metadata$WGS_SampAmpRate = SampleAmplRate(WGS_all_samples_rGenome, update = F)


## Upload metadata ----

terra_all_samples_metadata = read.table('~/Documents/Github/DataManagment_NeafseyLab/Metadata/Terra_gcloud/Terra_metadata_harmonized.tsv', sep = '\t', header = T)

# Update metadata----

terra_all_samples_metadata = left_join(terra_all_samples_metadata,
                                       PvGTSeq_rGenome@metadata %>% select(Sample_id, PvGTSeq_SampAmpRate),
                                       by = join_by('entity.sample_id' == 'Sample_id')
)

terra_all_samples_metadata = left_join(terra_all_samples_metadata,
                                       PvAmpliSeq_rGenome@metadata %>% select(Sample_id, PvAmpliSeq_SampAmpRate),
                                       by = join_by('entity.sample_id' == 'Sample_id')
)

terra_all_samples_metadata = left_join(terra_all_samples_metadata,
                                       rhAmpSeq_rGenome@metadata %>% select(Sample_id, rhAmpSeq_SampAmpRate),
                                       by = join_by('entity.sample_id' == 'Sample_id')
)

terra_all_samples_metadata = left_join(terra_all_samples_metadata,
                                       WGS_all_samples_rGenome@metadata %>% select(Sample_id, WGS_SampAmpRate),
                                       by = join_by('entity.sample_id' == 'Sample_id')
)


terra_all_samples_metadata %>%
  summarise(
    nSamples = n(),
    nSamples75 = sum(PvGTSeq_SampAmpRate >= .75 & 
                               PvAmpliSeq_SampAmpRate >= .75 &
                               rhAmpSeq_SampAmpRate >= .75 &
                               WGS_SampAmpRate >= .75, na.rm = T),
    nSamples75PvGTSeq = sum(PvGTSeq_SampAmpRate >= .75, na.rm = T),
    nSamples75PvAmpliSeq = sum(PvAmpliSeq_SampAmpRate >= .75, na.rm = T),
    nSamples75rhAmpSeq = sum(rhAmpSeq_SampAmpRate >= .75, na.rm = T),
    nSamples75WGS = sum(WGS_SampAmpRate >= .75, na.rm = T),
            .by = c(site_of_collection_world_region,
                    site_of_collection_snl0)
            ) %>%
  arrange(site_of_collection_world_region, site_of_collection_snl0)

write.table(terra_all_samples_metadata,
      '~/Documents/Github/DataManagment_NeafseyLab/Metadata/Terra_gcloud/Terra_metadata_harmonized.tsv', sep = '\t',
      quote = F,
      row.names = F
      )


## Add metadata----

PvGTSeq_rGenome@metadata = left_join(PvGTSeq_rGenome@metadata,
                                     terra_all_samples_metadata,
                                     by = join_by('Sample_id' == 'entity.sample_id'))


write_rGenome(PvGTSeq_rGenome, '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_all_samples_PvGTSeq_rGenome', 
              format = 'tsv',
              sep = '\t')


## Add metadata----

PvAmpliSeq_rGenome@metadata = left_join(PvAmpliSeq_rGenome@metadata,
                                        terra_all_samples_metadata,
                                        by = join_by('Sample_id' == 'entity.sample_id'))


write_rGenome(PvAmpliSeq_rGenome, '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_all_samples_AmpliSeq_rGenome', 
              format = 'tsv',
              sep = '\t')

## Add metadata----

rhAmpSeq_rGenome@metadata = left_join(rhAmpSeq_rGenome@metadata,
                                      terra_all_samples_metadata,
                                      by = join_by('Sample_id' == 'entity.sample_id'))


write_rGenome(rhAmpSeq_rGenome, '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/all_all_samples_rhAmpSeq_rGenome', 
              format = 'tsv',
              sep = '\t')

## Add metadata----
WGS_all_samples_rGenome@metadata = left_join(WGS_all_samples_rGenome@metadata,
                                             terra_all_samples_metadata,
                                             by = join_by('Sample_id' == 'entity.sample_id'))

write_rGenome(WGS_all_samples_rGenome, 
              format = 'tsv',
              name = '~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/WGS_all_samples_rGenome',
              sep = '\t')
