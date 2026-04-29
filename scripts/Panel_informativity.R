
# Call of functions and libraries----

source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_required_libraries.R')
source('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/amplseq_functions.R')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/hmmloglikelihood.cpp')
sourceCpp('~/Documents/Github/MHap-Analysis/docs/functions_and_libraries/Rcpp_functions.cpp')


markers = read.csv('~/Documents/Github/MHap-Analysis/docs/reference/Pviv_P01/PvGTSeq249_markersTable.csv')
markers$amplicon = gsub('(-|\\.|/|:)', '_', markers$amplicon)

q75_markers = markers$amplicon[markers$q75]
drs_markers = markers$amplicon[grepl('DRS', markers$use)]
problematic_markers = markers$amplicon[!grepl('None', markers$Comments)]

# Upload PvGTSeq and WGS data----

PvGTSeq_ampseq_wgs_masked_filtered75 = read_ampseq('~/Documents/Github/PvGTSeq_paper_draft/Pv_Amplicon_data/AmpSeq_Objects/merged_PvGTSeq_WGS_ampseq_masked_filtered75',
                                                   format = 'tsv',
                                                   sep = '\t')

# FIX ERROR in the pipeline of filtering and merging data. Check later the source of the error
#PvGTSeq_ampseq_wgs_masked_filtered75@gt[grepl('^:', PvGTSeq_ampseq_wgs_masked_filtered75@gt)] = NA

PvGTSeq_ampseq_masked_filtered75 = PvGTSeq_ampseq_wgs_masked_filtered75

# Upload three letters code for countries according to world bank----
world_region_codes = read.csv('~/Documents/Github/MHap-Analysis/docs/GADMTools/world-regions-according-to-the-world-bank.csv')

## Keep only markers with good amplification----
PvGTSeq_ampseq_masked_filtered75_q75 = filter_loci(
  PvGTSeq_ampseq_masked_filtered75,
  PvGTSeq_ampseq_masked_filtered75@markers$amplicon %in% q75_markers)

nrow(PvGTSeq_ampseq_masked_filtered75_q75@metadata)

## Keep markers for geographic subdivision----
PvGTSeq_ampseq_masked_filtered75_q75_geo = filter_loci(
  PvGTSeq_ampseq_masked_filtered75_q75,
  !(PvGTSeq_ampseq_masked_filtered75_q75@markers$amplicon %in% drs_markers))

## Remove replicated samples----
PvGTSeq_ampseq_masked_filtered75_q75_geo = remove_replicates(PvGTSeq_ampseq_masked_filtered75_q75_geo, v = 'Old_Sample_id')

nrow(PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata) # 2335 unique samples (PvGTSeq + WGS)
ncol(PvGTSeq_ampseq_masked_filtered75_q75_geo@gt) # 167 amplicons for geographic attribution

## Create an ampseq object for PvGTSeq data----

PvCRiSP_ampseq = filter_loci(PvGTSeq_ampseq_masked_filtered75_q75_geo, 
                    PvGTSeq_ampseq_masked_filtered75_q75_geo@markers$amplicon %in% c('CG2_related',
              'RIPR',
              'VPS11',
              'PIGM'
              )
            )

## Convert data to loci object----
PvGTSeq_loci_masked_filtered75_q75_geo = ampseq2loci(ampseq_object = PvGTSeq_ampseq_masked_filtered75_q75_geo)

PvCRiSP_loci = ampseq2loci(ampseq_object = PvCRiSP_ampseq)

## Run pairwise distance for PvCRiSP-----
if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvCRiSP_2335.csv')){
  pairwise_hamming_distance_PvCRiSP = NULL
  n= 200
  for(w in 1:n){
    print(w)
    start_time = Sys.time()
    pairwise_hamming_distance_PvCRiSP = rbind(
      pairwise_hamming_distance_PvCRiSP,
      pairwise_euclidean(obj = PvCRiSP_loci,
                         w = w, n = n, parallel = T)
    )
    end_time = Sys.time()
    print(end_time -start_time)
  }
  
  write.csv(pairwise_hamming_distance_PvCRiSP,
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvCRiSP_2335.csv',
            row.names = F, quote = F
  )
}else{
  pairwise_hamming_distance_PvCRiSP = 
    read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvCRiSP_2335.csv')
}

### Transform Hamming distance to IBS (1 - Dist) -----

pairwise_hamming_distance_PvCRiSP %<>%
  mutate(rhat = 1 - euDist)


## Run pairwise distance for PvGTSeq -----
if(!file.exists('~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvGTSeq_2335.csv')){
  pairwise_hamming_distance_PvGTSeq = NULL
  n= 1000
  for(w in 501:1000){
    print(w)
    start_time = Sys.time()
    pairwise_hamming_distance_PvGTSeq = rbind(
      pairwise_hamming_distance_PvGTSeq,
      pairwise_euclidean(obj = PvGTSeq_loci_masked_filtered75_q75_geo,
                         w = w, n = n, parallel = T)
    )
    end_time = Sys.time()
    print(end_time -start_time)
  }
  
  write.csv(pairwise_hamming_distance_PvGTSeq,
            '~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvGTSeq_2335.csv',
            row.names = F, quote = F
  )
}else{
  pairwise_hamming_distance_PvGTSeq = 
    read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/pairwise_hamming_distance_PvGTSeq_2335.csv')
}

### Transform Hamming distance to IBS (1 - Dist) -----

pairwise_hamming_distance_PvGTSeq %<>%
  mutate(rhat = 1 - euDist)


## Define Polyclonal infections----

PvGTSeq_ampseq_masked_filtered75_q75_geo = get_polygenomic(PvGTSeq_ampseq_masked_filtered75_q75_geo, strata = 'site_of_collection_snl0', update_popsummary = T)

Fws_q33 = 0.975
Frac_HetLoci_q66 = 0.05


Monoclonals = PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata %>%
  filter(Fws > Fws_q33, Frac_HetLoci < Frac_HetLoci_q66) %>% select(Sample_id) %>% unlist

Polyclonals = PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata %>%
  filter(!(Fws > Fws_q33 & Frac_HetLoci < Frac_HetLoci_q66)) %>% select(Sample_id) %>% unlist


PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata %<>%
  mutate(Clonality = case_when(
    Sample_id %in% Monoclonals ~ 'Monoclonal',
    Sample_id %in% Polyclonals ~ 'Polyclonal'
  ))


PvGTSeq_ampseq_masked_filtered75_q75_geo_mon = 
  filter_samples(PvGTSeq_ampseq_masked_filtered75_q75_geo,
                 PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata$Sample_id %in% Monoclonals)


# Calculate the proportion of samples that has IBS less than 1 using PvGTSeq----

IBS1_PvGTSeq = plot_frac_highly_related(
  pairwise_relatedness = pairwise_hamming_distance_PvGTSeq,
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
  )


Dist0_PvGTSeq = IBS1_PvGTSeq$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])


Populations_with_100pwcomparisons = 
  Dist0_PvGTSeq %>% filter(n >= 45) %>%
  select(Pop_comparison) %>% unlist

# Distribution of background IBS----

IBS_distr_PvGTSeq = plot_relatedness_distribution(
  pairwise_relatedness = pairwise_hamming_distance_PvGTSeq, 
  ncol = 1,
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
  Population = 'site_of_collection_snl0',
  type_pop_comparison = 'within'
)

# IBS_distr_PvCRiSP = plot_relatedness_distribution(
#   pairwise_relatedness = pairwise_hamming_distance_PvCRiSP, 
#   ncol = 1,
#   metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
#   Population = 'site_of_collection_snl0',
#   type_pop_comparison = 'within'
# )


data_plot2 = left_join(IBS_distr_PvGTSeq$relatedness,
                       world_region_codes,
                       by = join_by('Pop_comparison' == 'Entity'))%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% Populations_with_100pwcomparisons)


summ_data_plot2_1 = data_plot2 %>%
  summarise(median_ibs = median(rhat),
            mean_ibs = mean(rhat), .by = Pop_comparison) %>%
  arrange(mean_ibs)

summ_data_plot2_2 = data_plot2 %>%
  summarise(median_ibs = median(rhat),
            mean_ibs = mean(rhat), .by = Pop_comparison) %>%
  arrange(desc(mean_ibs))


summ_data_plot2_2 = left_join(summ_data_plot2_2, 
                              world_region_codes,
                              by = join_by('Pop_comparison' == 'Entity'))


# Calculate the proportion of samples that has IBS less than 1 using PvCRiSP----

IBS1_PvCRiSP = plot_frac_highly_related(
  pairwise_relatedness = pairwise_hamming_distance_PvCRiSP,
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
)


Dist0_PvCRiSP = IBS1_PvCRiSP$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])



Plot_1 = rbind(data.frame(Dist0_PvGTSeq, Set = 'PvGTSeq'),
               data.frame(Dist0_PvCRiSP, Set = 'PvCRiSP')) %>%
  filter(Pop_comparison %in% summ_data_plot2_1$Pop_comparison)%>%
  mutate(Pop_comparison = factor(Pop_comparison, levels = summ_data_plot2_1$Pop_comparison))%>%
  ggplot(aes(y = Pop_comparison, x = prop, shape = Set, color = Set))+
  geom_point(size = 4, position = position_dodge(width = .5))+
  geom_errorbar(aes(xmin = lower, xmax = upper, width = .2), position = position_dodge(width = .5))+
  scale_color_manual(values = c('firebrick3', 'dodgerblue3'))+
  theme_minimal()+
  labs(x = 'Probability two random samples are distinct',
       color = 'Panel',
       shape = 'Panel'
  )+
  theme(legend.position = 'inside',
        legend.position.inside = c(.2, .2),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
  )




Plot_2 = left_join(IBS_distr_PvGTSeq$relatedness,
          world_region_codes,
          by = join_by('Pop_comparison' == 'Entity')
          )%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% unique(Plot_1$data$Pop_comparison)
         )%>%
  mutate(Pop_comparison = factor(Code, levels = summ_data_plot2_2$Code))%>%
  ggplot(aes(x = rhat, fill = Code)) +
  geom_density(position = "stack", alpha = .7)+
  theme_bw()+
  facet_wrap(.~factor(Code,
                      levels = summ_data_plot2_2$Code
                      ), ncol = 1,
             scales = "free_y", switch = 'y')+
  labs(x = "Identity by state")+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
        )




Figure_6 = ggdraw() +
  draw_plot(Plot_1,
            x = 0, width = .7,
            y = 0, height = 1
            ) +
  draw_plot(Plot_2,
            x = 0.7, width = .3,
            y = 0, height = 1)




countries_diffIBS_greater95_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop >= .95)
countries_diffIBS_less95_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop < .95)
#countries_diffIBS_less90_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop < .9)

Honduras_diffIBS_less95_PvGTSeq_prop = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Honduras') %>% select(prop) %>% unlist %>% round(3)

Honduras_diffIBS_less95_PvGTSeq_lower = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Honduras') %>% select(lower) %>% unlist %>% round(3)

Honduras_diffIBS_less95_PvGTSeq_upper = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Honduras') %>% select(upper) %>% unlist %>% round(3)


Peru_diffIBS_less95_PvGTSeq_prop = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Peru') %>% select(prop) %>% unlist %>% round(3)

Peru_diffIBS_less95_PvGTSeq_lower = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Peru') %>% select(lower) %>% unlist %>% round(3)

Peru_diffIBS_less95_PvGTSeq_upper = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Peru') %>% select(upper) %>% unlist %>% round(3)


Malaysia_diffIBS_less95_PvGTSeq_prop = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Malaysia') %>% select(prop) %>% unlist %>% round(3)

Malaysia_diffIBS_less95_PvGTSeq_lower = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Malaysia') %>% select(lower) %>% unlist %>% round(3)

Malaysia_diffIBS_less95_PvGTSeq_upper = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Malaysia') %>% select(upper) %>% unlist %>% round(3)


Panama_diffIBS_less95_PvGTSeq_prop = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(prop) %>% unlist %>% round(3)

Panama_diffIBS_less95_PvGTSeq_lower = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(lower) %>% unlist %>% round(3)

Panama_diffIBS_less95_PvGTSeq_upper = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(upper) %>% unlist %>% round(3)

countries_diffIBS_greater95_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop >= .95)
countries_diffIBS_90_95_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop >= .9 & prop < .95)
countries_diffIBS_less90_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop < .9) %>% arrange(desc(prop))
#countries_diffIBS_less80_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop < .8)


countries_diffIBS_less90_PvCRiSP_text = 
  paste0(
    paste(countries_diffIBS_less90_PvCRiSP$Pop_comparison[
      -length(countries_diffIBS_less90_PvCRiSP$Pop_comparison)], 
      collapse = ', '),
    ' and ',
    paste(countries_diffIBS_less90_PvCRiSP$Pop_comparison[
      length(countries_diffIBS_less90_PvCRiSP$Pop_comparison)]
    ))

Informativity = Plot_1$data

Informativity = left_join(Informativity, summ_data_plot2_1,
                          by = "Pop_comparison")


Informativity %<>% mutate(Clonal_propagation = case_when(
  Pop_comparison %in% c(#'Ethiopia',
                        'Brazil', 
                        'Colombia',
                        'Panama',
                        #'Mexico',
                        'Peru',
                        'Honduras') ~ TRUE,
  .default = FALSE
))


Informativity %>%
  ggplot(aes(x = mean_ibs, y = log((prop - 0.0001)/(1 - (prop - 0.0001))), color = Set, shape = Clonal_propagation)) +
  geom_point()+
  facet_grid(Set~., scales = 'free')

model1 = lm(formula = log((prop - 0.0001)/(1 - (prop - 0.0001))) ~ mean_ibs + Clonal_propagation + Set, data = Informativity)

model1_summ = summary(model1)



model1_summ$adj.r.squared


# non-clonal and clonal comparisons----

# File Monoclonal_pairs_by_country_dist.csv was generated using the script non_clonal_pairwise_comparisons.R
monoclonal_pairwise_comparisons = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/Monoclonal_pairs_by_country_dist.csv')
#clonal_pairwise_comparisons = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/clonal_pairwise_comparisons.csv')
#non_clonal_pairwise_comparisons = read.csv('~/Documents/Github/PvGTSeq_paper_draft/Outputs/non_clonal_pairwise_comparisons2.csv')

## Define clonal groups and comparisons not resolved by WGS----
monoclonal_pairwise_comparisons %<>%
  mutate(Clonal_group = case_when(
    rhat >= 0.999 ~ 'Identical Genomes',
    rhat < 0.999 & rhat >= 0.99 ~ 'Clonal group',
    rhat < 0.99 ~ 'No clonal group'
  ))

## Add metadata----
All_broad_Pv4_PvGTSeq_metadata = read.table('~/Documents/Github/DataManagment_NeafseyLab/Metadata/All_Broad_PvSamples/All_broad_Pv4_PvWGS_metadata.tsv', sep = '\t', header = T)

names(monoclonal_pairwise_comparisons) = paste0(names(monoclonal_pairwise_comparisons), '_WGS')

## Edit Sample_ids of Yi and Yj samples----
monoclonal_pairwise_comparisons = 
  left_join(monoclonal_pairwise_comparisons,
            All_broad_Pv4_PvGTSeq_metadata[, c('Sample_id', 'PvGTSeq_id')],
            by = join_by('Yi_WGS' == 'Sample_id')
  )

monoclonal_pairwise_comparisons$Yi_WGS = monoclonal_pairwise_comparisons$PvGTSeq_id
monoclonal_pairwise_comparisons$PvGTSeq_id = NULL

monoclonal_pairwise_comparisons = 
  left_join(monoclonal_pairwise_comparisons,
            All_broad_Pv4_PvGTSeq_metadata[, c('Sample_id', 'PvGTSeq_id')],
            by = join_by('Yj_WGS' == 'Sample_id')
  )

monoclonal_pairwise_comparisons$Yj_WGS = monoclonal_pairwise_comparisons$PvGTSeq_id
monoclonal_pairwise_comparisons$PvGTSeq_id = NULL


monoclonal_pairwise_comparisons %>% filter(site_of_collection_snl0_WGS == 'Honduras')

## Add PvGTSeq and PvCRiSP data----

### Add to PvGTSeq data----
pairwise_hamming_distance_PvGTSeq2 = pairwise_hamming_distance_PvGTSeq

names(pairwise_hamming_distance_PvGTSeq2) = paste0(names(pairwise_hamming_distance_PvGTSeq2), '_PvGTSeq')


unique(PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata$Seq_Source)

PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata$Sample_id[grepl('Nova', PvGTSeq_ampseq_masked_filtered75_q75_geo@metadata$Sample_id)]


pairwise_hamming_distance_PvGTSeq2$Yi_PvGTSeq = 
  gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', pairwise_hamming_distance_PvGTSeq2$Yi_PvGTSeq)

pairwise_hamming_distance_PvGTSeq2$Yj_PvGTSeq = 
  gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', pairwise_hamming_distance_PvGTSeq2$Yj_PvGTSeq)


pairwise_hamming_distance_PvGTSeq3 = 
  left_join(pairwise_hamming_distance_PvGTSeq2,
            monoclonal_pairwise_comparisons,
            by = join_by('Yi_PvGTSeq' == 'Yi_WGS', 'Yj_PvGTSeq' == 'Yj_WGS')
  )

### Add to PvCRiSP data----
pairwise_hamming_distance_PvCRiSP2 = pairwise_hamming_distance_PvCRiSP

names(pairwise_hamming_distance_PvCRiSP2) = paste0(names(pairwise_hamming_distance_PvCRiSP2), '_PvCRiSP')

pairwise_hamming_distance_PvCRiSP2$Yi_PvCRiSP = 
  gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', pairwise_hamming_distance_PvCRiSP2$Yi_PvCRiSP)

pairwise_hamming_distance_PvCRiSP2$Yj_PvCRiSP = 
  gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', pairwise_hamming_distance_PvCRiSP2$Yj_PvCRiSP)


pairwise_hamming_distance_PvCRiSP3 = 
  left_join(pairwise_hamming_distance_PvCRiSP2,
            monoclonal_pairwise_comparisons,
            by = join_by('Yi_PvCRiSP' == 'Yi_WGS', 'Yj_PvCRiSP' == 'Yj_WGS')
  )


### Add PvGTSeq----

monoclonal_pairwise_comparisons = 
  left_join(monoclonal_pairwise_comparisons,
            pairwise_hamming_distance_PvGTSeq2,
            by = join_by('Yi_WGS' == 'Yi_PvGTSeq', 'Yj_WGS' == 'Yj_PvGTSeq')
  )

### Add PvCRiSP----
monoclonal_pairwise_comparisons = 
  left_join(monoclonal_pairwise_comparisons,
            pairwise_hamming_distance_PvCRiSP2,
            by = join_by('Yi_WGS' == 'Yi_PvCRiSP', 'Yj_WGS' == 'Yj_PvCRiSP')
  )



unique(monoclonal_pairwise_comparisons$Clonal_group_WGS)


monoclonal_pairwise_comparisons %>%
  summarise(min_dist = min(nDiff_WGS),
            max_dist = max(nDiff_WGS), .by = Clonal_group_WGS)


## Non-clonal comparisons----

monoclonal_pairwise_comparisons %>% 
  mutate(Yi = Yi_WGS,
         Yj = Yj_WGS,
         rhat = rhat_PvGTSeq) %>%
  filter(!is.na(rhat), site_of_collection_snl0_WGS == 'Peru') %>% View()

IBS1_PvGTSeq_nonclonal = plot_frac_highly_related(
  pairwise_relatedness = monoclonal_pairwise_comparisons %>% 
    mutate(Yi = Yi_WGS,
           Yj = Yj_WGS,
           rhat = rhat_PvGTSeq) %>%
    filter(!is.na(rhat) & Clonal_group_WGS %in% c('No clonal group')),
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata %>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
)


Dist0_PvGTSeq_nonclonal = IBS1_PvGTSeq_nonclonal$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])


Populations_with_100pwcomparisons_nonclonal = 
  Dist0_PvGTSeq_nonclonal %>% filter(n >= 10) %>%
  select(Pop_comparison) %>% unlist



# Distribution of background IBS----

IBS_distr_PvGTSeq_nonclonal = plot_relatedness_distribution(
  pairwise_relatedness = monoclonal_pairwise_comparisons %>% 
    mutate(Yi = Yi_WGS,
           Yj = Yj_WGS,
           rhat = rhat_WGS) %>%
    #rhat = rhat_WGS) %>%
    filter(!is.na(rhat)), 
  ncol = 1,
  metadata =  PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata %>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),
  Population = 'site_of_collection_snl0',
  type_pop_comparison = 'within'
)

# IBS_distr_PvCRiSP = plot_relatedness_distribution(
#   pairwise_relatedness = pairwise_hamming_distance_PvCRiSP, 
#   ncol = 1,
#   metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
#   Population = 'site_of_collection_snl0',
#   type_pop_comparison = 'within'
# )


data_plot2_nonclonal = left_join(IBS_distr_PvGTSeq_nonclonal$relatedness,
                                 world_region_codes,
                                 by = join_by('Pop_comparison' == 'Entity'))%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% Populations_with_100pwcomparisons_nonclonal)


summ_data_plot2_1_nonclonal = data_plot2_nonclonal %>%
  summarise(median_ibs = median(rhat),
            mean_ibs = mean(rhat), .by = Pop_comparison) %>%
  arrange(mean_ibs)

summ_data_plot2_2_nonclonal = data_plot2_nonclonal %>%
  summarise(median_ibs = median(rhat),
            mean_ibs = mean(rhat), .by = Pop_comparison) %>%
  arrange(desc(mean_ibs))


summ_data_plot2_2_nonclonal = left_join(summ_data_plot2_2_nonclonal, 
                                        world_region_codes,
                                        by = join_by('Pop_comparison' == 'Entity'))


# Calculate the proportion of samples that has IBS less than 1 using PvCRiSP----

IBS1_PvCRiSP_nonclonal = plot_frac_highly_related(
  pairwise_relatedness = monoclonal_pairwise_comparisons %>% 
    mutate(Yi = Yi_WGS,
           Yj = Yj_WGS,
           rhat = rhat_PvCRiSP) %>%
    filter(!is.na(rhat) & Clonal_group_WGS %in% c('No clonal group')),
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata%>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
)


Dist0_PvCRiSP_nonclonal = IBS1_PvCRiSP_nonclonal$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])



Plot_1_nonclonal = rbind(data.frame(Dist0_PvGTSeq_nonclonal, Set = 'PvGTSeq'),
                         data.frame(Dist0_PvCRiSP_nonclonal, Set = 'PvCRiSP')) %>%
  filter(Pop_comparison %in% summ_data_plot2_1_nonclonal$Pop_comparison)%>%
  mutate(Pop_comparison = factor(Pop_comparison, levels = summ_data_plot2_1_nonclonal$Pop_comparison))%>%
  ggplot(aes(y = Pop_comparison, x = prop, shape = Set, color = Set))+
  geom_point(size = 4, position = position_dodge(width = .5))+
  geom_errorbar(aes(xmin = lower, xmax = upper, width = .2), position = position_dodge(width = .5))+
  scale_color_manual(values = c('firebrick3', 'dodgerblue3'))+
  theme_minimal()+
  labs(title = 'A) Discrimination power (D)',
       x = 'Probability that two random samples have different genotypes\ngiven that they do not belong to the same clonal group',
       color = 'Panel',
       shape = 'Panel'
  )+
  theme(legend.position = 'inside',
        legend.position.inside = c(.2, .2),
        title = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
  )


Dist0_PvGTSeq_nonclonal
Dist0_PvCRiSP_nonclonal

Plot_2_nonclonal = left_join(IBS_distr_PvGTSeq_nonclonal$relatedness,
                             world_region_codes,
                             by = join_by('Pop_comparison' == 'Entity')
)%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% unique(Plot_1_nonclonal$data$Pop_comparison)
  )%>%
  mutate(Pop_comparison = factor(Code, levels = summ_data_plot2_2_nonclonal$Code))%>%
  ggplot(aes(x = rhat, fill = Code)) +
  geom_density(position = "stack", alpha = .7)+
  geom_vline(xintercept = .99, linetype = 2)+
  theme_bw()+
  facet_wrap(.~factor(Code,
                      levels = summ_data_plot2_2_nonclonal$Code
  ), ncol = 1,
  scales = "free_y", switch = 'y')+
  labs(title = 'B) IBS distribution',
       x = "Identity by state")+
  theme(title = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )




Figure_6_nonclonal = ggdraw() +
  draw_plot(Plot_1_nonclonal,
            x = 0, width = .7,
            y = 0, height = 1
  ) +
  draw_plot(Plot_2_nonclonal,
            x = 0.7, width = .3,
            y = 0, height = 1)


ggsave('~/Documents/Github/PvGTSeq_paper_draft/Files_to_upload/Fig_4.pdf', 
       Figure_6_nonclonal,
       device = 'pdf',
       #units = 'cm',
       width = 7,
       height = 6.5
)

## Clonal comparisons----


# Calculate the proportion of samples that has IBS less than 1 using PvGTSeq----

pairwise_hamming_distance_PvGTSeq3 %<>%
  mutate(Clonal_group_WGS = 
           case_when(
             is.na(Clonal_group_WGS) ~ 'PvGTSeq',
             .default = Clonal_group_WGS
             ))

pairwise_hamming_distance_PvGTSeq3 %>%
  ggplot(aes(x = rhat_PvGTSeq, y = rhat_WGS))+
  geom_point()

unique(pairwise_hamming_distance_PvGTSeq3$Clonal_group_WGS)



IBS1_PvGTSeq = plot_frac_highly_related(
  pairwise_relatedness = pairwise_hamming_distance_PvGTSeq3 %>%
    mutate(Yi = Yi_PvGTSeq,
           Yj = Yj_PvGTSeq,
           rhat = rhat_PvGTSeq) %>%
    filter(Clonal_group_WGS != "Identical Genomes"),
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata %>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),#HERE
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
)


Dist0_PvGTSeq = IBS1_PvGTSeq$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])


Dist0_PvGTSeq %>% View()

Populations_with_100pwcomparisons = 
  Dist0_PvGTSeq %>% filter(n >= 45) %>%
  select(Pop_comparison) %>% unlist

# Distribution of background IBS----

IBS_distr_PvGTSeq = plot_relatedness_distribution(
  pairwise_relatedness = pairwise_hamming_distance_PvGTSeq3 %>%
    mutate(Yi = Yi_PvGTSeq,
           Yj = Yj_PvGTSeq,
           rhat = rhat_WGS), 
  ncol = 1,
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata %>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),
  Population = 'site_of_collection_snl0',
  type_pop_comparison = 'within'
)

# IBS_distr_PvCRiSP = plot_relatedness_distribution(
#   pairwise_relatedness = pairwise_hamming_distance_PvCRiSP, 
#   ncol = 1,
#   metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata,
#   Population = 'site_of_collection_snl0',
#   type_pop_comparison = 'within'
# )


data_plot2 = left_join(IBS_distr_PvGTSeq$relatedness,
                       world_region_codes,
                       by = join_by('Pop_comparison' == 'Entity'))%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% Populations_with_100pwcomparisons)


data_plot2%>% filter(Pop_comparison == 'Panama') %>% View()

summ_data_plot2_1 = data_plot2 %>%
  summarise(median_ibs = median(rhat, na.rm = T),
            mean_ibs = mean(rhat, na.rm = T), .by = Pop_comparison) %>%
  arrange(mean_ibs)

summ_data_plot2_2 = data_plot2 %>%
  summarise(median_ibs = median(rhat, na.rm = T),
            mean_ibs = mean(rhat, na.rm = T), .by = Pop_comparison) %>%
  arrange(desc(mean_ibs))


summ_data_plot2_2 = left_join(summ_data_plot2_2, 
                              world_region_codes,
                              by = join_by('Pop_comparison' == 'Entity'))


# Calculate the proportion of samples that has IBS less than 1 using PvCRiSP----

pairwise_hamming_distance_PvCRiSP3 %<>%
  mutate(Clonal_group_WGS = 
           case_when(
             is.na(Clonal_group_WGS) ~ 'PvGTSeq',
             .default = Clonal_group_WGS
           ))

IBS1_PvCRiSP = plot_frac_highly_related(
  pairwise_relatedness = pairwise_hamming_distance_PvCRiSP3%>%
    mutate(Yi = Yi_PvCRiSP,
           Yj = Yj_PvCRiSP,
           rhat = rhat_PvCRiSP) %>%
    filter(Clonal_group_WGS != "Identical Genomes"),
  metadata = PvGTSeq_ampseq_masked_filtered75_q75_geo_mon@metadata %>%
    mutate(Sample_id = gsub('(_iSeq|_MiSeq|_NovaSeq1|_NovaSeq2|_broadWGS)', '', Sample_id)),
  Population = 'site_of_collection_snl0',
  threshold = 1,
  type_pop_comparison = 'within'
)


Dist0_PvCRiSP = IBS1_PvCRiSP$highly_related_table %>% 
  group_by(Pop_comparison) %>%
  mutate(prop = binconf(n - freq,
                        n,
                        alpha = 0.05,
                        method = "exact")[1],
         lower = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[2],
         upper = binconf(n - freq,
                         n,
                         alpha = 0.05,
                         method = "exact")[3])



Plot_1 = rbind(data.frame(Dist0_PvGTSeq, Set = 'PvGTSeq'),
               data.frame(Dist0_PvCRiSP, Set = 'PvCRiSP')) %>%
  filter(Pop_comparison %in% summ_data_plot2_1$Pop_comparison)%>%
  mutate(Pop_comparison = factor(Pop_comparison, levels = summ_data_plot2_1$Pop_comparison))%>%
  ggplot(aes(y = Pop_comparison, x = prop, shape = Set, color = Set))+
  geom_point(size = 4, position = position_dodge(width = .5))+
  geom_errorbar(aes(xmin = lower, xmax = upper, width = .2), position = position_dodge(width = .5))+
  scale_color_manual(values = c('firebrick3', 'dodgerblue3'))+
  theme_minimal()+
  labs(title = 'A) Discrimination power (D)',
       x = 'Probability that two random samples have different genotypes\ngiven that they are different at genome level',
       color = 'Panel',
       shape = 'Panel'
  )+
  theme(legend.position = 'inside',
        legend.position.inside = c(.2, .2),
        title = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 9),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10)
  )




Plot_2 = left_join(IBS_distr_PvGTSeq$relatedness,
                   world_region_codes,
                   by = join_by('Pop_comparison' == 'Entity')
)%>%
  filter(Type_of_comparison == 'Within', 
         Pop_comparison %in% unique(Plot_1$data$Pop_comparison)
  )%>%
  mutate(Pop_comparison = factor(Code, levels = summ_data_plot2_2$Code))%>%
  ggplot(aes(x = rhat, fill = Code)) +
  geom_density(position = "stack", alpha = .7)+
  geom_vline(xintercept = .999, linetype = 2)+
  theme_bw()+
  facet_wrap(.~factor(Code,
                      levels = summ_data_plot2_2$Code
  ), ncol = 1,
  scales = "free_y", switch = 'y')+
  labs(title = 'B) IBS distribution', x = "Identity by state")+
  theme(title = element_text(size  = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )




Figure_6 = ggdraw() +
  draw_plot(Plot_1,
            x = 0, width = .7,
            y = 0, height = 1
  ) +
  draw_plot(Plot_2,
            x = 0.7, width = .3,
            y = 0, height = 1)

ggsave('~/Documents/Github/PvGTSeq_paper_draft/Files_to_upload/Fig_S10.pdf', 
       Figure_6,
       device = 'pdf',
       #units = 'cm',
       width = 7,
       height = 6.5
)



countries_diffIBS_greater95_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop >= .95)
countries_diffIBS_less95_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop < .95)
#countries_diffIBS_less90_PvGTSeq = Plot_1$data %>% filter(Set == 'PvGTSeq', prop < .9)

# Honduras_diffIBS_less95_PvGTSeq_prop = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Honduras') %>% select(prop) %>% unlist %>% round(3)
# 
# Honduras_diffIBS_less95_PvGTSeq_lower = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Honduras') %>% select(lower) %>% unlist %>% round(3)
# 
# Honduras_diffIBS_less95_PvGTSeq_upper = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Honduras') %>% select(upper) %>% unlist %>% round(3)
# 
# 
# Peru_diffIBS_less95_PvGTSeq_prop = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Peru') %>% select(prop) %>% unlist %>% round(3)
# 
# Peru_diffIBS_less95_PvGTSeq_lower = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Peru') %>% select(lower) %>% unlist %>% round(3)
# 
# Peru_diffIBS_less95_PvGTSeq_upper = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Peru') %>% select(upper) %>% unlist %>% round(3)
# 
# 
# Malaysia_diffIBS_less95_PvGTSeq_prop = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Malaysia') %>% select(prop) %>% unlist %>% round(3)
# 
# Malaysia_diffIBS_less95_PvGTSeq_lower = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Malaysia') %>% select(lower) %>% unlist %>% round(3)
# 
# Malaysia_diffIBS_less95_PvGTSeq_upper = 
#   countries_diffIBS_less95_PvGTSeq %>% 
#   filter(Pop_comparison == 'Malaysia') %>% select(upper) %>% unlist %>% round(3)
# 

Panama_diffIBS_less95_PvGTSeq_prop = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(prop) %>% unlist %>% round(3)

Panama_diffIBS_less95_PvGTSeq_lower = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(lower) %>% unlist %>% round(3)

Panama_diffIBS_less95_PvGTSeq_upper = 
  countries_diffIBS_less95_PvGTSeq %>% 
  filter(Pop_comparison == 'Panama') %>% select(upper) %>% unlist %>% round(3)

countries_diffIBS_greater95_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop >= .95)
countries_diffIBS_90_95_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop >= .9 & prop < .95)
countries_diffIBS_less90_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop < .9) %>% arrange(desc(prop))
#countries_diffIBS_less80_PvCRiSP = Plot_1$data %>% filter(Set == 'PvCRiSP', prop < .8)


countries_diffIBS_less90_PvCRiSP_text = 
  paste0(
    paste(countries_diffIBS_less90_PvCRiSP$Pop_comparison[
      -length(countries_diffIBS_less90_PvCRiSP$Pop_comparison)], 
      collapse = ', '),
    ' and ',
    paste(countries_diffIBS_less90_PvCRiSP$Pop_comparison[
      length(countries_diffIBS_less90_PvCRiSP$Pop_comparison)]
    ))

Informativity = Plot_1$data

Informativity = left_join(Informativity, summ_data_plot2_1,
                          by = "Pop_comparison")


Informativity %<>% mutate(Clonal_propagation = case_when(
  Pop_comparison %in% c(#'Ethiopia',
    'Brazil', 
    'Colombia',
    'Panama',
    #'Mexico',
    'Peru',
    'Honduras') ~ TRUE,
  .default = FALSE
))


Informativity %>%
  ggplot(aes(x = mean_ibs, y = log((prop - 0.0001)/(1 - (prop - 0.0001))), color = Set, shape = Clonal_propagation)) +
  geom_point()+
  facet_grid(Set~., scales = 'free')

model1 = lm(formula = log((prop - 0.0001)/(1 - (prop - 0.0001))) ~ mean_ibs + Clonal_propagation + Set, data = Informativity)

model1_summ = summary(model1)



model1_summ$adj.r.squared
