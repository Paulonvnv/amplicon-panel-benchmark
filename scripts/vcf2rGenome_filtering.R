#!/bin/r env

print('R tool started')

library(argparse)
library(stringr)

parser = ArgumentParser()

parser$add_argument("-wd", "--wd", 
                    help="Path to input to files and output files")

parser$add_argument("-fd", "--fd", default = 'NaN',
                    help="Path to function files and reference files")

parser$add_argument("-vcf", "--vcf_file",  default = 'NaN',
                    help="name of the vcf file")

parser$add_argument("-o", "--out",
                    help="Prefix of output files")

parser$add_argument("-HardFilteringReport", "--HardFilteringReport", default = FALSE,
                    help="Boolean that indicate to generate or not data for the HardFilteringReport")

parser$add_argument("-remove_filtered_all", "--remove_filtered_all", default = FALSE,
                    help="Boolean that indicate to remove site that don't pass the hard filter")

parser$add_argument("-include_from_gff", "--include_from_gff", default = 'NaN',
                    help="regular expression pattern to be used to identify the desired feature to SELECT from the gff file")

parser$add_argument("-exclude_from_gff", "--exclude_from_gff", default = 'NaN',
                    help="regular expression pattern to be used to identify the desired feature to EXCLUDE from the gff file")

parser$add_argument("-gff", "--ref_gff", default = 'NaN',
                    help="name of .gff file containing coordinates of genomic regions")

parser$add_argument("-bed", "--bed", default = 'NaN',
                    help="name of .bed file containing core genomic regions to keep from the VCF file")

parser$add_argument("-ebed", "--exclude_bed", default = 'NaN',
                    help="name of .bed file containing no core genomic regions to exclude from the VCF file")

parser$add_argument("-rkeep", "--keep_regexp", default = 'NaN',
                    help="Regular expression to identify samples of interest")

parser$add_argument("-positions", "--positions", default = 'NaN',
                    help="name of .bed file containing positions to keep from the VCF file")

parser$add_argument("-non_ref_ac_any", "--non_ref_ac_any", default = 'NaN',
                    help = "integer that indicate to remove site sites with any Non-Reference (ALT) Counts (ac) within the range specified")

parser$add_argument("-t", "--nTasks", default = 20,
                    help="Number of Tasks arrays to split the convertion of the vcf to the rGenome format")

parser$add_argument("-do_vcf2rGenome", "--do_vcf2rGenome", default = FALSE,
                    help="Boolean that indicate to run or not vcf2rGenome")

parser$add_argument("-coverage_data", "--coverage_data", default = FALSE,
                    help="Boolean that indicate to generate or not coverage data")


parser$add_argument("-metadata", "--metadata", default = 'NaN',
                    help="string with the file name of the metadata, Sample ID's should be labeled as Sample_id")

parser$add_argument("-join_by", "--join_by", default = 'NaN',
                    help="string that indicates the variable where Sample ID's are stored. Just if  Sample ID's are not labeled as Sample_id")

parser$add_argument("-pop", "--pop", default = 'NaN',
                    help="Name of that variable in the metadata that indicates how to stratify the analysis")


parser$add_argument("-merge_rgenome", "--merge_rgenome", default = FALSE,
                    help="Boolean that indicate to merge or not rGenome chunks")

parser$add_argument("-tid", "--Task_id", default = 'NaN',
                    help="Tasks ID")

parser$add_argument("-n", "--nchunks", default = 100,
                    help="Number of iterations to split the vcf file")

parser$add_argument("-RDthres", "--ReadDepthThreshold", default = 5,
                    help="Minimun read depth to call an allele")

parser$add_argument("-post_filtering", "--post_filtering", default = FALSE,
                    help="Boolean that indicate to perform the postfiltering")

parser$add_argument("-type_of_polymorphism_to_remove", "--type_of_polymorphism_to_remove", default = 'NaN',
                    help="type of polymorphism to be removed")

parser$add_argument("-samprate", "--sample_ampl_rate", default = .75,
                    help="Minimun amplification rate of the sample")

parser$add_argument("-lamprate", "--locus_ampl_rate", default = .75,
                    help="Minimun amplification rate of the locus")

parser$add_argument("-ohetq", "--ObsHet_quantile", default = .95,
                    help="Observed Heterozygosity quantile")

parser$add_argument("-snpdq", "--SNP_density_quantile", default = .95,
                    help="SNP density quantile")

parser$add_argument("-alignment_filter", "--alignment_filter", default = "or",
                    help="Filter to remove alignment errors, options: 'het', 'snp_dens', 'and', 'or'")

parser$add_argument("-mask_formula", "--mask_formula", default = "h_ij>=0.66&h_ijminor>=0.9",
                    help="Fraction of heterozygous samples per alternativ allele threshold")

# Variables----

print("starting to parse variables")

args = parser$parse_args()

# Working directory
wd = args$wd

print(paste0('wd: ', wd))

setwd(wd) 

# Tools or functions directory
fd = args$fd

print(paste0('fd: ', fd))

# Starting vcf file
vcf_file = args$vcf_file

print(paste0('vcf_file: ', vcf_file))

# first filter: HardFilteringReport (optional)
HardFilteringReport = as.logical(args$HardFilteringReport)

print(paste0('HardFilteringReport: ', HardFilteringReport))

# first filter: remove_filtered_all (optional)
remove_filtered_all = as.logical(args$remove_filtered_all)

print(paste0('remove_filtered_all: ', remove_filtered_all))

# Second filter: Coding or non-coding regions (optional)
include_from_gff = args$include_from_gff
include_from_gff = ifelse(include_from_gff == 'NaN', NA, include_from_gff)

print(paste0('include_from_gff: ', include_from_gff))

# Second filter: CodingRegionsOnly (optional)
exclude_from_gff = args$exclude_from_gff
exclude_from_gff = ifelse(exclude_from_gff == 'NaN', NA, exclude_from_gff)

print(paste0('exclude_from_gff: ', exclude_from_gff))

ref_gff_file = args$ref_gff
ref_gff_file = ifelse(ref_gff_file == 'NaN', NA, file.path(fd, ref_gff_file))

print(paste0('ref_gff_file: ', ref_gff_file))

# Third filter: Core Regions Only

bed = args$bed
bed = ifelse(bed == 'NaN', NA, file.path(fd, bed))

print(paste0('bed: ', bed))

exclude_bed = args$exclude_bed
exclude_bed = ifelse(exclude_bed == 'NaN', NA, file.path(fd, exclude_bed))

print(paste0('exclude_bed: ', exclude_bed))


# Fourth filter: Subset of samples
keep_regexp = args$keep_regexp
keep_regexp = ifelse(keep_regexp == 'NaN', NA, keep_regexp)

print(paste0('keep_regexp: ', keep_regexp))

# Fifth filter: positions
positions = args$positions
positions = ifelse(positions == 'NaN', NA, file.path(fd, positions))

# Fifth filter: Only polymorphic sites
non_ref_ac_any = as.integer(args$non_ref_ac_any)
non_ref_ac_any = ifelse(non_ref_ac_any == 'NaN', NA, non_ref_ac_any)

print(paste0('non_ref_ac_any: ', non_ref_ac_any))

# Number of tasks arrays
nTasks = as.integer(args$nTasks)

print(paste0('nTasks: ', nTasks))

# Convert VCF to rGenome
do_vcf2rGenome = as.logical(args$do_vcf2rGenome)

print(paste0('do_vcf2rGenome: ', do_vcf2rGenome))

# Coverage data report
coverage_data = as.logical(args$coverage_data)

print(paste0('coverage_data: ', coverage_data))

# metadata
metadata = as.character(args$metadata)
metadata = ifelse(metadata == 'NaN', NA, metadata)

print(paste0('metadata: ', metadata))

# join metadata by
join_by = as.character(args$join_by)
join_by = ifelse(join_by == 'NaN', NA, join_by)

print(paste0('join_by: ', join_by))

# pop
pop = args$pop
pop = ifelse(pop == 'NaN', NA, pop)

print(paste0('pop: ', pop))

# merge_rgenome
merge_rgenome = as.logical(args$merge_rgenome)

print(paste0('merge_rgenome: ', merge_rgenome))

post_filtering = as.logical(args$post_filtering)

print(paste0('post_filtering: ', post_filtering))

# Task_id
Task_id = as.integer(args$Task_id)

print(paste0('Task_id: ', Task_id))

# Number of iterations to split the vcf file
nchunks = as.integer(args$nchunks)

print(paste0('nchunks: ', nchunks))

# Minimun read depth to call an allele
ReadDepthThreshold = as.numeric(args$ReadDepthThreshold)

print(paste0('ReadDepthThreshold: ', ReadDepthThreshold))

# Type of polymorphism to be removed

type_of_polymorphism_to_remove = args$type_of_polymorphism_to_remove

if(type_of_polymorphism_to_remove == 'NaN'){
  type_of_polymorphism_to_remove =  NA
}else{
  type_of_polymorphism_to_remove = str_split(type_of_polymorphism_to_remove, ';')[[1]]
}

print(paste0('type_of_polymorphism_to_remove: ', type_of_polymorphism_to_remove))

# sample_ampl_rate
sample_ampl_rate = as.numeric(args$sample_ampl_rate)

print(paste0('sample_ampl_rate: ', sample_ampl_rate))

# locus_ampl_rate
locus_ampl_rate  = as.numeric(args$locus_ampl_rate)

print(paste0('locus_ampl_rate: ', locus_ampl_rate))

# ObsHet_quantile
ObsHet_quantile = as.numeric(args$ObsHet_quantile)

print(paste0('ObsHet_quantile: ', ObsHet_quantile))

# SNP_density_quantile
SNP_density_quantile = as.numeric(args$SNP_density_quantile)

print(paste0('SNP_density_quantile: ', SNP_density_quantile))

# Alignment filter
alignment_filter = as.character(args$alignment_filter)

print(paste0('alignment_filter: ', alignment_filter))

# mask_formula filter
mask_formula = as.character(args$mask_formula)
mask_formula = gsub('"',"",mask_formula)

mask_formula = gsub('&'," & ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('\\|'," \\| ",mask_formula, ignore.case = TRUE)

if(grepl("\\w>\\d",mask_formula)){
  patterns = str_extract_all(mask_formula, "\\w>\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('>',' > ',pattern)
    mask_formula = gsub(pattern,
                        replacement,
                        mask_formula, ignore.case = TRUE)
  }
  
}

if(grepl("\\w<\\d",mask_formula)){
  patterns = str_extract_all(mask_formula, "\\w<\\d")[[1]]
  
  for(pattern in patterns){
    
    replacement = gsub('<',' < ',pattern)
    mask_formula = gsub(pattern,
                        replacement,
                        mask_formula, ignore.case = TRUE)
  }
  
}

mask_formula = gsub('>='," >= ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('<='," <= ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('=='," == ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('!='," != ",mask_formula, ignore.case = TRUE)

mask_formula = gsub('\\+'," \\+ ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('-'," - ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('\\*'," \\* ",mask_formula, ignore.case = TRUE)
mask_formula = gsub('/'," / ",mask_formula, ignore.case = TRUE)


print(paste0('mask_formula: ', mask_formula))

# output pattern
output = args$out

print(paste0('output: ', output))

# rGenome_object name
if(nTasks > 1){
  rGenome_object_name = paste0(output, '_rGenome_Chunk', Task_id)

}else{
  rGenome_object_name = paste0(output, '_rGenome')

}

print(paste0('rGenome_object_name: ', rGenome_object_name))

# R image name

if(nTasks > 1){
  imagename = paste0('Chunks/',output, '_rGenome_Chunk', Task_id, '.RData')

}else{
  imagename = paste0(output, '_rGenome.RData')

}

print(paste0('imagename: ', imagename))

final_vcf_name1 = paste0(output, '_filtered_final.vcf')

print(paste0('final_vcf_name1: ', final_vcf_name1))

final_vcf_name2 = paste0(output, '_filtered_masked_final.vcf')

print(paste0('final_vcf_name2: ', final_vcf_name2))

# Check packages and functions----

source(file.path(fd,'load_libraries.R'))
source(file.path(fd,'functions.R'))

print('Functions loaded')
#sourceCpp(file.path(fd,'Rcpp_functions.cpp'))

# Section 1: filter of only polymorphic sites----

if(!is.na(non_ref_ac_any)){
  
  output = paste0(output, '_NonREFacAny', non_ref_ac_any)
  
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_non_ref_ac_any.sh',
               non_ref_ac_any = non_ref_ac_any,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  vcf_file = paste0(output, '.recode.vcf')
  
  print("NonREFacAny done")
  system('echo "NonREFacAny done"')
  
}

# Section 2: Hard Filtering Report----

if(HardFilteringReport){
  if(nTasks > 1){
    
    print('Starting HardFilteringReport')
    ### Load vcf data----  
    # Split the genome in equal peaces
    files_nrows = as.integer(system(paste0("grep -v '^#' ", vcf_file, " | wc -l"), intern = TRUE))
    s = round(seq(1, files_nrows + 1, length.out = nTasks + 1))
    
    print(paste0('there are ', files_nrows, ' positions in the file'))
    
    print(paste0('the windows are ', paste(s, collapse = ', ')))
    
    HardFilteringStats = NULL
    
    for(w in 1:nTasks){
      
      system(paste0('cp ', vcf_file, ' ', gsub('.vcf$', paste0('_', w, '.vcf'), vcf_file)))
      
      temp_vcf_file = gsub('.vcf$', paste0('_', w, '.vcf'), vcf_file)
      
      print(paste0('uploading window ', w, ' of ', nTasks))
      vcf_start = s[w]
      print(paste0('window ', w, ' starts at postion ', vcf_start))
      
      vcf_end = s[w + 1] - 1
      print(paste0('window ', w, ' ends at postion ', vcf_end))
      
      # Upload the VCF to R environment
      temp_vcf_object = load_vcf(vcf = temp_vcf_file, na.rm = TRUE, start = vcf_start, end = vcf_end)
      
      print(paste0('chunk ', w, ' of ', nTasks, ' uploaded'))
      
      system(paste0('rm ', temp_vcf_file))
      
      temp_vcf_object$Type_of_polymorphism = get_type_of_polymorphism(temp_vcf_object)
      
      HardFilteringStats = rbind(HardFilteringStats,
                                data.frame(
                                  temp_vcf_object[,c('CHROM', 'POS', 'REF', 'ALT', 'FILTER','Type_of_polymorphism')],
                                  QUAL = as.numeric(temp_vcf_object$QUAL),
                                  QD = as.numeric(gsub('QD=','',str_extract(temp_vcf_object$INFO, 'QD=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
                                  SOR = as.numeric(gsub('SOR=','',str_extract(temp_vcf_object$INFO, 'SOR=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
                                  FS = as.numeric(gsub('FS=','',str_extract(temp_vcf_object$INFO, 'FS=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
                                  MQ = as.numeric(gsub('MQ=','',str_extract(temp_vcf_object$INFO, 'MQ=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
                                  MQRankSum = as.numeric(gsub('MQRankSum=','',str_extract(temp_vcf_object$INFO, 'MQRankSum=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
                                  ReadPosRankSum = as.numeric(gsub('ReadPosRankSum=','',str_extract(temp_vcf_object$INFO, 'ReadPosRankSum=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?')))
                                ))
      
      rm(list = c('w', 'temp_vcf_object', 'vcf_start', 'vcf_end'))
      
    }
    
    HardFilteringStats %<>% pivot_longer(cols = all_of(c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
                                        values_to = 'Value', names_to = 'Stat') %>%
      mutate(Stat = factor(Stat, levels = c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
             Type_of_polymorphism2 = case_when(
               grepl('SNP',Type_of_polymorphism) ~ 'SNP',
               !grepl('SNP',Type_of_polymorphism) ~ 'INDEL'),
             Thres = case_when(
               grepl('SNP',Type_of_polymorphism) & grepl('QUAL',Stat) ~ 30,
               grepl('SNP',Type_of_polymorphism) & grepl('QD',Stat) ~ 2,
               grepl('SNP',Type_of_polymorphism) & grepl('SOR',Stat) ~ 3,
               grepl('SNP',Type_of_polymorphism) & grepl('FS',Stat) ~ 60,
               grepl('SNP',Type_of_polymorphism) & grepl('MQ',Stat) ~ 40,
               grepl('SNP',Type_of_polymorphism) & grepl('MQRankSum',Stat) ~ -12.5,
               grepl('SNP',Type_of_polymorphism) & grepl('ReadPosRankSum',Stat) ~ 8,
               
               grepl('INDEL',Type_of_polymorphism) & grepl('QUAL',Stat) ~ 30,
               grepl('INDEL',Type_of_polymorphism) & grepl('QD',Stat) ~ 2,
               grepl('INDEL',Type_of_polymorphism) & grepl('SOR',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('FS',Stat) ~ 200,
               grepl('INDEL',Type_of_polymorphism) & grepl('MQ',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('MQRankSum',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('ReadPosRankSum',Stat) ~ 8
             )
      )
    
    plot_HardFilteringStats_density =  HardFilteringStats %>%
      ggplot(aes(x = Value, fill = Type_of_polymorphism2, color = Type_of_polymorphism2)) +
      geom_density(alpha = .4) +
      geom_vline(data = data.frame(Type_of_polymorphism2 = rep(c('SNP', 'INDEL'), each = 7),
                                   Stat = factor(rep(c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum'), 2),
                                                 levels = c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
                                   Thres = c(30, 2, 3, 60, 40, -12.5, -8, 
                                             30, 2, NA, 200, NA, NA, -20)),
                 aes(xintercept = Thres, color = Type_of_polymorphism2)
      ) + 
      facet_wrap(Stat ~., scale = 'free', ncol = 2) +
      theme_bw()+
      theme(legend.position = c(.7,.1))
    
    save(file = 'HardFilteringReport.RData', list = c('HardFilteringStats', 'plot_HardFilteringStats_density'))
    
  }else{
    
    ### Load vcf data----  
    temp_vcf_object = load_vcf(vcf = vcf_file, na.rm = TRUE)
    
    temp_vcf_object$Type_of_polymorphism = get_type_of_polymorphism(temp_vcf_object)
    
    HardFilteringStats = data.frame(
      temp_vcf_object[,c('CHROM', 'POS', 'REF', 'ALT', 'FILTER','Type_of_polymorphism')],
      QUAL = as.numeric(temp_vcf_object$QUAL),
      QD = as.numeric(gsub('QD=','',str_extract(temp_vcf_object$INFO, 'QD=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
      SOR = as.numeric(gsub('SOR=','',str_extract(temp_vcf_object$INFO, 'SOR=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
      FS = as.numeric(gsub('FS=','',str_extract(temp_vcf_object$INFO, 'FS=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
      MQ = as.numeric(gsub('MQ=','',str_extract(temp_vcf_object$INFO, 'MQ=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
      MQRankSum = as.numeric(gsub('MQRankSum=','',str_extract(temp_vcf_object$INFO, 'MQRankSum=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?'))),
      ReadPosRankSum = as.numeric(gsub('ReadPosRankSum=','',str_extract(temp_vcf_object$INFO, 'ReadPosRankSum=-?(0|[1-9][0-9]*)(\\.[0-9]+)?([eE][+-]?[0-9]+)?')))
    )
    
    rm('temp_vcf_object')
    
    HardFilteringStats %<>% pivot_longer(cols = all_of(c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
                                        values_to = 'Value', names_to = 'Stat') %>%
      mutate(Stat = factor(Stat, levels = c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
             Type_of_polymorphism2 = case_when(
               grepl('SNP',Type_of_polymorphism) ~ 'SNP',
               !grepl('SNP',Type_of_polymorphism) ~ 'INDEL'),
             Thres = case_when(
               grepl('SNP',Type_of_polymorphism) & grepl('QUAL',Stat) ~ 30,
               grepl('SNP',Type_of_polymorphism) & grepl('QD',Stat) ~ 2,
               grepl('SNP',Type_of_polymorphism) & grepl('SOR',Stat) ~ 3,
               grepl('SNP',Type_of_polymorphism) & grepl('FS',Stat) ~ 60,
               grepl('SNP',Type_of_polymorphism) & grepl('MQ',Stat) ~ 40,
               grepl('SNP',Type_of_polymorphism) & grepl('MQRankSum',Stat) ~ -12.5,
               grepl('SNP',Type_of_polymorphism) & grepl('ReadPosRankSum',Stat) ~ 8,
               
               grepl('INDEL',Type_of_polymorphism) & grepl('QUAL',Stat) ~ 30,
               grepl('INDEL',Type_of_polymorphism) & grepl('QD',Stat) ~ 2,
               grepl('INDEL',Type_of_polymorphism) & grepl('SOR',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('FS',Stat) ~ 200,
               grepl('INDEL',Type_of_polymorphism) & grepl('MQ',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('MQRankSum',Stat) ~ NA,
               grepl('INDEL',Type_of_polymorphism) & grepl('ReadPosRankSum',Stat) ~ 8
             )
      )
    
    plot_HardFilteringStats_density =  HardFilteringStats %>%
      ggplot(aes(x = Value, fill = Type_of_polymorphism2, color = Type_of_polymorphism2)) +
      geom_density(alpha = .4) +
      geom_vline(data = data.frame(Type_of_polymorphism2 = rep(c('SNP', 'INDEL'), each = 7),
                                   Stat = factor(rep(c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum'), 2),
                                                 levels = c('QUAL', 'QD', 'SOR', 'FS', 'MQ', 'MQRankSum', 'ReadPosRankSum')),
                                   Thres = c(30, 2, 3, 60, 40, -12.5, -8, 
                                             30, 2, NA, 200, NA, NA, -20)),
                 aes(xintercept = Thres, color = Type_of_polymorphism2)
      ) + 
      facet_wrap(Stat ~., scale = 'free', ncol = 2) +
      theme_bw()+
      theme(legend.position = c(.7,.1))
    
    save.image('HardFilteringReport.RData', list = c('HardFilteringStats', 'plot_HardFilteringStats_density'))
    
  }
}


# Section 3: Pre filtering vcf file ----

## Step 1: First filter remove_filtered_all ----

if(remove_filtered_all){
  
  # rename output vcf
  output = paste0(output, '_filtered')
  
  # apply filter using vcftools
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_remove_filtered_all.sh',
               remove_filtered_all = remove_filtered_all,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  # rename the new input vcf file for further steps
  vcf_file = paste0(output, '.recode.vcf')
  
  print("removed_filtered_all done")
  system('echo "removed_filtered_all done"')
  
  }

## Step 2: Second filter of coding or non-coding regions----

if(!is.na(include_from_gff)){
  
  # rename output vcf
  
  output = paste0(output, '_PATTERN', include_from_gff, 'SELECTED')
  
  # get coding regions coordinates from the gff file
  ref_gff = ape::read.gff(ref_gff_file)
  
  ref_gff = ref_gff[grepl(include_from_gff, ref_gff$type)&
                      !grepl('^Transfer',ref_gff$seqid),]
  
  ref_gff = cbind(ref_gff, as.data.frame(t(sapply(1:nrow(ref_gff), function(gene){
    attributes = strsplit(ref_gff[gene,][['attributes']], ';')[[1]]
    c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
      gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
  }))))
  
  genomic_regions = ref_gff[,c('seqid', 'start', 'end')]
  
  genomic_regions = genomic_regions[order(genomic_regions$seqid),][order(genomic_regions$start),]
  
  rownames(genomic_regions) = 1:nrow(genomic_regions)
  
  write.table(genomic_regions, 'genomic_regions.bed', sep = '\t', quote = FALSE, row.names = FALSE)
  
  # apply filter using vcftools
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_FilteredGenomicRegions.sh',
               bed = 'genomic_regions.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  # rename the new input vcf file for further steps
  vcf_file = paste0(output, '.recode.vcf')
  
  # remove coordinates
  system('rm genomic_regions.bed')
  
  print("include_from_gff done")
  system('echo "include_from_gff done"')
  
}

if(!is.na(exclude_from_gff)){
  
  # rename output vcf
  output = paste0(output, '_PATTERN', exclude_from_gff, 'REMOVED')
  
  # get coding regions coordinates from the gff file
  ref_gff = ape::read.gff(ref_gff_file)
  
  ref_gff = ref_gff[grepl(exclude_from_gff, ref_gff$type)&
                      !grepl('^Transfer',ref_gff$seqid),]
  
  ref_gff = cbind(ref_gff, as.data.frame(t(sapply(1:nrow(ref_gff), function(gene){
    attributes = strsplit(ref_gff[gene,][['attributes']], ';')[[1]]
    c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
      gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
  }))))
  
  genomic_regions = ref_gff[,c('seqid', 'start', 'end')]
  
  genomic_regions = genomic_regions[order(genomic_regions$seqid),][order(genomic_regions$start),]
  
  rownames(genomic_regions) = 1:nrow(genomic_regions)
  
  write.table(genomic_regions, 'genomic_regions.bed', sep = '\t', quote = FALSE, row.names = FALSE)
  
  # apply filter using vcftools
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_FilteredGenomicRegions.sh',
               exclude_bed = 'genomic_regions.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  # rename the new input vcf file for further steps
  vcf_file = paste0(output, '.recode.vcf')
  
  # remove coordinates
  system('rm genomic_regions.bed')
  
  print("exclude_from_gff done")
  system('echo "exclude_from_gff done"')
  
}

## Step 3: Third filter of Core regions Only----

if(!is.na(bed)){
  
  print("Starting selection from bed")
  
  output = paste0(output, '_SelectedFromBed')
  
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_keep_bed.sh',
               bed = bed,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  vcf_file = paste0(output, '.recode.vcf')
  
  print("keep_bed done")
  system('echo "keep_bed done"')
  
}

if(!is.na(exclude_bed)){
  
  output = paste0(output, '_RemovedFromBed')
  
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_exclude_bed.sh',
               exclude_bed = exclude_bed,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  vcf_file = paste0(output, '.recode.vcf')
  
  print("exclude_bed done")
  system('echo "exclude_bed done"')
  
}

## Step 4: Fourth filter of only a subset of samples----

if(!is.na(keep_regexp)){
  
  output = paste0(output, '_SelectedSamplesOnly')
  
  run_vcftools(vcf = vcf_file,
               out = output,
               bash_file = 'run_vcf_SelectedSamplesOnly.sh',
               keep_regexp = keep_regexp,
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  vcf_file = paste0(output, '.recode.vcf')
  
  system('rm samples.indv')
  print("SelectedSamplesOnly done")
  system('echo "SelectedSamplesOnly done"')
  
}

## Step 4: Fifth filter of only a subset of Positions----

if(!is.na(positions)){
  
  print("Starting selection from positions")
  system('echo "Starting selection from positions"')
  
  system(paste0("head ", positions))
  
  output = paste0(output, '_SelectedFromPositions')
  
  run_vcftools(vcf = vcf_file,
               bash_file = 'run_selected_positions_filtering.sh',
               out = output,
               positions = positions,
               recode = TRUE,
               recode_INFO_all = TRUE
               )
  
  vcf_file = paste0(output, '.recode.vcf')
  
  print("positions done")
  system('echo "positions done"')
}



# Section 4: vcf2rGenome----

## If the number of task arrays is greater than one----

if(do_vcf2rGenome){
  
  if(nTasks > 1){ # If task arrays greater than one
    
    system(paste0('cp ', vcf_file, ' ', gsub('.vcf$', paste0('_', Task_id, '.vcf'), vcf_file)))
    
    vcf_file = gsub('.vcf$', paste0('_', Task_id, '.vcf'), vcf_file)
    
    ### Load vcf data----  
    # Split the genome in equal peaces
    
    files_nrows = as.integer(system(paste0("grep -v '^#' ", vcf_file, " | wc -l"), intern = TRUE))
    
    print(paste0('there are ', files_nrows, ' positions in the file'))
    
    s = round(seq(1, files_nrows + 1, length.out = nTasks + 1))
    
    print(paste0('the windows are :', paste(s, collapse = ',')))
    
    print(paste0('uploading window ', Task_id, ' of ', nTasks))
    
    vcf_start = s[Task_id]
    print(paste0('window ', Task_id, ' starts at postion ', vcf_start))
    
    vcf_end = s[Task_id + 1] - 1
    print(paste0('window ', Task_id, ' starts at postion ', vcf_end))
    
    
    # Upload the VCF to R environment
    temp_vcf_object = load_vcf(vcf = vcf_file, na.rm = TRUE, start = vcf_start, end = vcf_end)
    
    print(paste0('window ', Task_id, ' uploaded'))
    
    system(paste0('rm ', vcf_file))
    
    ### Data for coverage report----
    # Generate report of read depth per position and percentage of genome coverage by sample for a range of min read depth from 1 to ReadDepthThreshold
    if(coverage_data){
      
      #### vcf2rGenome, threshold == 1----
      
      print(paste0('converting window ', Task_id, ' to rGenome'))
      
      temp_rGenome_object = vcf2rGenome(vcf = temp_vcf_object, n = nchunks, threshold = 1)
      
      print(paste0('window ', Task_id, ' converted to rGenome'))
      
      #### Include metadata if available----
      if(!is.na(metadata)){
        
        print('reading and adding metadata')
        external_metadata = read.csv(metadata)
        temp_rGenome_object@metadata = left_join(temp_rGenome_object@metadata, external_metadata, by = c('Sample_id' = join_by))
        
      }
      
      #### Read depth summary----
      
      print(paste0('measuring read depth for window ', Task_id))
      
      Read_Depth_Summ = NULL
      
      for(w in 1:nchunks){
        Read_Depth_Summ = rbind(Read_Depth_Summ, data.frame(Summarise_ReadDepth(obj = temp_rGenome_object, by = pop, w = w, n = nchunks), thres = 1))
        rm('w')
      }
      
      assign(paste0('Read_Depth_Summ_Chunk', Task_id), Read_Depth_Summ)
      rm('Read_Depth_Summ')
      
      # Calculate the proportion of amplified loci (amplification rate), with a coverage equals or greater than 1, 2, 3, 4 and 5, of the samples.
      
      #### Sample amplification rate by different coverage thresholds----
      
      print(paste0('measuring sample amplification rate for window ', Task_id))
      
      SampAmpRate_Summ = NULL
      
      for(thres in 0:(ReadDepthThreshold - 1)){
        temp_SampAmpRate = SampleAmplRate(temp_rGenome_object, update = FALSE, type = 'count', threshold = ifelse(thres == 0, NA, thres))
        temp_SampAmpRate = data.frame(Sample_id = names(temp_SampAmpRate),
                                      NumberOfAmplifiedLoci = temp_SampAmpRate,
                                      Threshold = thres + 1)
        
        if(!is.na(pop)){
          temp_SampAmpRate = left_join(temp_SampAmpRate,
                                       temp_rGenome_object@metadata[,c('Sample_id', pop)],
                                       by = 'Sample_id')
        }
        
        
        SampAmpRate_Summ = rbind(SampAmpRate_Summ, temp_SampAmpRate)
        
        rm(list = c('temp_SampAmpRate', 'thres'))
        
      }
      
      assign(paste0('SampAmpRate_Summ_Chunk', Task_id), SampAmpRate_Summ)
      rm('SampAmpRate_Summ')
      
      #### Locus amplification rate by different coverage thresholds----
      
      print(paste0('measuring locus amplification rate for window ', Task_id))
      
      LocusAmpRate_Summ = NULL
      
      for(thres in 0:(ReadDepthThreshold - 1)){
        temp_LocusAmpRate = LocusAmplRate(obj = temp_rGenome_object, 
                                          update = FALSE, 
                                          by = pop, 
                                          threshold = ifelse(thres == 0, NA, thres))
        temp_LocusAmpRate = data.frame(Locus_id = ifelse(is.null(ncol(temp_LocusAmpRate)),
                                                         names(temp_LocusAmpRate),
                                                         rownames(temp_LocusAmpRate)),
                                       temp_LocusAmpRate,
                                       Threshold = thres + 1)
        
        LocusAmpRate_Summ = rbind(LocusAmpRate_Summ, temp_LocusAmpRate)
        
        rm(list = c('temp_LocusAmpRate', 'thres'))
        
      }
      
      assign(paste0('LocusAmpRate_Summ_Chunk', Task_id), LocusAmpRate_Summ)
      rm('LocusAmpRate_Summ')
      
      ### Mask alleles below the ReadDepthThreshold----
      
      print(paste0('pruning alleles bellow threshold for window ', Task_id))
      
      temp_rGenome_object@gt = prune_alleles(temp_rGenome_object, threshold = ReadDepthThreshold - 1, n = nchunks)
      
      assign(rGenome_object_name, temp_rGenome_object)
      rm('temp_rGenome_object')
      
      save(file = imagename, list = ls()[grep(paste0(rGenome_object_name,'|LocusAmpRate_Summ|SampAmpRate_Summ|Read_Depth_Summ'),ls())])
      
      print(paste0('image saved for window ', Task_id))
      
    }else{
      
      #### vcf2rGenome----
      temp_rGenome_object = vcf2rGenome(vcf = temp_vcf_object, n = nchunks, threshold = ReadDepthThreshold)
      assign(rGenome_object_name, temp_rGenome_object)
      rm('temp_rGenome_object')
      
      save(file = imagename, list = ls()[grep(rGenome_object_name,ls())])
      
      print(paste0('image saved for window ', Task_id))
    }
    
  }else{
    
    ## Only one task array----
    
    ### Load vcf data----  
    
    # Upload the VCF to R environment
    vcf_object = load_vcf(vcf = vcf_file, na.rm = TRUE)
    
    ### Data for coverage report----
    # Generate report of read depth per position and percentage of genome coverage by sample for a range of min read depth from 1 to ReadDepthThreshold
    if(coverage_data){
      
      #### vcf2rGenome, threshold == 1----
      rGenome_object = vcf2rGenome(vcf = vcf_object, n = nchunks, threshold = 1)
      
      #### Include metadata if available----
      if(!is.na(metadata)){
        external_metadata = read.csv(metadata)
        rGenome_object@metadata = left_join(rGenome_object@metadata, external_metadata, by = c('Sample_id' = join_by))
        
      }
      
      #### Read depth summary----
      Read_Depth_Summ = NULL
      
      for(w in 1:nchunks){
        Read_Depth_Summ = rbind(Read_Depth_Summ, data.frame(Summarise_ReadDepth(obj = rGenome_object, by = pop, w = w, n = nchunks), thres = 1))
        rm('w')
      }
      
      # Calculate the proportion of amplified loci (amplification rate), with a coverage equals or greater than 1, 2, 3, 4 and 5, of the samples.
      
      
      #### Sample amplification rate by different coverage thresholds----
      
      SampAmpRate_Summ = NULL
      
      for(thres in 0:(ReadDepthThreshold - 1)){
        temp_SampAmpRate = SampleAmplRate(rGenome_object, update = FALSE, type = 'count', threshold = ifelse(thres == 0, NA, thres))
        temp_SampAmpRate = data.frame(Sample_id = names(temp_SampAmpRate),
                                      NumberOfAmplifiedLoci = temp_SampAmpRate,
                                      Threshold = thres + 1)
        
        if(!is.na(pop)){
          temp_SampAmpRate = left_join(temp_SampAmpRate,
                                       rGenome_object@metadata[,c('Sample_id', pop)],
                                       by = 'Sample_id')
        }

        SampAmpRate_Summ = rbind(SampAmpRate_Summ, temp_SampAmpRate)
        
        rm(list = c('temp_SampAmpRate', 'thres'))
        
      }
      
      ## Plot SampAmpRate distribution----
      
      if(!is.na(pop)){
        
        ### Edit SampAmpRate_Summ
        
        colnames(SampAmpRate_Summ) = c('Sample_id', 'NumberOfAmplifiedLoci', 'Threshold', 'Pop')
      
        SampAmpRate_Summ %<>% 
          select(Sample_id, Threshold, NumberOfAmplifiedLoci, Pop)
        
        SampAmpRate_Summ %<>% mutate(AmpRate = NumberOfAmplifiedLoci/nrow(rGenome_object@loci_table))
        
        plot_amplifiedsamples_amplifiedloci_curve = SampAmpRate_Summ %>% 
          group_by(Threshold, Pop)%>%
          summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
                    AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
                    AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
                    AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
                    AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
                    AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
                    AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
                    AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
                    AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
                    AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
                    AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
                    AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
                    AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
                    AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
                    AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
                    AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
                    AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
                    AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
                    AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
                    AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1)
                    ) %>%
          pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                       values_to = 'Percentage',
                       names_to = 'AmpRate') %>%
          mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
          ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
          geom_line() +
          geom_vline(xintercept = sample_ampl_rate*100, linetype = 2) +
          facet_grid(.~Pop) +
          theme_bw() +
          labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
        
        
        
      }else{
        
        SampAmpRate_Summ %<>% 
          select(Sample_id, Threshold, NumberOfAmplifiedLoci)
        
        SampAmpRate_Summ %<>% mutate(AmpRate = NumberOfAmplifiedLoci/nrow(rGenome_object@loci_table))
        
        plot_amplifiedsamples_amplifiedloci_curve = SampAmpRate_Summ %>% 
          group_by(Threshold)%>%
          summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
                    AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
                    AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
                    AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
                    AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
                    AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
                    AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
                    AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
                    AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
                    AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
                    AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
                    AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
                    AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
                    AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
                    AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
                    AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
                    AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
                    AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
                    AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
                    AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1)
                    ) %>%
          pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                       values_to = 'Percentage',
                       names_to = 'AmpRate') %>%
          mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
          ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
          geom_line() +
          geom_vline(xintercept = sample_ampl_rate*100, linetype = 2) +
          theme_bw() +
          labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Coverage')
      }
      
      #### Locus amplification rate by different coverage thresholds----
      
      LocusAmpRate_Summ = NULL
      
      for(thres in 0:(ReadDepthThreshold - 1)){
        temp_LocusAmpRate = LocusAmplRate(obj = rGenome_object, update = FALSE, by = pop, threshold = ifelse(thres == 0, NA, thres))
        temp_LocusAmpRate = data.frame(Locus_id = ifelse(is.null(ncol(temp_LocusAmpRate)),
                                                         names(temp_LocusAmpRate),
                                                         rownames(temp_LocusAmpRate)),
                                       temp_LocusAmpRate,
                                       Threshold = thres + 1)
        
        LocusAmpRate_Summ = rbind(LocusAmpRate_Summ, temp_LocusAmpRate)
        
        rm(list = c('temp_LocusAmpRate', 'thres'))
        
      }
      
      ### Mask alleles below the ReadDepthThreshold----
      
      rGenome_object@gt = prune_alleles(rGenome_object, threshold = ReadDepthThreshold - 1, n = nchunks)
      
      assign(rGenome_object_name, rGenome_object)
      rm('rGenome_object')
      
      save(file = imagename, list = ls()[grep(paste0(rGenome_object_name,'|LocusAmpRate_Summ|SampAmpRate_Summ|Read_Depth_Summ'),ls())])
      
      print(paste0('image saved'))
      
    }else{
      
      ### vcf2rGenome----
      rGenome_object = vcf2rGenome(vcf = vcf_object, n = nchunks, threshold = ReadDepthThreshold)
      assign(rGenome_object_name, rGenome_object)
      rm('rGenome_object')
      
      save(file = imagename, list = ls()[grep(rGenome_object_name,ls())])
      
      print(paste0('image saved'))
      
    }
    
  }
  
}


# Section 5: Merge chunks----

if(merge_rgenome){
  
  ## Load rGenome chunks----
  rGenome_objects = NULL
  LocusAmpRate_Summ = NULL
  Read_Depth_Summ = NULL
  SampAmpRate_Summ = NULL
  for(file in list.files('Chunks/', pattern = '.RData')){
    load(file.path('Chunks/', file))
    print(paste0('enviroment for file', file, ' loaded'))
    
    #print(paste(ls(), collapse = ', '))
    
    rGenome_objects[[gsub('Chunks/|\\.RData', '', file)]] = get(gsub('Chunks/|\\.RData', '', file))
    print('rGenome_object stored')
    
    if(paste0('LocusAmpRate_Summ', str_extract(file, '_Chunk\\d+')) %in% ls()){
      LocusAmpRate_Summ = rbind(LocusAmpRate_Summ, get(paste0('LocusAmpRate_Summ', str_extract(file, '_Chunk\\d+'))))
      print('LocusAmpRate_Summ merged')
    }else{
      print('LocusAmpRate_Summ not found')
    }
    
    if(paste0('Read_Depth_Summ', str_extract(file, '_Chunk\\d+')) %in% ls()){
      Read_Depth_Summ = rbind(Read_Depth_Summ, get(paste0('Read_Depth_Summ', str_extract(file, '_Chunk\\d+'))))
      print('Read_Depth_Summ merged')
    }else{
      print('Read_Depth_Summ not found')
    }
    
    if(paste0('SampAmpRate_Summ', str_extract(file, '_Chunk\\d+')) %in% ls()){
      
      temp_SampAmpRate_Summ = get(paste0('SampAmpRate_Summ', str_extract(file, '_Chunk\\d+')))
      
      temp_SampAmpRate_Summ %<>% mutate(row_id = paste(Sample_id, Threshold, sep = '_'))
      
      if(!is.na(pop)){
        temp_SampAmpRate_Summ = temp_SampAmpRate_Summ[,c("row_id", "NumberOfAmplifiedLoci", pop)]
        colnames(temp_SampAmpRate_Summ) = c("row_id", "NumberOfAmplifiedLoci", "Pop")
        
      }else{
        temp_SampAmpRate_Summ = temp_SampAmpRate_Summ[,c("row_id", "NumberOfAmplifiedLoci")]
      }
      
      
      if(is.null(SampAmpRate_Summ)){
        
        SampAmpRate_Summ = temp_SampAmpRate_Summ
        print(paste(colnames(SampAmpRate_Summ), collapse = ', '))
        print('SampAmpRate_Summ merged')
        
      }else{
        
        if(!is.na(pop)){
          SampAmpRate_Summ = rbind(SampAmpRate_Summ, temp_SampAmpRate_Summ) %>%
            group_by(row_id, Pop) %>%
            dplyr::summarise(NumberOfAmplifiedLoci = sum(NumberOfAmplifiedLoci))
          
          print(paste(colnames(SampAmpRate_Summ), collapse = ', '))
          
          SampAmpRate_Summ %<>% dplyr::select(row_id, NumberOfAmplifiedLoci, Pop)
          
          print('SampAmpRate_Summ merged')
          
        }else{
          
          SampAmpRate_Summ = rbind(SampAmpRate_Summ, temp_SampAmpRate_Summ) %>% 
            group_by(row_id)%>%
            dplyr::summarise(NumberOfAmplifiedLoci = sum(NumberOfAmplifiedLoci))
          
          print(paste(colnames(SampAmpRate_Summ), collapse = ', '))
          SampAmpRate_Summ %<>% dplyr::select(row_id, NumberOfAmplifiedLoci)
          print('SampAmpRate_Summ merged')
        }
        
      }
      
    }else{
      print('SampAmpRate_Summ not found')
    }
    
    print(paste0(file, ' uploaded and merged'))
    
    rm(list = gsub('Chunks/|\\.RData', '', file))
    rm(list = paste0('LocusAmpRate_Summ', str_extract(file, '_Chunk\\d+')))
    rm(list = paste0('Read_Depth_Summ', str_extract(file, '_Chunk\\d+')))
    rm(list = paste0('SampAmpRate_Summ', str_extract(file, '_Chunk\\d+')))
    rm(list = c('file', 'temp_SampAmpRate_Summ'))
  }
  
  ## Combine rGenome objects----
  
  rGenome_object = rGenome(gt = NULL,
                           loci_table = NULL,
                           metadata = NULL)
  
  rGenome_object@metadata = rGenome_objects[[1]]@metadata
  
  for(obj in names(rGenome_objects)){
    rGenome_object@gt = rbind(rGenome_object@gt, rGenome_objects[[obj]]@gt)
    rGenome_object@loci_table = rbind(rGenome_object@loci_table, rGenome_objects[[obj]]@loci_table)
    rGenome_objects[[obj]] = NULL
    rm('obj')
  }
  
  print('rGenome_objects merged')
  
  rm('rGenome_objects') 
  
  ## Plot SampAmpRate distribution----

  
  if(!is.na(pop) & !is.null(SampAmpRate_Summ)){
    
    ### Edit SampAmpRate_Summ
    
    SampAmpRate_Summ %<>% mutate(Sample_id = gsub('_\\d+$', '', row_id),
                                 Threshold = str_extract(row_id, '\\d+$')) %>% select(Sample_id, Threshold, NumberOfAmplifiedLoci, Pop)
    
    SampAmpRate_Summ %<>% mutate(AmpRate = NumberOfAmplifiedLoci/nrow(rGenome_object@loci_table))
    
    plot_amplifiedsamples_amplifiedloci_curve = SampAmpRate_Summ %>% 
      group_by(Threshold, Pop)%>%
      summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
                AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
                AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
                AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
                AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
                AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
                AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
                AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
                AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
                AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
                AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
                AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
                AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
                AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
                AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
                AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
                AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
                AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
                AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
                AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1)
                ) %>%
      pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                   values_to = 'Percentage',
                   names_to = 'AmpRate') %>%
      mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
      ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
      geom_line() +
      geom_vline(xintercept = sample_ampl_rate*100, linetype = 2) +
      facet_grid(.~Pop) +
      theme_bw() +
      labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Read Depth')
    
    SampAmpRate_Summ = SampAmpRate_Summ[,c('Sample_id', 'Threshold', 'NumberOfAmplifiedLoci', 'AmpRate', 'Pop')]
    colnames(SampAmpRate_Summ) = c('Sample_id', 'Threshold', 'NumberOfAmplifiedLoci', 'AmpRate', pop)
    
  }else if(!is.null(SampAmpRate_Summ)){
    
    ### Edit SampAmpRate_Summ
    
    SampAmpRate_Summ %<>% mutate(Sample_id = gsub('_\\d+$', '', row_id),
                                 Threshold = str_extract(row_id, '\\d+$')) %>% select(Sample_id, Threshold, NumberOfAmplifiedLoci)
    
    SampAmpRate_Summ %<>% mutate(AmpRate = NumberOfAmplifiedLoci/nrow(rGenome_object@loci_table))
    
    plot_amplifiedsamples_amplifiedloci_curve = SampAmpRate_Summ %>% 
      group_by(Threshold)%>%
      summarise(AmpRate5 = round(100*sum(AmpRate >= .05)/n(), 1),
                AmpRate10 = round(100*sum(AmpRate >= .10)/n(), 1),
                AmpRate15 = round(100*sum(AmpRate >= .15)/n(), 1),
                AmpRate20 = round(100*sum(AmpRate >= .20)/n(), 1),
                AmpRate25 = round(100*sum(AmpRate >= .25)/n(), 1),
                AmpRate30 = round(100*sum(AmpRate >= .30)/n(), 1),
                AmpRate35 = round(100*sum(AmpRate >= .35)/n(), 1),
                AmpRate40 = round(100*sum(AmpRate >= .40)/n(), 1),
                AmpRate45 = round(100*sum(AmpRate >= .45)/n(), 1),
                AmpRate50 = round(100*sum(AmpRate >= .50)/n(), 1),
                AmpRate55 = round(100*sum(AmpRate >= .55)/n(), 1),
                AmpRate60 = round(100*sum(AmpRate >= .60)/n(), 1),
                AmpRate65 = round(100*sum(AmpRate >= .65)/n(), 1),
                AmpRate70 = round(100*sum(AmpRate >= .70)/n(), 1),
                AmpRate75 = round(100*sum(AmpRate >= .75)/n(), 1),
                AmpRate80 = round(100*sum(AmpRate >= .80)/n(), 1),
                AmpRate85 = round(100*sum(AmpRate >= .85)/n(), 1),
                AmpRate90 = round(100*sum(AmpRate >= .90)/n(), 1),
                AmpRate95 = round(100*sum(AmpRate >= .95)/n(), 1),
                AmpRate100 = round(100*sum(AmpRate >= 1)/n(), 1)
                ) %>%
      pivot_longer(cols = paste0('AmpRate', seq(5, 100, 5)),
                   values_to = 'Percentage',
                   names_to = 'AmpRate') %>%
      mutate(AmpRate = as.numeric(gsub('AmpRate','', AmpRate)))%>%
      ggplot(aes(x = AmpRate, y = Percentage, color = as.factor(Threshold), group = as.factor(Threshold))) +
      geom_line() +
      geom_vline(xintercept = sample_ampl_rate*100, linetype = 2) +
      theme_bw() +
      labs(x = '% of amplified loci (amplification rate)', y = '% of Samples', color = 'Min Read Depth')
  }
  
  ## Save image----
  assign(rGenome_object_name, rGenome_object)
  rm('rGenome_object')
  
  save(file = imagename, list = ls()[grep(paste0(rGenome_object_name,'|LocusAmpRate_Summ|SampAmpRate_Summ|Read_Depth_Summ|plot_amplifiedsamples_amplifiedloci_curve'),ls())])
  
  print('RData image saved')
  
}

# Section 6: Post-Filtering----

if(post_filtering){
  
  load(imagename)
  
  print('image for post-filtering loaded')
  
  assign('rGenome_object', get(rGenome_object_name))
  rm(list = c(rGenome_object_name))
  
  ## Setp 1: Remove samples with less than sample_ampl_rate% of amplified loci----
  
  rGenome_object = SampleAmplRate(rGenome_object, update = TRUE)
  
  removed_samples_rGenome = filter_samples(obj = rGenome_object, 
                                  v = rGenome_object@metadata$SampleAmplRate < sample_ampl_rate)
  
  rGenome_object = filter_samples(obj = rGenome_object, 
                                  v = rGenome_object@metadata$SampleAmplRate >= sample_ampl_rate)
  
  print(paste0('removal of samples with amplification rate below ', sample_ampl_rate))
  
  ### Update allele counts----
  
  allele_counts = NULL
  
  for(w in 1:nchunks){
    allele_counts = rbind(allele_counts, get_AC(obj = rGenome_object, w = w, n = nchunks))
    rm(w)
  }
  
  rGenome_object@loci_table = cbind(rGenome_object@loci_table,
                                    allele_counts)
  rm('allele_counts')
  
  print('update allele counts and allele frequencies done')
  
  ### Filter monomorphic sites respect to the reference----
  
  removed_monomorphic_sites_rGenome = filter_loci(rGenome_object,
                               v = (rGenome_object@loci_table$Cardinality <= 1 &
                                       grepl('^0',rGenome_object@loci_table$Allele_Counts)))
  
  rGenome_object = filter_loci(rGenome_object,
                               v = !(rGenome_object@loci_table$Cardinality <= 1 &
                                       grepl('^0',rGenome_object@loci_table$Allele_Counts)))
  
  print('removal of no polymorphic sites done')
  
  ## Step 2: Remove loci with high missing data ----
  
  rGenome_object = LocusAmplRate(rGenome_object)
  
  ### Filter loci with amplification rate below locus_ampl_rate% ----
  
  removed_sites_LowAmplRate_rGenome = filter_loci(rGenome_object, v = rGenome_object@loci_table$LocusAmplRate < locus_ampl_rate)
  rGenome_object = filter_loci(rGenome_object, v = rGenome_object@loci_table$LocusAmplRate >= locus_ampl_rate)
  
  print(paste0('removal of loci with amplification rate below ', locus_ampl_rate))
  
  ### Update alternative alleles----
  
  rGenome_object@loci_table$ALT =
    ifelse(rGenome_object@loci_table$REF !=
             gsub(',([ATCG]|\\*)+',
                  '',
                  gsub(':\\d+',
                       '',
                       rGenome_object@loci_table$Alleles)),
           gsub(':\\d+',
                '',
                rGenome_object@loci_table$Alleles),
           gsub('^([ATCG]|\\*)+,',
                '',
                gsub(':\\d+', '', rGenome_object@loci_table$Alleles))
    )
  
  print('update list of alternative alleles')
  
  ## Step 3: Filter of potential PCR genotyping artifacts (Homopolymers & Short Tandem Repeats)----
  
  ### Classify SNPs, INDELS, Homopolymers and Short tandem repeats----
  
  Type_of_polymorphism = NULL
  
  for(w in 1:nchunks){
    Type_of_polymorphism = c(Type_of_polymorphism, get_type_of_polymorphism(rGenome_object, w = w, n = nchunks))
    rm(w)
  }
  
  rGenome_object@loci_table$Type_of_polymorphism = Type_of_polymorphism
  rm(Type_of_polymorphism)
  
  print('diferentiation of SNPs and INDELs done')
  
  ### Effect of PCR genotyping artifacts on heterozygosity----
  
  #### Proportion of Heterozygous samples per site
  
  ObsHet = NULL
  
  for(w in 1:nchunks){
    
    ObsHet = c(ObsHet, get_ObsHet(rGenome_object, by = 'loci', w = w, n = nchunks))
    rm(w)
  }
  rGenome_object@loci_table$ObsHet = ObsHet
  rm(ObsHet)
  
  #### Fraction of Heterozygous samples per Alternative alleles per site
  
  frac_ofHet_pAlts = NULL
  
  for(w in 1:nchunks){
    frac_ofHet_pAlts = c(frac_ofHet_pAlts, frac_ofHet_pAlt(rGenome_object, w = w, n = nchunks))
    rm(w)
  }
  
  rGenome_object@loci_table$frac_ofHet_pAlts = frac_ofHet_pAlts
  rm(frac_ofHet_pAlts)
  
  #### Distribution of the fraction of alternative alleles in polyclonal samples
  plot_frac_ofHet_pAlt_histogram = rGenome_object@loci_table %>% ggplot(aes(x = frac_ofHet_pAlts))+
    geom_histogram(binwidth = 0.01)+
    labs(x = 'Fraction of heterozygous samples per alternative allele per site',
         y = 'Number of sites (Loci)')+
    theme_bw()
  
  #### Classify sites based on the proportion of alternative alleles in polyclonal samples per site
  rGenome_object@loci_table %<>% mutate(ALT_FILTER =case_when(
    frac_ofHet_pAlts == 1 ~ '100%',
    frac_ofHet_pAlts < 1 & frac_ofHet_pAlts > .5 ~ '50 - 99%',
    frac_ofHet_pAlts <= .5 ~ '<=50%',
  ))
  
  print('Classification of sites based on the proportion of alternative alleles in polyclonal samples per site done')
  
  #### Distribution of observed heterozygosity per locus, by type of marker, by group (defined in previuos step 10.2)
  
  plot_obHet_histogram = rGenome_object@loci_table %>%
    ggplot(aes(x = ObsHet,
               fill = factor(ALT_FILTER,
                             levels = c('<=50%', '50 - 99%', '100%'))))+
    geom_histogram(binwidth = .01)+
    scale_fill_manual(values = c('dodgerblue3', 'gold3', 'firebrick3'))+
    facet_wrap(.~factor(Type_of_polymorphism, levels = c(
      'SNP',
      'INDEL',
      'INDEL:Homopolymer',
      'INDEL:Dinucleotide_STR',
      'INDEL:Trinucleotide_STR',
      'INDEL:Tetranucleotide_STR',
      'INDEL:Pentanucleotide_STR',
      'INDEL:Hexanucleotide_STR'
    )), scales = 'free_y', ncol = 4)+
    labs(y = 'Number of Loci',
         x = 'Observed Heterozygosity per locus',
         fill = 'Het/Alt')+
    theme_bw()
  
  #### Remove Homopolymers and DiSTRs----
  
  if(sum(!is.na(type_of_polymorphism_to_remove)) > 0){
    
    removed_polymorphism_rGenome = filter_loci(rGenome_object,
                                 v = (rGenome_object@loci_table$Type_of_polymorphism %in%
                                         type_of_polymorphism_to_remove))
    
    rGenome_object = filter_loci(rGenome_object,
                                 v = !(rGenome_object@loci_table$Type_of_polymorphism %in%
                                         type_of_polymorphism_to_remove))
  }
  
  
  print('removal of Homopolymers and Dinucleotide STRs done')
  
  ## Step 4: Filter of problematic genomic regions to map the reads ----
  
  ### Regions with excess of Heterozygosity----
  
  #### Distribution of proportion of heterozygous per site
  manhatan_plot_ObsHet_by_AltFilter = rGenome_object@loci_table %>%
    mutate(Type_of_polymorphism2 = case_when(
      Type_of_polymorphism == 'SNP' ~ 'SNP',
      Type_of_polymorphism != 'SNP' ~ 'INDELs'),
      CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
    ggplot(aes(x = POS, y = ObsHet, color = factor(ALT_FILTER, levels = c('<=50%', '50 - 99%', '100%')))) +
    geom_point(alpha = 0.5, size = .25) +
    scale_color_manual(values = c('dodgerblue3', 'gold3', 'firebrick3'))+
    facet_grid(Type_of_polymorphism2 ~ CHROM2, scales = 'free_x', space = 'free_x')+
    theme_bw()+
    labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'Het/Alt')+
    theme(legend.position = 'right',
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  #### Average of heterozygous samples per site per gene
  rGenome_object@loci_table = mean_ObsHet(rGenome_object, 
                                          gff = ref_gff_file, 
                                          Type_of_polymorphism = 'SNP')
  
  
  print('mean observed heterozygosity by gene calculation done')
  
  #### Define a threshold for trusted variants
  ObsHet_Threshold = quantile(rGenome_object@loci_table %>%
                                filter(ALT_FILTER == '<=50%', Type_of_polymorphism == 'SNP') %>%
                                group_by(gene_id) %>%
                                dplyr::summarise(mean_ObsHet = max(mean_ObsHet, na.rm = T)) %>%
                                select(mean_ObsHet) %>%
                                unlist, ObsHet_quantile)
  print('ObsHet threshold done')
  
  #### Distribution of Average Observed Heteozygosity in SNPs per genomic position
  
  plot_mean_obHet_histogram = rGenome_object@loci_table %>%
    filter(ALT_FILTER == '<=50%', Type_of_polymorphism == 'SNP') %>%
    group_by(gene_id) %>%
    dplyr::summarise(mean_ObsHet = max(mean_ObsHet, na.rm = T)) %>%
    ggplot(aes(x = mean_ObsHet))+
    geom_histogram(binwidth = .005)+
    geom_vline(xintercept = ObsHet_Threshold)+
    labs(y = 'Number of genomic regions', x = 'Average Observed Heteozygosity in SNPs') +
    theme_bw()
  
  #### Differentiate genomic regions that PASS the threshold----
  
  rGenome_object@loci_table %<>% mutate(ObsHet_Filter = case_when(
    mean_ObsHet > ObsHet_Threshold ~ 'OUT',
    mean_ObsHet <= ObsHet_Threshold ~ 'PASS'
  ))
  
  #### Visualize which genomic regions PASS the filter
  
  manhatan_plot_ObsHet_by_ObsHetFilter = rGenome_object@loci_table %>%
    mutate(Type_of_polymorphism2 = case_when(
      Type_of_polymorphism == 'SNP' ~ 'SNP',
      Type_of_polymorphism != 'SNP' ~ 'INDEL'),
      CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
    ggplot(aes(x = POS, y = ObsHet, color = ObsHet_Filter)) +
    geom_point(alpha = 0.5, size = .25) +
    geom_hline(yintercept = ObsHet_Threshold)+
    scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
    facet_grid(Type_of_polymorphism2 ~ CHROM2, scales = 'free_x', space = 'free_x')+
    theme_bw()+
    labs(y = 'Observed Heterozygosity',
         x = 'Chromosomal position',
         color = 'Mean ObsHet > Th')+
    theme(legend.position = 'right',
          axis.text.x = element_blank())
  
  print('identification of regions (genes) with excess of observed heterozygosity done')
  
  ### Regions with excess of density of Polymorphism----
  
  #### Calculate the density of SNPs per genomic region and define a threshold----
  rGenome_object@loci_table = SNP_density(rGenome_object, gff = ref_gff_file)
  
  print('SNP density calculation done')
  
  print(SNP_density_quantile)
  
  SNP_density_threshold = quantile(rGenome_object@loci_table %>%
                                     filter(Type_of_polymorphism == 'SNP') %>%
                                     group_by(gene_id) %>%
                                     dplyr::summarise(SNP_density = max(SNP_density, na.rm = T)) %>%
                                     select(SNP_density) %>%
                                     unlist(), probs = SNP_density_quantile)
  
  print('SNP density threshold done')
  #### Distribution of SNP density
  
  plot_SNPdensity_histogram = rGenome_object@loci_table %>%
    filter(Type_of_polymorphism == 'SNP') %>%
    group_by(gene_id) %>%
    dplyr::summarise(SNP_density = max(SNP_density, na.rm = T)) %>%
    ggplot(aes(x = SNP_density))+
    geom_histogram(binwidth = .001)+
    geom_vline(xintercept = SNP_density_threshold)+
    labs(y = 'Number of genomic regions', x = 'SNP density') +
    theme_bw()
  
  #### Differentiate/Filter genomic regions that PASS the threshold----
  
  rGenome_object@loci_table %<>% 
    mutate(SNP_density_Filter = case_when(
      SNP_density > SNP_density_threshold ~ 'OUT',
      SNP_density <= SNP_density_threshold ~ 'PASS'
    ))
  
  #### Visualize which genomic regions PASS the filter
  
  manhatan_plot_ObsHet_by_SNPdensityFilter = rGenome_object@loci_table %>%
    mutate(Type_of_polymorphism2 = case_when(
      Type_of_polymorphism == 'SNP' ~ 'SNP',
      Type_of_polymorphism != 'SNP' ~ 'INDEL'),
      CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
    ggplot(aes(x = POS, y = ObsHet, color = SNP_density_Filter)) +
    geom_point(alpha = 0.5, size = .25) +
    scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
    facet_grid(Type_of_polymorphism2 ~ CHROM2, scales = 'free_x', space = 'free_x')+
    theme_bw()+
    labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'SNP density')+
    theme(legend.position = 'right',
          axis.text.x = element_blank())
  
  print('identification of regions (genes) with excess of SNP density done')
  
  ### Differentiate of problematic genomic regions to map reads----
  
  if(alignment_filter == 'het'){
    
    rGenome_object@loci_table %<>%
      mutate(Alignment_Filter = case_when(
        ObsHet_Filter == 'OUT' ~ 'OUT',
        ObsHet_Filter == 'PASS' ~ 'PASS'
      ))
    
  }else if(alignment_filter == 'snp_dens'){
    
    rGenome_object@loci_table %<>%
      mutate(Alignment_Filter = case_when(
        SNP_density_Filter == 'OUT' ~ 'OUT',
        SNP_density_Filter == 'PASS' ~ 'PASS'
      ))
    
  }else if(alignment_filter == 'and'){
    
    rGenome_object@loci_table %<>%
      mutate(Alignment_Filter = case_when(
        ObsHet_Filter == 'OUT' & SNP_density_Filter == 'OUT' ~ 'OUT',
        !(ObsHet_Filter == 'OUT' & SNP_density_Filter == 'OUT') ~ 'PASS'
      ))
    
  }else if(alignment_filter == 'or'){
    
    rGenome_object@loci_table %<>%
      mutate(Alignment_Filter = case_when(
        ObsHet_Filter == 'OUT' | SNP_density_Filter == 'OUT' ~ 'OUT',
        !(ObsHet_Filter == 'OUT' | SNP_density_Filter == 'OUT') ~ 'PASS'
      ))
    
  }else if(alignment_filter == 'none'){
    
    rGenome_object@loci_table %<>%
      mutate(Alignment_Filter =  'PASS')
    
  }
  
  
  
  #### Visualize which genomic regions PASS the filter
  
  manhatan_plot_ObsHet_by_ObsHetSNPdensityFilter = rGenome_object@loci_table %>%
    mutate(Type_of_polymorphism2 = case_when(
      Type_of_polymorphism == 'SNP' ~ 'SNP',
      Type_of_polymorphism != 'SNP' ~ 'INDEL'),
      CHROM2  = gsub('(^(\\d|P|v)+_|_v1)', '',CHROM))%>%
    ggplot(aes(x = POS, y = ObsHet, color = Alignment_Filter)) +
    geom_point(alpha = 0.5, size = .25) +
    scale_color_manual(values = c('firebrick2', 'dodgerblue3'))+
    facet_grid(Type_of_polymorphism2 ~ CHROM2, scales = 'free_x', space = 'free_x')+
    theme_bw()+
    labs(y = 'Observed Heterozygosity', x = 'Chromosomal position', color = 'Alignment\nFilter')+
    theme(legend.position = 'right',
          axis.text.x = element_blank())
  
  #### Removed Genomic regions or genes ----
  
  removed_genomicregions = 
    filter_loci(rGenome_object,
                v = rGenome_object@loci_table$Alignment_Filter == 'OUT')
  
  print("removal of regions (genes) that didn't map correctly done")
  
  #### Remove problematic genomic regions to map reads ----
  
  rGenome_object_filtered = 
    filter_loci(rGenome_object,
                v = rGenome_object@loci_table$Alignment_Filter == 'PASS')
  
  ##  Step 5: Mask stochastic PCR errors and de-novo mutations----
  
  ### Calculate allele counts and frac_ofHet_pAlt by each allele----
  allele_count_frac_ofHet_pAlt = NULL
  
  for(w in 1:nchunks){

    allele_count_frac_ofHet_pAlt = rbind(allele_count_frac_ofHet_pAlt,
                                         frac_ofHet_pAlt_byAllele(
                                           rGenome_object_filtered,
                                           w = w,
                                           n = nchunks,
                                           add_variable = c('Type_of_polymorphism')))

  }

  
  
  plot_AlleleCount_HetCount = allele_count_frac_ofHet_pAlt %>%
    ggplot(aes(x = P_ij, 
               y = H_ij,
               color = h_ijminor))+
    geom_point(alpha = .05, size = .5)+
    facet_wrap(~factor(Type_of_polymorphism,
                       levels = c('SNP',
                                  'INDEL',
                                  'INDEL:Trinucleotide_STR',
                                  'INDEL:Tetranucleotide_STR',
                                  'INDEL:Pentanucleotide_STR',
                                  'INDEL:Hexanucleotide_STR')))+
    theme_bw()+
    scale_color_continuous(type = 'viridis')+
    labs(x = 'Samples with the alt. allele (P_ij)',
         y = 'Het. samples with the alt. allele (H_ij)',
         color = 'h_ijminor')
  
  plot_pij_hij_hijminor = allele_count_frac_ofHet_pAlt %>%
    ggplot(aes(x = p_ij, 
               y = h_ij,
               color = h_ijminor))+
    geom_point(alpha = .1, size = 0.5)+
    facet_wrap(~factor(Type_of_polymorphism,
                       levels = c('SNP',
                                  'INDEL',
                                  'INDEL:Trinucleotide_STR',
                                  'INDEL:Tetranucleotide_STR',
                                  'INDEL:Pentanucleotide_STR',
                                  'INDEL:Hexanucleotide_STR')))+
    theme_bw()+
    scale_color_continuous(type = 'viridis')+
    labs(x = 'Alternative allele frequency (p_ij)',
         y = 'h_ij (H_ij/P_ij)',
         color = 'h_ijminor')
  
  
  plot_AlleleFreq_histogram = allele_count_frac_ofHet_pAlt %>%
    mutate(h_ijminor_cat = case_when(
      h_ijminor >= .9 ~ '0.9 - 1.0',
      h_ijminor < .9 & h_ijminor >= .8 ~ '0.8 - 0.9',
      h_ijminor < .8 & h_ijminor >= .7 ~ '0.7 - 0.8',
      h_ijminor < .7 & h_ijminor >= .6 ~ '0.6 - 0.7',
      h_ijminor < .6 & h_ijminor >= .5 ~ '0.5 - 0.6',
      h_ijminor < .5 & h_ijminor >= .4 ~ '0.4 - 0.5',
      h_ijminor < .4 & h_ijminor >= .3 ~ '0.3 - 0.4',
      h_ijminor < .3 & h_ijminor >= .2 ~ '0.2 - 0.3',
      h_ijminor < .2 & h_ijminor >= .1 ~ '0.1 - 0.2',
      h_ijminor < .1 ~ '0.0 - 0.1'
    ),
    h_ijminor_cat = factor(h_ijminor_cat, levels = c(
      '0.9 - 1.0',
      '0.8 - 0.9',
      '0.7 - 0.8',
      '0.6 - 0.7',
      '0.5 - 0.6',
      '0.4 - 0.5',
      '0.3 - 0.4',
      '0.2 - 0.3',
      '0.1 - 0.2',
      '0.0 - 0.1'
    ))
    )%>%
    ggplot(aes(x = p_ij, fill = h_ijminor_cat))+
    geom_histogram(binwidth = 0.01, position = 'stack')+
    facet_wrap(~factor(Type_of_polymorphism,
                       levels = c('SNP',
                                  'INDEL',
                                  'INDEL:Trinucleotide_STR',
                                  'INDEL:Tetranucleotide_STR',
                                  'INDEL:Pentanucleotide_STR',
                                  'INDEL:Hexanucleotide_STR')), scales = 'free_y')+
    theme_bw()+
    scale_fill_viridis_d(direction = -1)+
    labs(x = 'Alternative allele frequency (p_ij)',
         y = 'Number of alternative alleles')
  
  
  plot_hij_histogram = allele_count_frac_ofHet_pAlt %>%
    mutate(h_ijminor_cat = case_when(
      h_ijminor >= .9 ~ '0.9 - 1.0',
      h_ijminor < .9 & h_ijminor >= .8 ~ '0.8 - 0.9',
      h_ijminor < .8 & h_ijminor >= .7 ~ '0.7 - 0.8',
      h_ijminor < .7 & h_ijminor >= .6 ~ '0.6 - 0.7',
      h_ijminor < .6 & h_ijminor >= .5 ~ '0.5 - 0.6',
      h_ijminor < .5 & h_ijminor >= .4 ~ '0.4 - 0.5',
      h_ijminor < .4 & h_ijminor >= .3 ~ '0.3 - 0.4',
      h_ijminor < .3 & h_ijminor >= .2 ~ '0.2 - 0.3',
      h_ijminor < .2 & h_ijminor >= .1 ~ '0.1 - 0.2',
      h_ijminor < .1 ~ '0.0 - 0.1'
    ),
    h_ijminor_cat = factor(h_ijminor_cat, levels = c(
      '0.9 - 1.0',
      '0.8 - 0.9',
      '0.7 - 0.8',
      '0.6 - 0.7',
      '0.5 - 0.6',
      '0.4 - 0.5',
      '0.3 - 0.4',
      '0.2 - 0.3',
      '0.1 - 0.2',
      '0.0 - 0.1'
    ))
    )%>%
    ggplot(aes(x = h_ij, fill = h_ijminor_cat))+
    geom_histogram(binwidth = .01,alpha = .9, position = 'stack')+
    facet_wrap(~factor(Type_of_polymorphism,
                       levels = c('SNP',
                                  'INDEL',
                                  'INDEL:Trinucleotide_STR',
                                  'INDEL:Tetranucleotide_STR',
                                  'INDEL:Pentanucleotide_STR',
                                  'INDEL:Hexanucleotide_STR')), scales = 'free_y')+
    theme_bw()+
    scale_fill_viridis_d(direction = -1)+
    labs(x = 'h_ij (H_ij/P_ij)',
         y = 'Number of alternative alleles')
  
  
  
  ### Mask rGenome_object_filtered----
  rGenome_object_filtered_masked = rGenome_object_filtered
  
  gt_masked = NULL
  
  for(w in 1:nchunks){
    gt_masked = rbind(gt_masked,
                      mask_alt_alleles(
                        rGenome_object_filtered,
                        w = w,
                        n = nchunks,
                        mask_formula = mask_formula))
    print(paste0('Chunk ', w, ' masked'))
  }
  
  rGenome_object_filtered_masked@gt = gt_masked
  
  
  #### Update allele counts----

  
  allele_counts = NULL
  
  for(w in 1:nchunks){
    allele_counts = rbind(allele_counts, get_AC(obj = rGenome_object_filtered_masked, w = w, n = nchunks))
    print(w)
    rm(w)
  }
  
  rGenome_object_filtered_masked@loci_table$Cardinality = NULL
  rGenome_object_filtered_masked@loci_table$Alleles = NULL
  rGenome_object_filtered_masked@loci_table$Allele_Counts = NULL
  
  rGenome_object_filtered_masked@loci_table = cbind(rGenome_object_filtered_masked@loci_table,
                                    allele_counts)
  rm('allele_counts')
  
  print('update allele counts and allele frequencies done')
  
  #### Filter monomorphic sites respect to the reference----
  
  removed_monomorphic_sites_rGenome_filtered_masked = filter_loci(rGenome_object_filtered_masked,
                                                  v = (rGenome_object_filtered_masked@loci_table$Cardinality <= 1 &
                                                         grepl('^0',rGenome_object_filtered_masked@loci_table$Allele_Counts)))
  
  rGenome_object_filtered_masked = filter_loci(rGenome_object_filtered_masked,
                               v = !(rGenome_object_filtered_masked@loci_table$Cardinality <= 1 &
                                       grepl('^0',rGenome_object_filtered_masked@loci_table$Allele_Counts)))
  
  print('removal of no polymorphic sites done')
  
  #### Update alternative alleles----
  
  rGenome_object_filtered_masked@loci_table$ALT =
    ifelse(rGenome_object_filtered_masked@loci_table$REF !=
             gsub(',([ATCG]|\\*)+',
                  '',
                  gsub(':\\d+',
                       '',
                       rGenome_object_filtered_masked@loci_table$Alleles)),
           gsub(':\\d+',
                '',
                rGenome_object_filtered_masked@loci_table$Alleles),
           gsub('^([ATCG]|\\*)+,',
                '',
                gsub(':\\d+', '', rGenome_object_filtered_masked@loci_table$Alleles))
    )
  
  print('update list of alternative alleles')
  
  #### Update observed heterozygosity and Heterozygous samples per Alternative alleles per site----
  
  #### Proportion of Heterozygous samples per site
  
  ObsHet = NULL
  
  for(w in 1:nchunks){
    
    ObsHet = c(ObsHet, get_ObsHet(rGenome_object_filtered_masked, by = 'loci', w = w, n = nchunks))
    rm(w)
  }
  rGenome_object_filtered_masked@loci_table$ObsHet = ObsHet
  rm(ObsHet)
  
  #### Fraction of Heterozygous samples per Alternative alleles per site
  
  frac_ofHet_pAlts = NULL
  
  for(w in 1:nchunks){
    frac_ofHet_pAlts = c(frac_ofHet_pAlts, frac_ofHet_pAlt(rGenome_object_filtered_masked, w = w, n = nchunks))
    rm(w)
  }
  
  rGenome_object_filtered_masked@loci_table$frac_ofHet_pAlts = frac_ofHet_pAlts
  rm(frac_ofHet_pAlts)
  
  #### Classify sites based on the proportion of alternative alleles in polyclonal samples per site
  rGenome_object_filtered_masked@loci_table %<>% mutate(ALT_FILTER =case_when(
    frac_ofHet_pAlts == 1 ~ '100%',
    frac_ofHet_pAlts < 1 & frac_ofHet_pAlts > .5 ~ '50 - 99%',
    frac_ofHet_pAlts <= .5 ~ '<=50%',
  ))
  
  ### Re-classify SNPs, INDELS, Homopolymers and Short tandem repeats----
  
  Type_of_polymorphism = NULL
  
  for(w in 1:nchunks){
    Type_of_polymorphism = c(Type_of_polymorphism, get_type_of_polymorphism(rGenome_object_filtered_masked, w = w, n = nchunks))
    rm(w)
  }
  
  rGenome_object_filtered_masked@loci_table$Type_of_polymorphism = Type_of_polymorphism
  rm(Type_of_polymorphism)

  
  ## Export vcf file and rGenome----
  filtered_vcf = rGenome2vcf(rGenome_object_filtered)
  write.table(filtered_vcf, final_vcf_name1, quote = FALSE, row.names = FALSE, sep = '\t')
  
  print('the filtered vcf file has been exported')
  
  assign(paste0(rGenome_object_name, '_filtered'), rGenome_object_filtered)
  
  
  
  filtered_masked_vcf = rGenome2vcf(rGenome_object_filtered_masked)
  write.table(filtered_masked_vcf, final_vcf_name2, quote = FALSE, row.names = FALSE, sep = '\t')

  print('the filtered and masked vcf file has been exported')
  
  assign(paste0(rGenome_object_name, '_filtered_masked'), rGenome_object_filtered_masked)
  
  
  save(file = gsub('.RData','_filtered.RData', imagename), list = ls()[
    grep(
      paste(
        paste0(rGenome_object_name, '_filtered'),
        paste0(rGenome_object_name, '_filtered_masked'),
        'LocusAmpRate_Summ',
        'SampAmpRate_Summ',
        'Read_Depth_Summ',
        'plot_amplifiedsamples_amplifiedloci_curve',
        'plot_frac_ofHet_pAlt_histogram',
        'plot_obHet_histogram',
        'manhatan_plot_ObsHet_by_AltFilter',
        'plot_mean_obHet_histogram',
        'ObsHet_Threshold',
        'manhatan_plot_ObsHet_by_ObsHetFilter',
        'plot_SNPdensity_histogram',
        'SNP_density_threshold',
        'manhatan_plot_ObsHet_by_SNPdensityFilter',
        'manhatan_plot_ObsHet_by_ObsHetSNPdensityFilter',
        'removed_samples_rGenome',
        'removed_monomorphic_sites_rGenome',
        'removed_sites_LowAmplRate_rGenome',
        'removed_polymorphism_rGenome',
        'removed_genomicregions',
        'allele_count_frac_ofHet_pAlt',
        'plot_AlleleCount_HetCount',
        'plot_pij_hij_hijminor',
        'plot_AlleleFreq_histogram',
        'plot_hij_histogram',
        sep = '|'),ls())])
  
  print('RData image with plots of the post-filtering saved')
  
  ## Generate traditional filtered VCF----
  
  write.table(rGenome_object_filtered@loci_table[,c('CHROM', 'POS')],
              'selected_pos.bed',
              quote = F,
              row.names = F,
              sep = '\t'
              )
  
  run_vcftools(vcf = vcf_file,
               bash_file = 'run_final_filtering.sh',
               out = paste0(output, '_filtered_final2'),
               positions = 'selected_pos.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  print('Filtered traditional vcf generated')
  
  ## Generate traditional filtered masked VCF----
  
  write.table(rGenome_object_filtered_masked@loci_table[,c('CHROM', 'POS')],
              'selected_pos.bed',
              quote = F,
              row.names = F,
              sep = '\t'
  )
  
  run_vcftools(vcf = vcf_file,
               bash_file = 'run_final_filtering_masking.sh',
               out = paste0(output, '_filtered_masked_final2'),
               positions = 'selected_pos.bed',
               recode = TRUE,
               recode_INFO_all = TRUE)
  
  rm(list = c('rGenome_object', 'rGenome_object_filtered', 'rGenome_object_filtered_masked'))
  
  print('Filtered and masked traditional vcf generated')
  
}









