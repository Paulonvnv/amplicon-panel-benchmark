# run_vcftools ----
run_vcftools = function(vcf = NULL,
                        
                        bash_file = NULL, 
                        
                        out = NULL,
                        
                        keep_regexp = NULL, # regular expression pattern that identify samples
                        remove_regexp = NULL,
                        
                        keep = NULL,
                        remove = NULL,
                        
                        chr = NULL,
                        not_chr = NULL,
                        
                        bed = NULL,
                        exclude_bed = NULL,
                        
                        positions = NULL,
                        exclude_positions = NULL,
                        
                        keep_only_indels = FALSE,
                        remove_indels = FALSE,
                        
                        remove_filtered_all = FALSE,
                        
                        maf = NULL,
                        max_maf = NULL,
                        
                        non_ref_af = NULL,
                        max_non_ref_af = NULL,
                        non_ref_ac = NULL,
                        max_non_ref_ac = NULL,
                        
                        non_ref_af_any = NULL,
                        max_non_ref_af_any = NULL,
                        non_ref_ac_any = NULL,
                        max_non_ref_ac_any = NULL,
                        
                        mac = NULL,
                        max_mac = NULL,
                        
                        min_alleles = NULL,
                        max_alleles = NULL,
                        
                        # Output options
                        freq = FALSE,
                        counts = FALSE,
                        
                        depth = FALSE,
                        site_depth = FALSE,
                        site_mean_depth = FALSE,
                        geno_depth = FALSE,
                        
                        hap_r2 = FALSE,
                        geno_r2 = FALSE,
                        geno_chisq = FALSE,
                        hap_r2_positions = NULL,
                        geno_r2_positions = NULL,
                        ld_window = NULL,
                        ld_window_bp = NULL,
                        ld_window_min = NULL,
                        ld_window_bp_min = NULL,
                        min_r2 = NULL,
                        interchrom_hap_r2 = FALSE,
                        interchrom_geno_r2 = FALSE,
                        
                        TsTv = NULL,
                        TsTv_by_count = FALSE,
                        TsTv_by_qual = FALSE,
                        
                        site_pi = FALSE,
                        window_pi = NULL,
                        window_pi_step = NULL,
                        
                        weir_fst_pop = NULL,
                        fst_window_size = NULL,
                        fst_window_step = NULL,
                        
                        het = FALSE,
                        TajimaD = NULL,
                        relatedness = FALSE,
                        relatedness2 = FALSE,
                        
                        recode = FALSE,
                        recode_bcf = FALSE,
                        recode_INFO_all = FALSE
){
  
  print('Starting VCFTools')
  
  vcf_run_file = c('#!/bin/bash',
                   'source /broad/software/scripts/useuse',
                   'use .vcftools-0.1.14',
                   'use Tabix')
  
  vcf_arguments = 'vcftools'
  
  if(grepl('.vcf$', vcf)){
    vcf_arguments = paste0(vcf_arguments, ' --vcf ', vcf)
  }else if(grepl('.vcf.gz$', vcf)){
    vcf_arguments = paste0(vcf_arguments, ' --gzvcf ', vcf)
  }else if(grepl('.bcf$', vcf)){
      vcf_arguments = paste0(vcf_arguments, ' --bcf ', vcf)}
  
  # Filters
  
  if(!is.null(keep_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  vcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'samples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --keep samples.indv')
  }
  
  if(!is.null(remove_regexp)){
    system(paste0("zgrep '^#[A-Z]' ",  vcf, " > ", "temp_samples.indv"))
    temp_samples = as.character(read.csv("temp_samples.indv", header = FALSE, sep = '\t'))
    samples = temp_samples[grepl(keep_regexp,temp_samples)]
    write.table(samples, 'rsamples.indv', sep = '\t', quote = FALSE, row.names = FALSE, col.names = F)
    system(paste0('rm ', "temp_samples.indv"))
    vcf_arguments = paste0(vcf_arguments, ' --remove rsamples.indv')
  }
  
  if(!is.null(keep)){
    vcf_arguments = paste0(vcf_arguments, ' --keep ', keep)
  }
  
  if(!is.null(remove)){
    vcf_arguments = paste0(vcf_arguments, ' --remove ', remove)
  }
  
  if(!is.null(chr)){
    for(chromosome in chr){
      vcf_arguments = paste0(vcf_arguments, ' --chr ', chromosome)
    }
  }
  
  if(!is.null(not_chr)){
    for(chromosome in not_chr){
      vcf_arguments = paste0(vcf_arguments, ' --not-chr ', chromosome)
    }
  }
  
  if(!is.null(bed)){vcf_arguments = paste0(vcf_arguments, ' --bed ', bed)
  print('Adding bed argument')
  }
  if(!is.null(exclude_bed)){vcf_arguments = paste0(vcf_arguments, ' --exclude-bed ', exclude_bed)}
  
  if(!is.null(positions)){vcf_arguments = paste0(vcf_arguments, ' --positions ', positions)}
  if(!is.null(exclude_positions)){vcf_arguments = paste0(vcf_arguments, ' --exclude-positions ', exclude_positions)}
  
  if(keep_only_indels){vcf_arguments = paste0(vcf_arguments, ' --keep-only-indels')}
  if(remove_indels){vcf_arguments = paste0(vcf_arguments, ' --remove-indels')}
  
  if(remove_filtered_all){vcf_arguments = paste0(vcf_arguments, ' --remove-filtered-all')}
  
  if(!is.null(maf)){vcf_arguments = paste0(vcf_arguments, ' --maf ', maf)}
  if(!is.null(max_maf)){vcf_arguments = paste0(vcf_arguments, ' --max-maf ', max_maf)}
  
  if(!is.null(non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af ', non_ref_af)}
  if(!is.null(max_non_ref_af)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af ', max_non_ref_af)}
  if(!is.null(non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac ', non_ref_ac)}
  if(!is.null(max_non_ref_ac)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac ', max_non_ref_ac)}
  
  if(!is.null(non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-af-any ', non_ref_af_any)}
  if(!is.null(max_non_ref_af_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-af-any ', max_non_ref_af_any)}
  if(!is.null(non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --non-ref-ac-any ', non_ref_ac_any)}
  if(!is.null(max_non_ref_ac_any)){vcf_arguments = paste0(vcf_arguments, ' --max-non-ref-ac-any ', max_non_ref_ac_any)}
  if(!is.null(mac)){vcf_arguments = paste0(vcf_arguments, ' --mac ', mac)}
  if(!is.null(max_mac)){vcf_arguments = paste0(vcf_arguments, ' --max-mac ', max_mac)}
  
  if(!is.null(min_alleles)){vcf_arguments = paste0(vcf_arguments, ' --min-alleles ', min_alleles)}
  if(!is.null(max_alleles)){vcf_arguments = paste0(vcf_arguments, ' --max-alleles ', max_alleles)}
  
  # Output options
  
  if(freq){vcf_arguments = paste0(vcf_arguments, ' --freq')}
  if(counts){vcf_arguments = paste0(vcf_arguments, ' --counts')}
  
  if(depth){vcf_arguments = paste0(vcf_arguments, ' --depth')}
  if(site_depth){vcf_arguments = paste0(vcf_arguments, ' --site-depth')}
  if(site_mean_depth){vcf_arguments = paste0(vcf_arguments, ' --site-mean-depth')}
  if(geno_depth){vcf_arguments = paste0(vcf_arguments, ' --geno_depth')}
  
  if(hap_r2){vcf_arguments = paste0(vcf_arguments, ' --hap-r2')}
  if(geno_r2){vcf_arguments = paste0(vcf_arguments, ' --geno-r2')}
  if(geno_chisq){vcf_arguments = paste0(vcf_arguments, ' --geno-chisq')}
  if(!is.null(ld_window)){vcf_arguments = paste0(vcf_arguments, ' --ld-window ', ld_window)}
  if(!is.null(ld_window_bp)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp ', ld_window_bp)}
  if(!is.null(ld_window_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-min ', ld_window_min)}
  if(!is.null(ld_window_bp_min)){vcf_arguments = paste0(vcf_arguments, ' --ld-window-bp-min ', ld_window_bp_min)}
  if(!is.null(min_r2)){vcf_arguments = paste0(vcf_arguments, ' --min-r2 ', min_r2)}
  if(interchrom_hap_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-hap-r2')}
  if(interchrom_geno_r2){vcf_arguments = paste0(vcf_arguments, ' --interchrom-geno-r2')}
  
  if(!is.null(TsTv)){vcf_arguments = paste0(vcf_arguments, ' --TsTv ', TsTv)}
  if(TsTv_by_count){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-count')}
  if(TsTv_by_qual){vcf_arguments = paste0(vcf_arguments, ' --TsTv-by-qual')}
  
  if(site_pi){vcf_arguments = paste0(vcf_arguments, ' --site-pi')}
  if(!is.null(window_pi)){vcf_arguments = paste0(vcf_arguments, ' --window-pi ', window_pi)}
  if(!is.null(window_pi_step)){vcf_arguments = paste0(vcf_arguments, ' --window-pi-step ', window_pi_step)}
  
  if(!is.null(weir_fst_pop)){
    for(pop in weir_fst_pop){
      vcf_arguments = paste0(vcf_arguments, ' --weir-fst-pop ', pop)
    }
  }
  
  if(!is.null(fst_window_size)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-size ', fst_window_size)}
  if(!is.null(fst_window_step)){vcf_arguments = paste0(vcf_arguments, ' --fst-window-step ', fst_window_step)}
  
  if(het){vcf_arguments = paste0(vcf_arguments, ' --het')}
  if(!is.null(TajimaD)){vcf_arguments = paste0(vcf_arguments, ' --TajimaD ', TajimaD)}
  if(relatedness){vcf_arguments = paste0(vcf_arguments, ' --relatedness')}
  if(relatedness2){vcf_arguments = paste0(vcf_arguments, ' --relatedness2')}
  
  if(recode){vcf_arguments = paste0(vcf_arguments, ' --recode')}
  if(recode_bcf){vcf_arguments = paste0(vcf_arguments, ' --recode-bcf')}
  if(recode_INFO_all){vcf_arguments = paste0(vcf_arguments, ' --recode-INFO-all')}
  
  if(!is.null(out)){vcf_arguments = paste0(vcf_arguments, ' --out ', out)}
  
  vcf_run_file = c(vcf_run_file, vcf_arguments)
  
  print(vcf_arguments)
  
  write.table(vcf_run_file, bash_file, row.names = FALSE, quote = FALSE, col.names = FALSE)
  
  system(paste0('chmod 777 ', bash_file))
  system(paste0('./', bash_file))
  system(paste0('rm ', bash_file))
  
}

# load_vcf----
load_vcf = function(vcf = NULL,
                    gzvcf = NULL,
                    na.rm = TRUE,
                    start = NULL,
                    end = NULL
                    ){
  
  if(!is.null(vcf)){
    
    temp_tsv_name = gsub('vcf', 'tsv', vcf)
    system(paste0("grep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  if(!is.null(gzvcf)){
    system(paste0("zgrep -v '^##' ", vcf, " > ", temp_tsv_name))
  }
  
  col_names = strsplit(system(paste0("grep '^#' ", temp_tsv_name), intern = TRUE), '\t')[[1]]
  col_names[1] = 'CHROM'
  
  temp_tsv_name2 = gsub('.tsv$','_2.tsv',temp_tsv_name)
  
  if(is.null(start) | is.null(end)){
    system(paste0("grep -v '^#' ", temp_tsv_name, " > ", temp_tsv_name2))
  }else if(!is.null(start) & !is.null(end)){
    system(paste0("grep -v '^#' ", 
                  temp_tsv_name, 
                  " | sed -n ", 
                  as.character(as.integer(start)), 
                  ",", 
                  as.character(as.integer(end)), 
                  "p > ", 
                  temp_tsv_name2))
  }
  
  vcf = read.table(temp_tsv_name2, header = FALSE)
  names(vcf) = col_names
  
  if(na.rm){
    vcf %<>% filter(ALT != '.')
  }
  
  system(paste0('rm ', temp_tsv_name))
  system(paste0('rm ', temp_tsv_name2))
  
  return(vcf)
}


# rGenome S4class and vcf2rGenome----

## rGenome S4 class

setClass('rGenome', slots = c(
  gt = "ANY",
  loci_table = "ANY",
  metadata = "ANY"
))

## rGenome constructor
rGenome = function(gt = NULL,
                   loci_table = NULL,
                   metadata = NULL){
  obj = new('rGenome')
  obj@gt = gt
  obj@loci_table = loci_table
  obj@metadata = metadata
  
  return(obj)
}

## vcf2rGenome----

vcf2rGenome = function(vcf, n = 500, threshold = 5) {
            
            # Generate metadata
            metadata = data.frame(Sample_id = names(vcf)[-1:-9])
            rownames(metadata) = metadata[['Sample_id']]
            
            # generate loci_table
            loci_table = vcf[,1:9]
            rownames(loci_table) = paste(loci_table$CHROM, loci_table$POS, sep = '_')
            
            # generate a haplotype table (gt)
            
            gt = NULL
            for(w in 1:n){
              start = Sys.time()
              gt = rbind(gt, get_GTAD_matrix(vcf, w = w, n = n, threshold = threshold))
              end = Sys.time()
              print(w)
              print(end-start)
            }
            
            obj = rGenome(gt = gt, loci_table = loci_table, metadata = metadata)
          
            return(obj)
          }

## rGenome2vcf----
rGenome2vcf = function(rGenome_object){
  
  loci_table = rGenome_object@loci_table
  
  vcf_names1 = c("CHROM", "POS", "ID",  "REF", "ALT",  "QUAL", "FILTER", "INFO", "FORMAT")
  
  loci_table$INFO = paste(loci_table$INFO, apply(loci_table[, !(colnames(loci_table) %in% vcf_names1)], 1, function(x){paste(paste(names(x), x, sep = '='), collapse = ';')}), sep = ';')
  
  loci_table$FORMAT = 'GT:AD'

  loci_table = loci_table[,vcf_names1]
  
  colnames(loci_table) = c("#CHROM", "POS", "ID",  "REF", "ALT",  "QUAL", "FILTER", "INFO", "FORMAT")
  
  gt = rGenome_object@gt

  ad = matrix(gsub('\\d+:', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  ad = gsub('/', ',', ad)
  
  ad[is.na(ad)] = 0
  
  gt = matrix(gsub(':\\d+', '', gt), nrow = nrow(gt), ncol = ncol(gt),
                dimnames = list(rownames(gt), colnames(gt)))
  
  
  gt1 = gsub('/\\d+', '', gt)
  gt2 = gsub('\\d+/', '', gt)
  
  gt3 = matrix(paste(gt1, gt2 , sep = '/'), nrow = nrow(gt), ncol = ncol(gt),
               dimnames = list(rownames(gt), colnames(gt)))
  
  gt3 = gsub('NA', '.', gt3)
  
  final_gt = matrix(paste(gt3, ad , sep = ':'), nrow = nrow(gt), ncol = ncol(gt),
                    dimnames = list(rownames(gt), colnames(gt)))
  
  final_vcf = cbind(loci_table, final_gt)
  
  return(final_vcf)
  
}

# SampleAmplRate----

setGeneric("SampleAmplRate", function(obj, update = TRUE, type = 'freq', threshold = NA, n = 100) standardGeneric("SampleAmplRate"))

setMethod("SampleAmplRate", signature(obj = "rGenome"),
          
          function(obj, update = TRUE, type = 'freq', threshold = NA, n = 100) {
            
            obj2 = obj
            
            if(!is.na(threshold)){
              obj2@gt = prune_alleles(obj = obj2, threshold = threshold, n = n)
            }
            
            if(update){
              
              if(type == 'freq'){
                
                obj2@metadata[['SampleAmplRate']] =
                  colSums(!is.na(obj2@gt))/nrow(obj2@gt)
                
              }else if(type == 'count'){
                
                obj2@metadata[['NumberOfAmplifiedLoci']] =
                  colSums(!is.na(obj2@gt))
                
              }
              
              return(obj2)
              
            }else{
              
              if(type == 'freq'){
                
                result = colSums(!is.na(obj2@gt))/nrow(obj2@gt)
                
              }else if(type == 'count'){
                
                result = colSums(!is.na(obj2@gt))
                
              }
              
              return(result)
            }
            
          }
)

# LocusAmplRate----

setGeneric("LocusAmplRate", function(obj, update = TRUE, by = NA, threshold = NA, n = 100) standardGeneric("LocusAmplRate"))

setMethod("LocusAmplRate", signature(obj = "rGenome"),
          
          function(obj, update = TRUE, by = NA, threshold = NA, n = 100) {
            
            obj2 = obj
            
            if(!is.na(threshold)){
              obj2@gt = prune_alleles(obj = obj2, threshold = threshold, n = n)
            }
            
            if(update & is.na(by)){
              obj2@loci_table[['LocusAmplRate']] =
                1 - rowSums(is.na(obj2@gt))/ncol(obj2@gt)
              
              return(obj2)
              
            }else{
              
              if(!is.na(by)){
                
                result = NULL
                
                for(Pop in unique(obj2@metadata[,by])){
                  
                  temp_object = filter_samples(obj2, v = obj2@metadata[,by] == Pop)
                  temp_result = 1 - rowSums(is.na(temp_object@gt))/ncol(temp_object@gt)
                  
                  result = cbind(result, temp_result)
                  
                }
                
                colnames(result) = unique(obj2@metadata[,by])
                rownames(result) = rownames(obj2@gt)
              
              }else{
                
                result = 1 - rowSums(is.na(obj2@gt))/ncol(obj2@gt)
                
              }
              
              return(result)
            }
            
          }
)

# # filter_samples----
# 
# setGeneric("filter_samples", function(obj, v = NULL) standardGeneric("filter_samples"))
# 
# setMethod("filter_samples", signature(obj = "rGenome"),
#           
#           function(obj, v = NULL) {
#             
#             obj2 = rGenome(gt = obj@gt[,v],
#                            loci_table = obj@loci_table,
#                            metadata = obj@metadata[v,])
#             
#             return(obj2)
#           }
# )
# 
# # filter_loci----
# 
# setGeneric("filter_loci", function(obj, v = NULL) standardGeneric("filter_loci"))
# 
# setMethod("filter_loci", signature(obj = "rGenome"),
#           
#           function(obj, v = NULL) {
#             
#             obj2 = obj
#             obj2@gt = obj2@gt[v,]
#             obj2@loci_table = obj2@loci_table[v,]
#             
#             if(is.null(nrow(obj2@gt)) & ncol(obj@gt) > 1){
#               
#               obj2@gt = matrix(obj2@gt, nrow = 1, ncol = length(obj2@gt),
#                                dimnames = list(
#                                  rownames(obj2@loci_table),
#                                  names(obj2@gt)
#                                ))
#               
#             }else if(is.null(nrow(obj2@gt)) & ncol(obj@gt) == 1){
#               
#               obj2@gt = matrix(obj2@gt, ncol = 1, nrow = length(obj2@gt),
#                                dimnames = list(
#                                  rownames(obj2@loci_table),
#                                  colnames(obj@gt)
#                                ))
#               
#             }
#             
#             return(obj2)
#           }
# )

# handle_ploidy----

handle_ploidy = function(gt, monoclonals, polyclonals, w = 1, n = 1){
  
  s = round(seq(1,nrow(gt)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  if(sum(grepl(':', gt)) > 0){
      gt = matrix(gsub(':\\d+', '', gt[low:high,]), nrow = high - low + 1, ncol = ncol(gt),
                  dimnames = list(rownames(gt)[low:high], colnames(gt)))
  }
  
  # if(sum(grepl('/', gt)) == 0){
  #   
  #   gt3 = gt
  #   
  # }else
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    gt1 = gsub('/\\d+', '', gt)
    gt2 = gsub('\\d+/', '', gt)
    gt2[gt2==gt1] = NA
    gt3 = cbind(gt1, gt2)
  }else{
    gt_mono = matrix(gt[,monoclonals],
                     nrow = nrow(gt),
                     ncol = length(monoclonals),
                     dimnames = list(rownames(gt), monoclonals))
    
    gt_mono = gsub('/\\d+', '', gt_mono)
    
    
    gt_poly = matrix(gt[,polyclonals],
                     nrow = nrow(gt),
                     ncol = length(polyclonals),
                     dimnames = list(rownames(gt), polyclonals)
                     )
    gt_poly1 = gsub('/\\d+', '', gt_poly)
    
    gt_poly2 = gsub('\\d+/', '', gt_poly)
    
    if(!is.null(polyclonals)){
      colnames(gt_poly1) = paste(colnames(gt_poly1), 'C1', sep = '_')
      colnames(gt_poly2) = paste(colnames(gt_poly2), 'C2', sep = '_')
    }
    
    gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
  }
  
  return(gt3)
  
}

# get_AC----

setGeneric("get_AC", function(obj = NULL, w = 1, n = 1,
                              update_alleles = TRUE,
                              monoclonals = NULL, 
                              polyclonals = NULL
                              ) standardGeneric("get_AC"))

setMethod("get_AC", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1,
                  update_alleles = TRUE,
                  monoclonals = NULL, 
                  polyclonals = NULL){
  
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            gt3 = handle_ploidy(gt = gt, w = w, n = n, monoclonals = monoclonals, polyclonals = polyclonals)
            
            if(update_alleles){
              if('Alleles' %in% colnames(loci)){
                loci2 = loci[low:high,'Alleles']
                loci2 = matrix(loci2, ncol = 1,
                               dimnames = list(
                                 rownames(loci[low:high,]),
                                 'Alleles'
                               ))
              }else{
                loci2 = loci[low:high,c('REF', 'ALT')]
                
                loci2[['Alleles']] = apply(loci2, 1, function(variant_site){
                  alleles = c(variant_site['REF'], strsplit(variant_site['ALT'], ',')[[1]])
                  paste(paste(alleles, 0:(length(alleles) - 1), sep = ':'), collapse = ',')
                })
                
              }
            }
            
            
            
            alleles = t(sapply(1:nrow(gt3), function(locus){
              alleles = unique(gt3[locus,])
              alleles = sort(alleles[!is.na(alleles)])
              
              if(update_alleles){
                original_alleles = gsub(':\\d+', '', strsplit(loci2[locus,'Alleles'], ',')[[1]])
                names(original_alleles) = gsub('.+:', '', strsplit(loci2[locus,'Alleles'], ',')[[1]])
              }
              
              
              
              nalleles = length(alleles)
              AC = sapply(alleles, function(allele){
                sum(gt3[locus,] == allele, na.rm = T)
              })
              
              AC = paste(paste(alleles, AC, sep = ':'), collapse = ',')
              
              if(update_alleles){
                alleles = paste(
                  paste(
                    original_alleles[alleles], alleles, sep = ':'), collapse = ',')}
              
              if(update_alleles){
                c(nalleles, alleles, AC)}else{
                  c(nalleles, AC)
                }
            }))
            
            alleles = as.data.frame(alleles)
            
            if(update_alleles){
              colnames(alleles) = c('Cardinality', 'Alleles', 'Allele_Counts')
            }else{
              colnames(alleles) = c('Cardinality', 'Allele_Counts')
            }
            rownames(alleles) = rownames(gt3)
            alleles$Cardinality = as.integer(alleles$Cardinality)
            
            
            return(alleles)
            
          }
)


# get_ExpHet----

setGeneric("get_ExpHet", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL) standardGeneric("get_ExpHet"))

setMethod("get_ExpHet", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL){
  
  gt = obj@gt
  loci = obj@loci_table
  metadata = obj@metadata
  
  if(!is.null(by)){
    
    populations = t(table(metadata[[by]]))
    populations = data.frame(population = colnames(populations), nsamples = populations[1,])
    
    ExpHet = NULL
    
    for(pop in populations$population){
      
      if(populations[pop,][['nsamples']] >= 2){
        
        samples = metadata[metadata[[by]] == pop,][['Sample_id']]
        temp_pop = filter_samples(obj = obj, v = samples)
  
        temp_monoclonals = monoclonals[monoclonals %in% samples]
        if(length(temp_monoclonals) == 0){
          temp_monoclonals = NULL
        }
        temp_polyclonals = polyclonals[polyclonals %in% samples]
        if(length(temp_polyclonals) == 0){
          temp_polyclonals = NULL
        }
        
        temp_AC = get_AC(obj = temp_pop, w =1, n = 1, update_alleles = FALSE, monoclonals = temp_monoclonals, polyclonals = temp_polyclonals)
        
        temp_alleles_counts = temp_AC$Allele_Counts
        
        temp_AC = NULL
        for(x in 1:length(temp_alleles_counts)){
          alleles_count = strsplit(temp_alleles_counts[x], ',')[[1]]
          temp_AC[[x]] = gsub('^\\d+:', '', alleles_count)
        }
        
        temp_ExpHet = sapply(1:length(temp_AC), function(pos){
          n = sum(as.numeric(temp_AC[[pos]]))
          allele_counts = as.numeric(temp_AC[[pos]])
          allele_freq = allele_counts/n
          sp2 = sum(allele_freq^2)
          ExpHet = n*(1 - sp2)/(n-1)
        })
        
      }else{
        
        # fill with NA's if there is only one sample in the population
        temp_ExpHet = rep(NA, nrow(temp_pop@gt))
        
      }
      
      ExpHet = cbind(ExpHet, temp_ExpHet)
      
    }
    
    # Calculate Heterozygosity for the whole population
    
    AC = get_AC(obj = obj, w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
    alleles_counts = AC$Allele_Counts
    
    AC = NULL
    for(x in 1:length(alleles_counts)){
      alleles_count = strsplit(alleles_counts[x], ',')[[1]]
      AC[[x]] = gsub('^\\d+:', '', alleles_count)
    }
    
    Total = unlist(sapply(1:length(AC), function(pos){
      n = sum(as.numeric(AC[[pos]]))
      allele_counts = as.numeric(AC[[pos]])
      allele_freq = allele_counts/n
      sp2 = sum(allele_freq^2)
      ExpHet = n*(1 - sp2)/(n-1)
    }))
    
    ExpHet = cbind(ExpHet, Total)
    
    colnames(ExpHet) = c(populations$population, 'Total')
    rownames(ExpHet) = rownames(gt)
    
    ExpHet = cbind(loci[,c('CHROM','POS','Cardinality','Type_of_polymorphism')], ExpHet)
    
  }else{
    
    if(!update_AC){
      AC = sapply(1:nrow(loci), function(x){
        AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
        gsub('^\\d+:', '', AC)
      })
    }else{
      
      AC = get_AC(obj = obj, w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
      alleles_counts = AC$Allele_Counts
      
      AC = NULL
      for(x in 1:length(alleles_counts)){
        alleles_count = strsplit(alleles_counts[x], ',')[[1]]
        AC[[x]] = gsub('^\\d+:', '', alleles_count)
      }
    }
    
    ExpHet = base::sapply(1:length(AC), function(pos){
      n = sum(as.numeric(AC[[pos]]))
      allele_counts = as.numeric(AC[[pos]])
      allele_freq = allele_counts/n
      sp2 = sum(allele_freq^2)
      ExpHet = n*(1 - sp2)/(n-1)
    })
    
  }
  
  return(ExpHet)
})

# get_private_alleles----

setGeneric("get_private_alleles", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL, min_samp = 1, fixed = FALSE, tolerance = 0.25) standardGeneric("get_private_alleles"))

setMethod("get_private_alleles", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL, by = NULL, min_samp = 1, fixed = FALSE, tolerance = 0.25){
            
            if(is.null(min_samp)|is.na(min_samp)){
              min_samp = 1
            }
            
            gt = obj@gt
            loci = obj@loci_table
            metadata = obj@metadata
            
            if(!is.null(by)){
              
              populations = t(table(metadata[[by]]))
              populations = data.frame(population = colnames(populations), nsamples = populations[1,])
              
              Pop_alleles = matrix(NA,
                                   nrow = nrow(loci),
                                   ncol = nrow(populations),
                                   dimnames = list(rownames(loci),
                                                   populations$population))
              
              for(pop in populations$population){
                
                samples = metadata[metadata[[by]] == pop,][['Sample_id']]
                temp_pop = filter_samples(obj = obj, v = samples)
                
                temp_monoclonals = monoclonals[monoclonals %in% samples]
                if(length(temp_monoclonals) == 0){
                  temp_monoclonals = NULL
                }
                temp_polyclonals = polyclonals[polyclonals %in% samples]
                if(length(temp_polyclonals) == 0){
                  temp_polyclonals = NULL
                }
                
                temp_AC = get_AC(obj = temp_pop, w =1, n = 1, update_alleles = FALSE, monoclonals = temp_monoclonals, polyclonals = temp_polyclonals)
                
                
                Pop_alleles[,pop] = temp_AC$Allele_Counts
                
              }
              
              Pop_alleles[Pop_alleles == ''] = NA
              
              if(fixed){
                
                maxn_amplified_samples = gsub('^.+:', '', gsub(',.+$', '', Pop_alleles))
                
                maxn_amplified_samples = apply(maxn_amplified_samples, 2, function(pop){
                  max(as.integer(pop), na.rm = T)
                })
                
                if(tolerance < 1){
                  
                  maxn_amplified_samples = ceiling(maxn_amplified_samples*(1-tolerance))
                  
                }else{
                  maxn_amplified_samples = maxn_amplified_samples - tolerance
                }
                
              }
              
              private_alleles = matrix(NA,
                                   nrow = nrow(loci),
                                   ncol = nrow(populations),
                                   dimnames = list(rownames(loci),
                                                   populations$population))
              
              for(variant_site in 1:nrow(private_alleles)){
                for(pop in colnames(private_alleles)){
                  
                  tested_pop_alleles = unlist(strsplit(gsub(':\\d+','',Pop_alleles[variant_site, pop]), ','))
                  
                  tested_pop_alleles_counts = as.integer(unlist(strsplit(gsub('\\d+:','',Pop_alleles[variant_site, pop]), ',')))
                  
                  if(fixed){
                    
                    tested_pop_alleles = tested_pop_alleles[tested_pop_alleles_counts >= maxn_amplified_samples[pop]]
                    
                  }else{
                    
                    tested_pop_alleles = tested_pop_alleles[tested_pop_alleles_counts >= min_samp]
                    
                  }
                  
                  
                  
                  tested_pop_alleles = tested_pop_alleles[!is.na(tested_pop_alleles)]
                  
                  if(length(tested_pop_alleles) == 1){
                    
                    if(!is.na(tested_pop_alleles) & !is.null(tested_pop_alleles)){
                      
                      overall_pop_alleles = unique(unlist(strsplit(gsub(':\\d+','',Pop_alleles[variant_site, !(colnames(Pop_alleles) %in% pop)]), ',')))
                      
                      if(sum(!(tested_pop_alleles %in% overall_pop_alleles)) > 0){
                        
                        private_alleles[variant_site, pop] = paste(tested_pop_alleles[!(tested_pop_alleles %in% overall_pop_alleles)], collapse = ',')
                        
                      }
                      
                    }
                    
                  }else if(length(tested_pop_alleles) > 1){
                    
                    overall_pop_alleles = unique(unlist(strsplit(gsub(':\\d+','',Pop_alleles[variant_site, !(colnames(Pop_alleles) %in% pop)]), ',')))
                    
                    if(sum(!(tested_pop_alleles %in% overall_pop_alleles)) > 0){
                      
                      private_alleles[variant_site, pop] = paste(tested_pop_alleles[!(tested_pop_alleles %in% overall_pop_alleles)], collapse = ',')
                      
                    }
                    
                  }
                }
              }
              
                
              
            }else{
              
              stop("Categorical Variable should be provided using 'by = ' argument")
              
            }
            
            return(list(private_alleles = private_alleles,
                        Pop_alleles = Pop_alleles))
          })


# get_foreign_alleles_bySample----


get_foreign_alleles_bySample = function(obj, 
                                          list_of_private_alleles,
                                          filter_variant_sites){
  
  list_of_private_alleles = gsub(',','(,|$)|', list_of_private_alleles)
  
  list_of_private_alleles = list_of_private_alleles[filter_variant_sites]
  
  gt = obj@gt
  
  gt = gsub(':\\d+', '', gt)
  
  gt = gt[filter_variant_sites,]
  
  foreign_alleles = NULL
  
  for(sample in colnames(gt)){
    
    foreign_alleles_by_sample = sapply(1:length(list_of_private_alleles), function(variant_site){
      grepl(list_of_private_alleles[variant_site], gt[variant_site,sample])
    })
    
    foreign_alleles = rbind(foreign_alleles,
                              data.frame(Sample_id = sample,
                                         nforeign_alleles = sum(foreign_alleles_by_sample, na.rm = T),
                                         nscreened_variant_sites = sum(!is.na(gt[,sample])),
                                         tforeign_variant_sites = sum(!is.na(list_of_private_alleles)),
                                         nforeign_variant_sites_in_sample = sum(!is.na(gt[,sample]) & !is.na(list_of_private_alleles)),
                                         frac_foreign_alleles = sum(foreign_alleles_by_sample, na.rm = T)/sum(!is.na(gt[,sample]) & !is.na(list_of_private_alleles)),
                                         variant_sites_w_foreign_alleles = paste(rownames(gt)[which(foreign_alleles_by_sample)], collapse = ',')
                              ))
    
  }
  
  return(foreign_alleles)
  
}









# get_nuc_div----

setGeneric("get_nuc_div", function(obj = NULL, monoclonals = NULL, polyclonals = NULL, gff = NULL, type_of_region = NULL, window = NULL, by = NULL, min_samp_size = 2) standardGeneric("get_nuc_div"))

setMethod("get_nuc_div", signature(obj = "rGenome"),
          function(obj = NULL, monoclonals = NULL, polyclonals = NULL, gff = NULL, type_of_region = NULL,  window = NULL, by = NULL, min_samp_size = 2){
            
            loci = obj@loci_table
            metadata = obj@metadata
            
            # Define DNA regions to calculate their nucleotide diversity
            # It could be by gene or by fixed size windows
            
            # By gene region specified by a gff object
            if(is.object(gff) & is.null(window)){
              
              dna_regions = gff
              
              # By gene region specified by a gff file
            }else if(file.exists(gff) & is.null(window)){
              
              ref_gff = ape::read.gff(gff)
              dna_regions = ref_gff[grepl(type_of_region, ref_gff$type)&
                                      !grepl('^Transfer',ref_gff$seqid),
                                    c('seqid', 'start', 'end', 'attributes')]
              
              dna_regions = dna_regions[order(dna_regions$start),]
              dna_regions = dna_regions[order(dna_regions$seqid),]
              rownames(dna_regions) = 1:nrow(dna_regions)
              
              dna_regions = cbind(dna_regions, as.data.frame(t(sapply(1:nrow(dna_regions), function(gene){
                attributes = strsplit(dna_regions[gene,][['attributes']], ';')[[1]]
                c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
                  gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
              }))))
              
              dna_regions = dna_regions[,c('seqid', 'start', 'end', 'gene_id', 'gene_description')]
              
              # By fixed size window
            }else if(!is.null(window)){
              
              ref_gff = ape::read.gff(gff)
              ref_gff = ref_gff[grepl(type_of_region, ref_gff$type)&
                                  !grepl('^Transfer',ref_gff$seqid),
                                c('seqid', 'start', 'end', 'attributes')]
              
              ref_gff = ref_gff[order(ref_gff$start),]
              ref_gff = ref_gff[order(ref_gff$seqid),]
              rownames(ref_gff) = 1:nrow(ref_gff)
              
              ref_gff = cbind(ref_gff, as.data.frame(t(sapply(1:nrow(ref_gff), function(gene){
                attributes = strsplit(ref_gff[gene,][['attributes']], ';')[[1]]
                c(gene_id = gsub('^ID=','',attributes[grep('^ID=', attributes)]),
                  gene_description = gsub('^description=','',attributes[grep('^description=', attributes)]))
              }))))
              
              chrom_length = loci %>% group_by(CHROM) %>% summarise(length = max(POS))
              
              chrom_intervals = sapply(chrom_length$length, function(chrom){
                seq(1, chrom, window)
              })
              
              
              dna_regions = NULL
              
              for(chrom in 1:length(chrom_intervals)){
                dna_regions = rbind(dna_regions, data.frame(seqid = chrom_length[chrom,][['CHROM']],
                                                            start = chrom_intervals[[chrom]],
                                                            end = chrom_intervals[[chrom]] - 1 + window))
                
              }
              
              
              dna_regions$gene_ids = sapply(1:nrow(dna_regions),function(bin){
                paste(ref_gff[ref_gff[['seqid']] == dna_regions[bin,][['seqid']] &
                                ((ref_gff[['start']] > dna_regions[bin,][['start']] &
                                    ref_gff[['start']] < dna_regions[bin,][['end']])|
                                   
                                   (ref_gff[['end']] > dna_regions[bin,][['start']] &
                                      ref_gff[['end']] < dna_regions[bin,][['end']])|
                                   
                                   (dna_regions[bin,][['start']] > ref_gff[['start']] &
                                      dna_regions[bin,][['start']] < ref_gff[['end']])|
                                   
                                   (dna_regions[bin,][['end']] > ref_gff[['start']] &
                                      dna_regions[bin,][['end']] < ref_gff[['end']])
                                ),][['gene_id']], collapse = ',')}, simplify = T)
              
              ###### Add gene_description
              dna_regions$genes_description = sapply(1:nrow(dna_regions),function(bin){
                paste(ref_gff[ref_gff[['seqid']] == dna_regions[bin,][['seqid']] &
                                ((ref_gff[['start']] > dna_regions[bin,][['start']] &
                                    ref_gff[['start']] < dna_regions[bin,][['end']])|
                                   
                                   (ref_gff[['end']] > dna_regions[bin,][['start']] &
                                      ref_gff[['end']] < dna_regions[bin,][['end']])|
                                   
                                   (dna_regions[bin,][['start']] > ref_gff[['start']] &
                                      dna_regions[bin,][['start']] < ref_gff[['end']])|
                                   
                                   (dna_regions[bin,][['end']] > ref_gff[['start']] &
                                      dna_regions[bin,][['end']] < ref_gff[['end']])
                                ),][['gene_description']], collapse = ',')}, simplify = T)
              
            }else if(is.null(gff) & is.null(window)){
              
              print('You must provide a gff file or define a window size')
              
            }
            
            
            # If calculation have to be done individually by each population
            
            if(!is.null(by)){
              
              populations = t(table(metadata[[by]]))
              populations = data.frame(population = colnames(populations), nsamples = populations[1,])
              
              # for each provided population calculates pi and pi_var
              for(pop in populations$population){
                
                # if population has at least two samples (min_samp_size can be modified)
                if(populations[pop,][['nsamples']] >= min_samp_size){
                  
                  samples = metadata[metadata[[by]] == pop,][['Sample_id']]
                  temp_pop = filter_samples(obj = obj, v = samples)
                  
                  temp_monoclonals = monoclonals[monoclonals %in% samples]
                  if(length(temp_monoclonals) == 0){
                    temp_monoclonals = NULL
                  }
                  temp_polyclonals = polyclonals[polyclonals %in% samples]
                  if(length(temp_polyclonals) == 0){
                    temp_polyclonals = NULL
                  }
                  
                  gt = temp_pop@gt
                  
                  gt3 = handle_ploidy(gt, monoclonals = temp_monoclonals, polyclonals = temp_polyclonals)
                  gt3 = as.data.frame(gt3)
                  
                  pi = NULL
                  var = NULL
                  
                  for(region in 1:nrow(dna_regions)){
                    positions = paste(dna_regions[region, ][['seqid']], dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                    
                    region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                    
                    temp_gt = gt3[rownames(gt3) %in% positions,]
                    
                    n = ncol(temp_gt)
                    
                    if(nrow(temp_gt)>0){
                      
                      temp_loci = loci[rownames(loci) %in% positions,]
                      
                      temp_indels = temp_loci[temp_loci$Type_of_polymorphism != 'SNP',]
                      
                      temp_gt = temp_gt[temp_loci$Type_of_polymorphism == 'SNP',] # Keep only SNPs
                      
                      if(nrow(temp_gt)>0){
                        
                        if(length(temp_indels$ALT) > 0){
                          
                          ALTs = temp_indels$ALT
                          
                          ALTs = strsplit(ALTs, ',')
                          
                          gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                            max(nchar(alt))
                          })
                          
                          gaps[gaps < 0] = 0
                          
                          region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                          
                        }
                        
                        haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})),
                                                     maxsum = ncol(temp_gt))
                        
                        if(length(haplotypes_counts) > 1){
                          
                          names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                          
                          haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                          haplotypes = names(haplotypes_counts)
                          
                          combinations = combn(1:length(haplotypes),2)
                          
                          temp_pi = NULL
                          
                          for(comb in 1:ncol(combinations)){
                            
                            x_i = haplotypes_freqs[combinations[1,comb]]
                            x_j = haplotypes_freqs[combinations[2,comb]]
                            
                            seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                            seq_i[seq_i == '_'] = NA
                            seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                            seq_j[seq_j == '_'] = NA
                            
                            pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                            temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                          }
                          
                          
                          pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                          var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                          
                        }else{
                          
                          pi = c(pi, 0)
                          var = c(var, 0)
                        }
                        
                      }else{
                        
                        pi = c(pi, NA)
                        var = c(var, NA)
                      }
                      
                    }else{
                      
                      pi = c(pi, NA)
                      var = c(var, NA)
                      
                    }
                    
                  }
                  
                }else{
                  
                  print(paste0(pop, ' will be excluded because sample size is less than 2 individuals')) 
                  
                }
                
                temp_dna_regions_pi = data.frame(pi, var)
                
                names(temp_dna_regions_pi) = c(paste0(pop, '_pi'),
                                               paste0(pop, '_pi_var'))
                
                dna_regions = cbind(dna_regions, temp_dna_regions_pi)
              }
              
              gt = obj@gt
              
              gt3 = handle_ploidy(gt, monoclonals = monoclonals, polyclonals = polyclonals)
              gt3 = as.data.frame(gt3)
              
              pi = NULL
              var = NULL
              
              for(region in 1:nrow(dna_regions)){
                positions = paste(dna_regions[region, ][['seqid']],dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                
                region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                
                temp_gt = gt3[rownames(gt3) %in% positions,]
                
                n = ncol(temp_gt)
                
                if(nrow(temp_gt)>0){
                  
                  temp_loci = loci[rownames(loci) %in% positions,]
                  
                  temp_indels = temp_loci[temp_loci$Type_of_polymorphism != 'SNP',]
                  
                  temp_gt = temp_gt[temp_loci$Type_of_polymorphism == 'SNP',] # Keep only SNPs
                  
                  if(nrow(temp_gt)>0){
                    
                    if(length(temp_indels$ALT) > 0){
                      
                      ALTs = temp_indels$ALT
                      
                      ALTs = strsplit(ALTs, ',')
                      
                      gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                        max(nchar(alt))
                      })
                      
                      gaps[gaps < 0] = 0
                      
                      region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                      
                    }
                    
                    haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})),
                                                 maxsum = ncol(temp_gt))
                    
                    if(length(haplotypes_counts) > 1){
                      
                      names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                      
                      haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                      haplotypes = names(haplotypes_counts)
                      
                      combinations = combn(1:length(haplotypes),2)
                      
                      temp_pi = NULL
                      
                      for(comb in 1:ncol(combinations)){
                        
                        x_i = haplotypes_freqs[combinations[1,comb]]
                        x_j = haplotypes_freqs[combinations[2,comb]]
                        
                        seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                        seq_i[seq_i == '_'] = NA
                        seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                        seq_j[seq_j == '_'] = NA
                        
                        pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                        temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                      }
                      
                      
                      pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                      var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                      
                    }else{
                      
                      pi = c(pi, 0)
                      var = c(var, 0)
                    }
                    
                  }else{
                    
                    pi = c(pi, NA)
                    var = c(var, NA)
                  }
                  
                }else{
                  
                  pi = c(pi, NA)
                  var = c(var, NA)
                  
                }
                
              }
              
              dna_regions = cbind(dna_regions, data.frame(Total_pi = pi, Total_pi_var = var))
              
              
            }else{
              
              gt = obj@gt
              
              gt3 = handle_ploidy(gt, monoclonals = monoclonals, polyclonals = polyclonals)
              gt3 = as.data.frame(gt3)
              
              pi = NULL
              var = NULL
              
              for(region in 1:nrow(dna_regions)){
                positions = paste(dna_regions[region, ][['seqid']],dna_regions[region, ][['start']]:dna_regions[region, ][['end']], sep = '_')
                
                region_length = dna_regions[region, ][['end']] - dna_regions[region, ][['start']] + 1
                
                temp_gt = gt3[rownames(gt3) %in% positions,]
                
                n = ncol(temp_gt)
                
                if(nrow(temp_gt)>0){
                  
                  temp_loci = loci[rownames(loci) %in% positions,]
                  
                  temp_indels = temp_loci[temp_loci$Type_of_polymorphism != 'SNP',]
                  
                  temp_gt = temp_gt[temp_loci$Type_of_polymorphism == 'SNP',] # Keep only SNPs
                  
                  if(nrow(temp_gt)>0){
                    
                    if(length(temp_indels$ALT) > 0){
                      
                      ALTs = temp_indels$ALT
                      
                      ALTs = strsplit(ALTs, ',')
                      
                      gaps =  nchar(temp_indels$REF) - sapply(ALTs, function(alt){
                        max(nchar(alt))
                      })
                      
                      gaps[gaps < 0] = 0
                      
                      region_length = region_length - sum(gaps) # remove deletions from the total size of the region
                      
                    }
                    
                    haplotypes_counts =  summary(as.factor(sapply(1:ncol(temp_gt), function(x){paste(temp_gt[,x], collapse = '')})),
                                                 maxsum = ncol(temp_gt))
                    
                    if(length(haplotypes_counts) > 1){
                      
                      names(haplotypes_counts) = gsub('NA', "_", names(haplotypes_counts))
                      
                      haplotypes_freqs = haplotypes_counts/sum(haplotypes_counts)
                      haplotypes = names(haplotypes_counts)
                      
                      combinations = combn(1:length(haplotypes),2)
                      
                      temp_pi = NULL
                      
                      for(comb in 1:ncol(combinations)){
                        
                        x_i = haplotypes_freqs[combinations[1,comb]]
                        x_j = haplotypes_freqs[combinations[2,comb]]
                        
                        seq_i = unlist(strsplit(haplotypes[combinations[1,comb]], ''))
                        seq_i[seq_i == '_'] = NA
                        seq_j = unlist(strsplit(haplotypes[combinations[2,comb]], ''))
                        seq_j[seq_j == '_'] = NA
                        
                        pi_ij = sum(seq_i != seq_j, na.rm = T)/(region_length - sum(is.na(seq_i != seq_j)))
                        temp_pi = c(temp_pi, 2*x_i*x_j*pi_ij)
                      }
                      
                      
                      pi = c(pi, (n/(n-1))*(sum(temp_pi)))
                      var = c(var, (n + 1)*sum(temp_pi)/(3*(n - 1)*region_length) + 2*(n^2 + n + 3)*sum(temp_pi)^2/(9*n*(n - 1)))
                      
                    }else{
                      
                      pi = c(pi, 0)
                      var = c(var, 0)
                    }
                    
                  }else{
                    
                    pi = c(pi, NA)
                    var = c(var, NA)
                  }
                  
                }else{
                  
                  pi = c(pi, NA)
                  var = c(var, NA)
                  
                }
                
              }
              
              dna_regions$pi = pi
              dna_regions$pi_var = var
              
            }
            
            return(dna_regions)
          })


# get_EffCard----

setGeneric("get_EffCard", function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL) standardGeneric("get_EffCard"))

setMethod("get_EffCard", signature(obj = "rGenome"),
          function(obj = NULL, update_AC = FALSE, monoclonals = NULL, polyclonals = NULL){
  
            gt = obj@gt
            loci = obj@loci_table
            
            if(!update_AC){
              AC = sapply(1:nrow(loci), function(x){
                AC = strsplit(loci[x,'Allele_Counts'], ',')[[1]]
                gsub('^\\d+:', '', AC)
              })
            }else{
              AC = get_AC(gt = obj,w =1, n = 1, update_alleles = FALSE, monoclonals = monoclonals, polyclonals = polyclonals)
              AC = AC$Allele_Counts
              AC = sapply(AC, function(x){
                AC = strsplit(x, ',')[[1]]
                gsub('^\\d+:', '', AC)
              })
            }
  
            EffCard = sapply(1:length(AC), function(pos){
              n = sum(as.numeric(AC[[pos]]))
              allele_counts = as.numeric(AC[[pos]])
              allele_freq = allele_counts/n
              sp2 = sum(allele_freq^2)
              ExpHet = n * (1 - sp2)/(n - 1)
              EffCard = 1/(1-ExpHet)
            })
            
            return(EffCard)
            }
          )

# get_ObsHet----

setGeneric("get_ObsHet", function(obj = NULL, by = 'loci', w = 1, n = 1) standardGeneric("get_ObsHet"))

setMethod("get_ObsHet", signature(obj = "rGenome"),
          function(obj = NULL,
                   by = 'loci', #sample
                   w = 1, n = 1){
            
            gt = obj@gt
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            if(by == 'loci'){
              gt = gt[low:high,]
              ObsHet = rowSums(matrix(grepl('/', gt), nrow = nrow(gt), ncol = ncol(gt)))/ncol(gt)
            }else if(by == 'sample'){
              ObsHet = colSums(matrix(grepl('/', gt), nrow = nrow(gt), ncol = ncol(gt)))/nrow(gt)
            }else{
              print('by argument should be loci or sample')
            }
            
            
            return(ObsHet)
          }
  
)

# get_type_of_polymorphism----

get_type_of_polymorphism = function(obj, w = 1, n = 1){
  
  if(class(obj) == 'rGenome'){
    loci_table = obj@loci_table
  }else{
    loci_table = obj
  }
  
  s = round(seq(1,nrow(loci_table)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  loci_df = loci_table[low:high,c('REF', 'ALT')]
  
  alleles = apply(loci_df, 1, function(x){paste(x[1], x[2], sep = ',')})
  
  logical_vector = sapply(alleles,
                          function(site){
                            prod(nchar(strsplit(site, ',')[[1]]) == 1)
                          })
  
  del_vector = sapply(alleles,
                      function(site){
                        grepl('\\*', site)
                      })
  
  type_of_marker = ifelse(logical_vector == 1 & del_vector == FALSE, 'SNP', 'INDEL')
  
  indels = alleles[type_of_marker == 'INDEL']
  
  homopolymers = sapply(indels, function(site){
    ((length(strsplit(paste(gsub('^.','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]]) > 1)*
       length(unique(strsplit(paste(gsub('^.','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]])) == 1|
       (length(strsplit(paste(gsub('.$','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]]) > 1)*
       length(unique(strsplit(paste(gsub('.$','',strsplit(site, ',')[[1]]), collapse = ''), '')[[1]])) == 1)
  })
  
  type_of_marker[type_of_marker == 'INDEL'] = ifelse(homopolymers, 'INDEL:Homopolymer', 'INDEL')
  
  # short tandem repeats
  
  indels = alleles[type_of_marker == 'INDEL']
  
  STRs = sapply(indels, function(site){
    str = strsplit(site, ',')[[1]]
    
    str = gsub('^.', '', str[which.max(sapply(str, function(allele) {nchar(allele)}))])
    str2 = gsub('.$', '', str[which.max(sapply(str, function(allele) {nchar(allele)}))])
    
    indel = NULL
    
    if(nchar(str) < 2){
      indel = 'INDEL'
    }
    
    if(nchar(str) > 2 & nchar(str)%%2 == 0 & is.null(indel)){
      dinucletide = substring(str, seq(1, nchar(str), 2), seq(1, nchar(str), 2) + 1)
      dinucletide2 = substring(str2, seq(1, nchar(str2), 2), seq(1, nchar(str2), 2) + 1)
      if(prod(grepl(dinucletide[1], dinucletide)) == 1 | prod(grepl(dinucletide2[1], dinucletide2)) == 1){
        indel = 'INDEL:Dinucleotide_STR'
      }
    }
    
    if(nchar(str) > 3 & nchar(str)%%3 == 0 & is.null(indel)){
      trinucletide = substring(str, seq(1, nchar(str), 3), seq(1, nchar(str), 3) + 2)
      trinucletide2 = substring(str2, seq(1, nchar(str2), 3), seq(1, nchar(str2), 3) + 2)
      if(prod(grepl(trinucletide[1], trinucletide)) == 1| prod(grepl(trinucletide2[1], trinucletide2)) == 1)
        indel = 'INDEL:Trinucleotide_STR'
    }
    
    if(nchar(str) > 4 & nchar(str)%%4 == 0 & is.null(indel)){
      tetranucletide = substring(str, seq(1, nchar(str), 4), seq(1, nchar(str), 4) + 3)
      
      tetranucletide2 = substring(str2, seq(1, nchar(str2), 4), seq(1, nchar(str2), 4) + 3)
      if(prod(grepl(tetranucletide[1], tetranucletide)) == 1| prod(grepl(tetranucletide2[1], tetranucletide2)) == 1){
        indel = 'INDEL:Tetranucleotide_STR'
      }
    }
    
    if(nchar(str) > 5 & nchar(str)%%5 == 0 & is.null(indel)){
      pentanucletide = substring(str, seq(1, nchar(str), 5), seq(1, nchar(str), 5) + 4)
      pentanucletide2 = substring(str, seq(1, nchar(str2), 5), seq(1, nchar(str2), 5) + 4)
      if(prod(grepl(pentanucletide[1], pentanucletide)) == 1 | prod(grepl(pentanucletide2[1], pentanucletide2)) == 1){
        indel = 'INDEL:Pentanucleotide_STR'
      }
    }
    
    if(nchar(str) > 6 & nchar(str)%%6 == 0 & is.null(indel)){
      hexanucletide = substring(str, seq(1, nchar(str), 6), seq(1, nchar(str), 6) + 5)
      hexanucletide2 = substring(str2, seq(1, nchar(str2), 6), seq(1, nchar(str2), 6) + 5)
      if(prod(grepl(hexanucletide[1], hexanucletide)) == 1 | prod(grepl(hexanucletide2[1], hexanucletide2)) == 1){
        indel = 'INDEL:Hexanucleotide_STR'
      }
    }
    
    if(is.null(indel)){
      indel = 'INDEL'
    }
    
    indel
    
  })
  
  type_of_marker[type_of_marker == 'INDEL'] = STRs
  
  return(type_of_marker)
  
}

# Average Fraction of heterozygous samples per alternative allele per site----

setGeneric("frac_ofHet_pAlt", function(obj = NULL, w = 1, n = 1) standardGeneric("frac_ofHet_pAlt"))

setMethod("frac_ofHet_pAlt", signature(obj = "rGenome"),
          
          function(obj = NULL, w = 1, n = 1){
            
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w + 1] - 1
            
            alt = gsub('^\\d+,', '', gsub('(\\w+|\\*):', '', loci[low:high, 'Alleles']))
            gt = gsub(':\\d+', '',gt[low:high,])
            
            # Heterozygous positions
            HetPos = matrix(grepl('/', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            # For each variant site
            frac_ofHet_pAlt = sapply(1:nrow(gt), function(variant) {
              temp_gts = gt[variant,] # genotypes observed in that site
              alleles = strsplit(alt[variant], ',')[[1]] # alternative alleles observed in that site
              
              # Vector of presence or absence of each alternative allele
              haplotypes = sapply(alleles,
                                 function(allele){
                                   haplotypes = grepl(allele, temp_gts)})
              
              # Samples where alternative alleles are present and the site is heterozygous
              het_haplotypes = (haplotypes == 1 & HetPos[variant,] == 1)    
              
              sum(het_haplotypes, na.rm = T)/sum(haplotypes, na.rm = T)
            })
            
            return(frac_ofHet_pAlt)
            
          }
          )

# Fraction of heterozygous samples per alternative allele per site----
setGeneric("frac_ofHet_pAlt_byAllele", function(obj = NULL, w = 1, n = 1, add_variable = NULL) standardGeneric("frac_ofHet_pAlt_byAllele"))

setMethod("frac_ofHet_pAlt_byAllele", signature(obj = "rGenome"),
          
          function(obj = NULL, w = 1, n = 1, add_variable = NULL){
            
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w + 1] - 1
            
            loci = loci[low:high,]
            alt = gsub('^0,', '', gsub('(\\w+|\\*):', '', loci[, 'Alleles']))
            
            gt = gsub(':\\d+', '',gt[low:high,])
            
            # Heterozygous positions
            HetPos = matrix(grepl('/', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            allele_count_frac_ofHet_pAlt = NULL
            
            for(variant_sie in 1:nrow(gt)){
              temp_gts = gt[variant_sie,] # genotypes observed in that site
              alleles = strsplit(alt[variant_sie], ',')[[1]] # alternative alleles observed in that site
              
              # Vector of presence or absence of each alternative allele
              h_ij = t(sapply(alleles,
                              function(allele){
                                P_ij = grepl(allele, temp_gts)
                                
                                H_ijminor = grepl(paste0('/',allele), temp_gts)
                                
                                # Samples where alternative alleles are present and the site is heterozygous
                                H_ij = (P_ij == 1 & HetPos[variant_sie,] == 1)    
                                
                                H_ijminor = (H_ijminor == 1 & HetPos[variant_sie,] == 1)  
                                
                                c(sum(P_ij, na.rm = T), 
                                  sum(H_ij, na.rm = T), 
                                  sum(H_ijminor, na.rm = T),
                                  sum(H_ij, na.rm = T)/sum(P_ij, na.rm = T),
                                  ifelse(is.na(sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)), 
                                         0,
                                         sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T))
                                )
                                
                              }))
              
              
              allele_count_frac_ofHet_pAlt_temp = as.data.frame(cbind(alleles, h_ij))
              
              names(allele_count_frac_ofHet_pAlt_temp) = c('Allele', 'P_ij', 'H_ij', 'H_ijminor', 'h_ij', 'h_ijminor')
              
              allele_count_frac_ofHet_pAlt_temp[['P_ij']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['P_ij']])
              allele_count_frac_ofHet_pAlt_temp[['H_ij']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['H_ij']])
              allele_count_frac_ofHet_pAlt_temp[['H_ijminor']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['H_ijminor']])
              
              allele_count_frac_ofHet_pAlt_temp[['h_ij']] = as.numeric(allele_count_frac_ofHet_pAlt_temp[['h_ij']])
              allele_count_frac_ofHet_pAlt_temp[['h_ijminor']] = as.numeric(allele_count_frac_ofHet_pAlt_temp[['h_ijminor']])
              
              
              allele_count_frac_ofHet_pAlt_temp[['p_ij']] = allele_count_frac_ofHet_pAlt_temp[['P_ij']]/ncol(gt)
              
              allele_count_frac_ofHet_pAlt_temp[['VariantSite_id']] = rownames(gt)[variant_sie]
              
              if(!is.null(add_variable)){
                for(variable in add_variable){
                  allele_count_frac_ofHet_pAlt_temp[[variable]] = loci[variant_sie, ][[variable]]
                }
              }
              
              allele_count_frac_ofHet_pAlt = rbind(allele_count_frac_ofHet_pAlt, allele_count_frac_ofHet_pAlt_temp) 
              
              
            }
            
            return(allele_count_frac_ofHet_pAlt)
            
          }
)


# Mask alternative alleles----
setGeneric("mask_alt_alleles", function(obj = NULL, w = 1, n = 1, mask_formula = "h_ij >= 0.66 & h_ijminor >= 0.9 & p_ij <= 1") standardGeneric("mask_alt_alleles"))

setMethod("mask_alt_alleles", signature(obj = "rGenome"),
          
          function(obj = NULL, w = 1, n = 1, mask_formula = "h_ij >= 0.66 & h_ijminor >= 0.9 & p_ij <= 1"){
            
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w + 1] - 1
            
            loci = loci[low:high,]
            alt = gsub('^0,', '', gsub('(\\w+|\\*):', '', loci[, 'Alleles']))
            
            gt_masked = gt[low:high,]
            
            gt = gsub(':\\d+', '',gt[low:high,])
            
            # Heterozygous positions
            HetPos = matrix(grepl('/', gt), ncol = ncol(gt), nrow = nrow(gt))
            
            # Check if formula is correct
            
            if(grepl("(h_ij|h_ijminor|p_ij|P_ij|H_ij|H_ijminor)(<|>|!|=)+", mask_formula)){
              stop("All mathematical and logical operators must be separated by blank spaces in mask_formula")
            }
            
            # modify mask_formula
            
            if(grepl("h_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "h_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "h_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("h_ij ", "allele_count_frac_ofHet_pAlt_temp[['h_ij']] ", mask_formula)
              }else{
                stop("Filter h_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("h_ijminor ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "h_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "h_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("h_ijminor ", "allele_count_frac_ofHet_pAlt_temp[['h_ijminor']] ", mask_formula)
              }else{
                stop("Filter h_ijminor is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("p_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "p_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "p_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("p_ij ", "allele_count_frac_ofHet_pAlt_temp[['p_ij']] ", mask_formula)
              }else{
                stop("Filter p_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("P_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "P_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "P_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("P_ij ", "allele_count_frac_ofHet_pAlt_temp[['P_ij']] ", mask_formula)
              }else{
                stop("Filter P_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("H_ij ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "H_ij (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "H_ij (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("H_ij ", "allele_count_frac_ofHet_pAlt_temp[['H_ij']] ", mask_formula)
              }else{
                stop("Filter H_ij is been called but there are spelling issues in this part of the formula")
              }
            }
            
            if(grepl("H_ijminor ", mask_formula)){
              
              mask_filter = str_extract(mask_formula, "H_ijminor (=|!|>|<)+ (\\d+\\.?\\d*|\\d*\\.?\\d+)")
              
              if(!is.na(mask_filter)){
                print(paste0('Filter ', str_extract(mask_formula, "H_ijminor (=|!|>|<)* (\\d+\\.?\\d*|\\d*\\.?\\d+)"), ' will be applied'))
                mask_formula = gsub("H_ijminor ", "allele_count_frac_ofHet_pAlt_temp[['H_ijminor']] ", mask_formula)
              }else{
                stop("Filter H_ijminor is been called but there are spelling issues in this part of the formula")
              }
            }
            
            mask_formula_check = str_split(mask_formula, "&|\\|")[[1]]
            mask_formula_check  = mask_formula_check[!grepl("allele_count_frac_ofHet_pAlt_temp", mask_formula_check)]
            
            
            if(length(mask_formula_check) > 0){
              for(wrong_filter in mask_formula_check){
                print(paste0("Spelling error with filter ", wrong_filter))
              }
              stop("Execution halted, revise mask_filter argument.\nPossible filters are:\nh_ij, h_ijminor, p_ij, P_ij, H_ij, H_ijminor")
            }
            
            for(variant_site in 1:nrow(gt)){
              temp_gts = gt[variant_site,] # genotypes observed in that site
              alleles = strsplit(alt[variant_site], ',')[[1]] # alternative alleles observed in that site
              
              # Vector of presence or absence of each alternative allele
              h_ij = t(sapply(alleles,
                                         function(allele){
                                           P_ij = grepl(allele, temp_gts)
                                           
                                           H_ijminor = grepl(paste0('/',allele), temp_gts)
                                           
                                           # Samples where alternative alleles are present and the site is heterozygous
                                           H_ij = (P_ij == 1 & HetPos[variant_site,] == 1)    
                                           
                                           H_ijminor = (H_ijminor == 1 & HetPos[variant_site,] == 1)  
                                           
                                           c(sum(P_ij, na.rm = T), 
                                             sum(H_ij, na.rm = T), 
                                             sum(H_ijminor, na.rm = T),
                                             sum(H_ij, na.rm = T)/sum(P_ij, na.rm = T),
                                             ifelse(is.na(sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T)), 
                                                    0,
                                                    sum(H_ijminor, na.rm = T)/sum(H_ij, na.rm = T))
                                           )
                                           
                                         }))
              
              
              allele_count_frac_ofHet_pAlt_temp = as.data.frame(cbind(alleles, h_ij))
              
              names(allele_count_frac_ofHet_pAlt_temp) = c('Allele', 'P_ij', 'H_ij', 'H_ijminor', 'h_ij', 'h_ijminor')
              
              allele_count_frac_ofHet_pAlt_temp[['P_ij']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['P_ij']])
              allele_count_frac_ofHet_pAlt_temp[['H_ij']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['H_ij']])
              allele_count_frac_ofHet_pAlt_temp[['H_ijminor']] = as.integer(allele_count_frac_ofHet_pAlt_temp[['H_ijminor']])
              
              allele_count_frac_ofHet_pAlt_temp[['h_ij']] = as.numeric(allele_count_frac_ofHet_pAlt_temp[['h_ij']])
              allele_count_frac_ofHet_pAlt_temp[['h_ijminor']] = as.numeric(allele_count_frac_ofHet_pAlt_temp[['h_ijminor']])
              
              
              allele_count_frac_ofHet_pAlt_temp[['p_ij']] = allele_count_frac_ofHet_pAlt_temp[['P_ij']]/ncol(gt)
              
              allele_count_frac_ofHet_pAlt_temp[['VariantSite_id']] = rownames(gt)[variant_site]
              
              # Identify alleles below thresholds
              
              
              removed_alleles = allele_count_frac_ofHet_pAlt_temp[
                eval(parse(text = mask_formula)),][['Allele']]
              
              
              if(length(removed_alleles) > 0){
                removed_pattern =paste('/?(', paste(removed_alleles, collapse = '|'), '):\\d+/?', sep = '')
                # Mask alleles below threshold
                gt_masked[variant_site,] = gsub(removed_pattern, '', gt_masked[variant_site, ])
                
                gt_masked[variant_site, gt_masked[variant_site, ] == ''] = NA
                
              }
              
            }
            
            return(gt_masked)
            
          }
)


# Get haplotype matrix----

get_GT_matrix = function(vcf, w = 1, n = 1, threshold = 5){
  
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-8]
  
  gt = t(sapply(1:nrow(w_gt_table), function(variant) {
    gt_pos = grep('GT',strsplit(w_gt_table[variant,1], ':')[[1]])
    ad_pos = grep('AD',strsplit(w_gt_table[variant,1], ':')[[1]])
    temp_gts = w_gt_table[variant,-1]
    
    sapply(1:length(temp_gts), function(sample){
      gt = strsplit(strsplit(as.character(temp_gts[sample]), ':')[[1]][gt_pos], '/')[[1]]
      ad = strsplit(strsplit(as.character(temp_gts[sample]), ':')[[1]][ad_pos], ',')[[1]]
      
      if(gt[1] != '.'){
        ad = ad[as.numeric(gt) + 1]
        gt_df = data.frame(gt = gt, ad = as.numeric(ad))
        gt_df = gt_df[order(gt_df$ad, decreasing = T),]
        gt_df = gt_df[gt_df$ad >= threshold,]
        gt = gt_df$gt
        gt = unique(gt)
        
        if(length(gt) > 0){
          gt = paste(gt, collapse = '/')
        }else{
          gt = NA
        }
        
      }else{
        gt = NA
      }
      
    })
  }))
  
  rownames(gt) = vcf[low:high,1:2] %>% mutate(Locus = paste(CHROM, POS, sep = '_')) %>% select(Locus) %>% unlist
  colnames(gt) = colnames(vcf)[-1:-9]
  
  return(gt)
  
}

get_GTAD_matrix = function(vcf, w = 1, n = 100, threshold = 5){
  
  #suppressWarnings()
  start = Sys.time()
  s = round(seq(1,nrow(vcf)+1, length.out=n+1))
  low = s[w]
  high = s[w+1]-1
  
  w_gt_table = vcf[low:high,-1:-9]
  
  gt_temp = gsub(':.+', '', as.matrix(w_gt_table))
  ad_temp = gsub(':.+', '', gsub('^\\d+/\\d+:', '', as.matrix(w_gt_table)))
  
  if(!is.null(dim(w_gt_table))){
    
    gt = t(sapply(1:nrow(gt_temp), function(variant) {
      sapply(1:ncol(gt_temp), function(sample){
        gt = strsplit(as.character(gt_temp[variant, sample]), '/')[[1]]
        ad = strsplit(as.character(ad_temp[variant, sample]), ',')[[1]]
        
        if(gt[1] != '.'){
          ad = ad[as.numeric(gt) + 1]
          gt_df = data.frame(gt = gt, ad = as.numeric(ad))
          gt_df = gt_df[order(gt_df$ad, decreasing = T),]
          gt_df = gt_df[gt_df$ad >= threshold,]
          gt = paste(gt_df$gt, gt_df$ad, sep = ':')
          gt = unique(gt)
          
          if(length(gt) > 0){
            gt = paste(gt, collapse = '/')
          }else{
            gt = NA
          }
          
        }else{
          gt = NA
        }
        
      })
    }))
    
  }else{
    
    gt = matrix(sapply(1:nrow(gt_temp), function(variant) {
      sapply(1:ncol(gt_temp), function(sample){
        gt = strsplit(as.character(gt_temp[variant, sample]), '/')[[1]]
        ad = strsplit(as.character(ad_temp[variant, sample]), ',')[[1]]
        
        if(gt[1] != '.'){
          ad = ad[as.numeric(gt) + 1]
          gt_df = data.frame(gt = gt, ad = as.numeric(ad))
          gt_df = gt_df[order(gt_df$ad, decreasing = T),]
          gt_df = gt_df[gt_df$ad >= threshold,]
          gt = paste(gt_df$gt, gt_df$ad, sep = ':')
          gt = unique(gt)
          
          if(length(gt) > 0){
            gt = paste(gt, collapse = '/')
          }else{
            gt = NA
          }
          
        }else{
          gt = NA
        }
        
      })
    }), ncol = 1)
    
  }

  
  
  end = Sys.time()
  end - start
  
  rownames(gt) = vcf[low:high,1:2] %>% mutate(Locus = paste(CHROM, POS, sep = '_')) %>% select(Locus) %>% unlist
  colnames(gt) = colnames(vcf)[-1:-9]
  
  return(gt)
  
}



# Calculate the average read depth by gene----

mean_ReadDepth = function(data, gff = 'genes.gff'){
  
  ref_gff = ape::read.gff('genes.gff')
  coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                   !grepl('^Transfer',ref_gff$seqid),c('seqid', 'start', 'end', 'attributes')]
  
  coding_regions$attributes = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=PVP01_([0-9]+|MIT[0-9]+|API[0-9]+)'))
  
  coding_regions = coding_regions[order(coding_regions$start),]
  coding_regions = coding_regions[order(coding_regions$seqid),]
  rownames(coding_regions) = 1:nrow(coding_regions)
  
  coding_regions$mean_DP = NA
  
  for(gene in 1:nrow(coding_regions)){
    coding_regions[gene, ][['mean_DP']] = mean(data[data$CHROM == coding_regions[gene, ][['seqid']]&
            data$POS >= coding_regions[gene, ][['start']]&
            data$POS <= coding_regions[gene, ][['end']],][['DP']])
    
  }
  
  coding_regions = coding_regions[!is.na(coding_regions$mean_DP),]
  
  data$mean_DP = NA
  data$gene_id = NA
  
  for(gene in 1:nrow(coding_regions)){
    data[data$CHROM == coding_regions[gene, ][['seqid']]&
            data$POS >= coding_regions[gene, ][['start']]&
            data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['attributes']]
    
    data[data$CHROM == coding_regions[gene, ][['seqid']]&
           data$POS >= coding_regions[gene, ][['start']]&
           data$POS <= coding_regions[gene, ][['end']],][['mean_DP']] = coding_regions[gene, ][['mean_DP']]
  }
  
  return(data)
}

# mean_ObsHet ----

setGeneric("mean_ObsHet", function(obj = NULL, gff = 'genes.gff', Type_of_polymorphism = 'both', filter_sites = NULL) standardGeneric("mean_ObsHet"))

setMethod("mean_ObsHet", signature(obj = "rGenome"),
          
          function(obj = NULL, gff = 'genes.gff', Type_of_polymorphism = 'both', filter_sites = NULL){ # both, 'SNP', 'INDEL'
  
            data = obj@loci_table
            
            if(!is.null(filter_sites)){
              obj_filtered = filter_loci(obj, v = filter_sites)
            }else{
              obj_filtered = obj
            }
            
            data_filtered = obj_filtered@loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer',ref_gff$seqid),
                                     c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            coding_regions$mean_ObsHet = NA
            
            # Calculate the mean observed heterozygosity by each gene
            for(gene in 1:nrow(coding_regions)){
              
              if(Type_of_polymorphism == 'SNP'){
                coding_regions[gene, ][['mean_ObsHet']] = mean(data_filtered[data_filtered$CHROM == coding_regions[gene, ][['seqid']] &
                                                                               data_filtered$POS + nchar(data_filtered$REF) - 1 >= coding_regions[gene, ][['start']] &
                                                                               data_filtered$POS <= coding_regions[gene, ][['end']] &
                                                                               data_filtered$Cardinality > 1 &
                                                                               data_filtered$Type_of_polymorphism == 'SNP',][['ObsHet']])
                
              }else if(Type_of_polymorphism == 'INDEL'){
                coding_regions[gene, ][['mean_ObsHet']] = mean(data_filtered[data_filtered$CHROM == coding_regions[gene, ][['seqid']] &
                                                                               data_filtered$POS + nchar(data_filtered$REF) - 1 >= coding_regions[gene, ][['start']] &
                                                                               data_filtered$POS <= coding_regions[gene, ][['end']] &
                                                                               data_filtered$Cardinality > 1 &
                                                                      grepl('INDEL',data_filtered$Type_of_polymorphism),][['ObsHet']])
                
              }else if(Type_of_polymorphism == 'both'){
                coding_regions[gene, ][['mean_ObsHet']] = mean(data_filtered[data_filtered$CHROM == coding_regions[gene, ][['seqid']]&
                                                                               data_filtered$POS + nchar(data_filtered$REF) - 1 >= coding_regions[gene, ][['start']]&
                                                                               data_filtered$POS <= coding_regions[gene, ][['end']] &
                                                                               data_filtered$Cardinality > 1,][['ObsHet']])
                
              }
              
            }
            
            # Genes not found in the loci_table
            coding_regions2 = coding_regions[is.na(coding_regions$mean_ObsHet),]
            # Genes found in the loci_table
            coding_regions = coding_regions[!is.na(coding_regions$mean_ObsHet),]
            
            
            # Add mean observed heterozygosity of the gene to each variant site
            data$mean_ObsHet = NA
            data$gene_id = NA
            
            for(gene in 1:nrow(coding_regions)){
              data[data$CHROM == coding_regions[gene, ][['seqid']]&
                     data$POS + nchar(data$REF) - 1 >= coding_regions[gene, ][['start']]&
                     data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['attributes']]
              
              data[data$CHROM == coding_regions[gene, ][['seqid']]&
                     data$POS + nchar(data$REF) - 1 >= coding_regions[gene, ][['start']]&
                     data$POS <= coding_regions[gene, ][['end']],][['mean_ObsHet']] = coding_regions[gene, ][['mean_ObsHet']]
            }
            
            
            
            if(sum(is.na(data$mean_ObsHet))>0){
              data[is.na(data$mean_ObsHet), ][['mean_ObsHet']] = 0
            }
            
            
            for(pos in rownames(data[is.na(data$gene_id), ])){
              
              attrib = coding_regions2[coding_regions2[['seqid']] == data[pos,][['CHROM']] &
                                         coding_regions2[['start']] <= data[pos,][['POS']] &
                                         coding_regions2[['end']] >= data[pos,][['POS']],][['attributes']]
              
              data[pos,][['gene_id']] = ifelse(length(attrib) == 0, NA, attrib)
              
            }
            
            return(data)
          }
  
)



# get_gene_description ----

get_gene_description_rGenome = function(obj = NULL, gff = 'genes.gff'){
            
            if(class(obj) == 'rGenome'){
              data = obj@loci_table
            }else{
              data = obj
            }
  
            
            data$POS_end = data$POS + nchar(data$REF) - 1
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer',ref_gff$seqid),
                                     c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$gene_id = gsub('ID=','',str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions$gene_description = gsub('description=','',str_extract(coding_regions$attributes, '(description=.+$|description=.+;)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            
            data$gene_id = NA
            data$gene_description = NA
            
            
            for(gene in 1:nrow(coding_regions)){
              
              if(nrow(data[data$CHROM == coding_regions[gene, ][['seqid']]&
                           data$POS >= coding_regions[gene, ][['start']]&
                           data$POS <= coding_regions[gene, ][['end']],]) != 0){
                
                data[data$CHROM == coding_regions[gene, ][['seqid']]&
                       data$POS >= coding_regions[gene, ][['start']]&
                       data$POS <= coding_regions[gene, ][['end']],][['gene_id']] = coding_regions[gene, ][['gene_id']]
                
                data[data$CHROM == coding_regions[gene, ][['seqid']]&
                       data$POS >= coding_regions[gene, ][['start']]&
                       data$POS <= coding_regions[gene, ][['end']],][['gene_description']] = coding_regions[gene, ][['gene_description']]
              }
              
            }
            
            # indels that start out of the gene region but end 
            
            for(pos in rownames(data[is.na(data$gene_id), ])){
              
              if(nrow(coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                     coding_regions[['start']] <= data[pos,][['POS_end']] &
                                     coding_regions[['end']] >= data[pos,][['POS_end']],]) != 0){
                data[pos,][['gene_id']] = coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                                           coding_regions[['start']] <= data[pos,][['POS_end']] &
                                                           coding_regions[['end']] >= data[pos,][['POS_end']],][['gene_id']]
                
                data[pos,][['gene_description']] = coding_regions[coding_regions[['seqid']] == data[pos,][['CHROM']] &
                                                                    coding_regions[['start']] <= data[pos,][['POS_end']] &
                                                                    coding_regions[['end']] >= data[pos,][['POS_end']],][['gene_description']]
              }
              
            }
            
            return(data[, c('gene_id', 'gene_description')])
          }
          

# SNP_density----

setGeneric("SNP_density", function(obj = NULL, gff = 'genes.gff', filter_sites = NULL) standardGeneric("SNP_density"))

setMethod("SNP_density", signature(obj = "rGenome"),
          
          function(obj = NULL, gff = 'genes.gff', filter_sites = NULL){
            
            data = obj@loci_table
            
            if(!is.null(filter_sites)){
              obj_filtered = filter_loci(obj, v = filter_sites)
            }else{
              obj_filtered = obj
            }
            
            data_filtered = obj_filtered@loci_table
            
            ref_gff = ape::read.gff(gff)
            coding_regions = ref_gff[grepl('gene', ref_gff$type)&
                                       !grepl('^Transfer', ref_gff$seqid), c('seqid', 'start', 'end', 'attributes')]
            
            coding_regions$attributes = gsub('ID=', '', str_extract(coding_regions$attributes, 'ID=(PVP01|PF3D7)_([0-9]+|MIT[0-9]+|API[0-9]+)'))
            
            coding_regions = coding_regions[order(coding_regions$start),]
            coding_regions = coding_regions[order(coding_regions$seqid),]
            rownames(coding_regions) = 1:nrow(coding_regions)
            
            coding_regions$SNP_density = NA
            
            for(gene in 1:nrow(coding_regions)){
              coding_regions[gene, ][['SNP_density']] = nrow(data_filtered[data_filtered$CHROM == coding_regions[gene, ][['seqid']] &
                                                                             data_filtered$POS + nchar(data_filtered$REF) - 1  >= coding_regions[gene, ][['start']] &
                                                                             data_filtered$POS <= coding_regions[gene, ][['end']] &
                                                                             data_filtered$Cardinality > 1 &
                                                                      data_filtered$Type_of_polymorphism == 'SNP',])/(coding_regions[gene, ][['end']] - coding_regions[gene, ][['start']] + 1)
              
            }
            
            data$SNP_density = NA
            
            for(gene in unique(data[!is.na(data$gene_id),][['gene_id']])){
              
              data[data$gene_id == gene &
                     !is.na(data$gene_id),][['SNP_density']] = coding_regions[coding_regions$attributes == gene, ][['SNP_density']]
            }
            
            if(sum(is.na(data$SNP_density)) > 0){
              data[is.na(data$SNP_density), ][['SNP_density']] = 0
            }
            
            return(data)
          }
)

# Summarise_ReadDepth----

Summarise_ReadDepth = function(obj, by = NA, w = 1, n = 100){
  
  s = round(seq(1, nrow(obj@gt) + 1, length.out = n + 1))
  low = s[w]
  high = s[w + 1] - 1
  
  metadata = obj@metadata
  
  if(!is.na(by)){
    
    populations = t(table(metadata[[by]]))
    populations = data.frame(population = colnames(populations), nsamples = populations[1, ])
    
    ReadDepth_Summ = NULL
    
    # For each subpopulation
    
    for(pop in populations$population){
      
      samples = metadata[metadata[[by]] == pop & !is.na(metadata[[by]]),][['Sample_id']]
      temp_pop = filter_samples(obj = obj, v = samples)
      
      ad = gsub('\\d+:', '', as.matrix(temp_pop@gt[low:high,]))
      
      ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
      ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
      
      ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
      
      temp_ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                       mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)}),
                                       sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                       median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                       quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                       quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                       iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                       )
      
      names(temp_ReadDepth_Summ) = paste(names(temp_ReadDepth_Summ), pop, sep = "_")
      
      if(is.null(ReadDepth_Summ)){
        ReadDepth_Summ = temp_ReadDepth_Summ
      }else{
        ReadDepth_Summ = cbind(ReadDepth_Summ, temp_ReadDepth_Summ)
      }
      
    }
    
    # For the total population

    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))
    
    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    
    ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
    
    temp_ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                     mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)}),
                                     sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                     median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                     quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                     quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                     iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                     )
    
    names(temp_ReadDepth_Summ) = paste(names(temp_ReadDepth_Summ), 'Total', sep = "_")
    
    ReadDepth_Summ = cbind(ReadDepth_Summ, temp_ReadDepth_Summ)
    rownames(ReadDepth_Summ) = rownames(ad)
    
  }else{
    
    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))
    
    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    
    ad3 = matrix(rowSums(cbind(c(ad1),c(ad2)), na.rm = T), nrow = nrow(ad), ncol = ncol(ad))
    
    ReadDepth_Summ = data.frame(Total_ReadDepth = rowSums(ad3, na.rm = T),
                                mean_ReadDepth = apply(ad3, 1, function(x){mean(x, na.rm = T)}),
                                sd_ReadDepth = apply(ad3, 1, function(x){sd(x, na.rm = T)}),
                                median_ReadDepth = apply(ad3, 1, function(x){median(x, na.rm = T)}),
                                quantile25_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .25, na.rm = T)}),
                                quantile75_ReadDepth = apply(ad3, 1, function(x){quantile(x, probs = .75, na.rm = T)}),
                                iqr_ReadDepth = apply(ad3, 1, function(x){IQR(x, na.rm = T)})
                                )
    
    rownames(ReadDepth_Summ) = rownames(ad)

  }
  
  return(ReadDepth_Summ)
  
}

# prune_alleles----

prune_alleles = function(obj, threshold = 4, n = 100){
  library(svMisc)
  s = round(seq(1,nrow(obj@gt)+1, length.out=n+1))
  gt6 = NULL
  
  for(w in 1:n){
    low = s[w]
    high = s[w+1]-1
  
    gt = gsub(':\\d+', '', as.matrix(obj@gt[low:high,]))
    gt1 = matrix(as.integer(gsub('/\\d+','',gt)), nrow = nrow(gt), ncol = ncol(gt))
    gt2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',gt)), nrow = nrow(gt), ncol = ncol(gt))
  
    ad = gsub('\\d+:', '', as.matrix(obj@gt[low:high,]))

    ad1 = matrix(as.integer(gsub('/\\d+','',ad)), nrow = nrow(ad), ncol = ncol(ad))
    ad2 = matrix(as.integer(gsub('^\\d+|^\\d+/','',ad)), nrow = nrow(ad), ncol = ncol(ad))
  
    # prune alleles
    gt1[ad1 <= threshold] = NA
    ad1[ad1 <= threshold] = NA
  
    gt2[ad2 <= threshold] = NA
    ad2[ad2 <= threshold] = NA
  
    gt3 = matrix(paste(gt1, ad1, sep = ':'), nrow = nrow(gt), ncol = ncol(gt))
    gt3[gt3 == 'NA:NA'] = NA
  
    gt4 = matrix(paste(gt2, ad2, sep = ':'), nrow = nrow(gt), ncol = ncol(gt))
    gt4[gt4 == 'NA:NA'] = NA
  
    gt5 = matrix(paste(gt3, gt4, sep = '/'), nrow = nrow(gt), ncol = ncol(gt))
    gt5[gt5 == 'NA/NA'] = NA
    gt5 = gsub('/NA', '', gt5)
  
    colnames(gt5) = colnames(gt)
    rownames(gt5) = rownames(gt)
    
    gt6 = rbind(gt6, gt5)
    
    progress(round(100*w/n))
    
  }
  
  return(gt6)
  
  }

# get_Fws_rGenome----

setGeneric("get_Fws_rGenome", function(obj = NULL, w = 1, n = 1) standardGeneric("get_Fws_rGenome"))

setMethod("get_Fws_rGenome", signature(obj = "rGenome"),
          function(obj = NULL, w = 1, n = 1){
  
            gt = obj@gt
            loci = obj@loci_table
            
            s = round(seq(1,nrow(gt)+1, length.out=n+1))
            low = s[w]
            high = s[w+1]-1
            
            gt = gt[low:high,]
            
            ExpHet = get_ExpHet(obj = obj, update_AC = TRUE)
            
            Hw = sapply(1:ncol(gt), function(sample){
              
              samp_alleles= gsub(':\\d+', '', gt[,sample])
              samp_allcounts = gsub('\\d+:', '', gt[,sample])
              
              samp_alleles1 = gsub('/\\d+$', '', samp_alleles)
              samp_alleles2 = gsub('^\\d+/', '', samp_alleles)
              
              samp_check = samp_alleles1 != samp_alleles2
              
              samp_alleles2[!samp_check] = NA
              
              samp_allcounts1 = gsub('/\\d+$', '', samp_allcounts)
              samp_allcounts2 = gsub('^\\d+/', '', samp_allcounts)
              
              samp_allcounts2[!samp_check] = NA
              
              samp_allcountsT = rowSums(cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2)), na.rm = T)
              
              samp_allfreq = cbind(as.integer(samp_allcounts1), as.integer(samp_allcounts2))/samp_allcountsT
              
              Hw = 1 - rowSums(samp_allfreq^2, na.rm = T)
              Hw[Hw==1] = NA
              Hw
            })
            
            Fws = 1 - (Hw/ExpHet)
            
            colMeans(Fws, na.rm = T)
            
          }
)

# get_genclon----

get_genclone = function(gt,
                    loci,
                    monoclonals,
                    polyclonals,
                    exclude_indels = T,
                    window = 150,
                    metadata){
  
  library(poppr)
  library(pegas)
  
  if(exclude_indels){
    gt = gt[loci$Type_of_polymorphism == 'SNP',]
    loci %<>% filter(Type_of_polymorphism == 'SNP')
    
  }
  
  gt = matrix(gsub(':\\d+', '', gt), nrow = nrow(gt), ncol = ncol(gt),
              dimnames = list(rownames(gt), colnames(gt)))
  
  if(is.null(monoclonals) & is.null(polyclonals)){
    gt1 = gsub('/\\d+', '', gt)
    gt2 = gsub('\\d+/', '', gt)
    gt2[gt2==gt1] = NA
    gt3 = cbind(gt1, gt2)
  }else{
    gt_mono = gt[,monoclonals]
    gt_mono = gsub('/\\d+', '', gt_mono)
    
    gt_poly = gt[,polyclonals]
    gt_poly1 = gsub('/\\d+', '', gt_poly)
    gt_poly2 = gsub('\\d+/', '', gt_poly)
    gt3 = cbind(gt_mono, gt_poly1, gt_poly2)
  }
  
  loci = loci[,c('CHROM', 'POS')]
  
  chrom_length = loci %>% group_by(CHROM) %>% summarise(length = max(POS))
  
  chrom_intervals = sapply(chrom_length$length, function(chrom){
    seq(1, chrom, window)
  })
  
  
  chrom_win = NULL
  
  for(chrom in 1:length(chrom_intervals)){
    chrom_win = rbind(chrom_win, data.frame(CHROM = chrom_length[chrom,][['CHROM']],
                                            start = chrom_intervals[[chrom]],
                                            end = chrom_intervals[[chrom]] - 1 + window))
    
  }
  options(scipen=999)
  
  window_matrix = NULL
  
  for(i in 1:window){
    window_matrix = cbind(window_matrix,paste(chrom_win[['CHROM']], as.character(chrom_win[['start']] + i - 1), sep = '_'))
  }
  
  options(scipen=0)
  
  chrom_win$nVar = rowSums(matrix(window_matrix %in% rownames(loci),
                                  ncol = ncol(window_matrix),
                                  nrow = nrow(window_matrix)), na.rm = T)
  
  chrom_win %<>% filter(nVar != 0)
  
  chrom_win %<>% mutate(Filter = nVar >=3)
  
  Filtered_pos  = rep(chrom_win[['Filter']], chrom_win[['nVar']])
  
  chrom_win%<>%filter(nVar >=3)
  
  loci = loci[Filtered_pos,]
  
  gt3 = gt3[Filtered_pos,]
  
  collapsed_gt3 = NULL
  
  for(hap in 1:ncol(gt3)){
    start = Sys.time()
    collapsed_gt3 = rbind(collapsed_gt3,sapply(1:nrow(chrom_win), function(win){
      
      if(win > 1){
        bin = 1:chrom_win[win,][['nVar']] + sum(chrom_win[1:win - 1,][['nVar']])
      }else{
        bin = 1:chrom_win[win,][['nVar']]
      }
      paste(gt3[bin,hap], collapse = '')
    }))
    
    end = Sys.time()
    print(paste(hap, 'in', end - start))
  }
  
  collapsed_gt3[grepl('NA',collapsed_gt3)] = NA
  
  colnames(collapsed_gt3) = paste(chrom_win$CHROM, chrom_win$start, chrom_win$end, sep = '_')
  rownames(collapsed_gt3) = colnames(gt3)
  
  loci_format = pegas::as.loci(collapsed_gt3)
  
  genind_format = loci2genind(loci_format, ploidy = 1)
  
  genclone_format = as.genclone(genind_format)
  
  genclone_format@strata = data.frame(sample = rownames(collapsed_gt3))
  
  genclone_format@strata = left_join(genclone_format@strata,
                                     metadata,
                                     by = 'sample')
  
  return(genclone_format)
  
}

# filter_gt_matrix ----

filter_gt_matrix = function(gt, # haplotype matrix
                            loci, # table with loci information
                            filter_table, # Table of filtered segments or positions
                            by = 'segments',
                            keep = TRUE,
                            exclude_indels = TRUE){
  
  if(exclude_indels){
    gt = gt[loci$Type_of_polymorphism == 'SNP',]
    loci %<>% filter(Type_of_polymorphism == 'SNP')
  }
  
  
  positions = NULL
  
  options(scipen=999)
  
  for(pos in 1:nrow(filter_table)){
    positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                   filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
    
  }
  
  options(scipen=0)
  
  rownames(loci) %in% positions
 
  gt = gt[rownames(loci) %in% positions,]
  
  return(gt)
  
}




# fastGRMcpp----
#' C++ implementation of a Genomic relationship matrix 'GRM'
#' 
#' @param X Matrix of the type 'MatrixXd' for which the GRM will be calculated.
#' 
#' @return Genomic relationship matrix (GRM).
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' 
#' @export
#' 

grm = function(X){
  grmCpp(X)
}


#' C++ implementation of a fast singular value decomposition (SVD)
#' 
#' @param X Symmetric matrix of the type 'MatrixXd' for which the SVD will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = grm(Data)
#' V = fastSVD(X, 2)
#' 
#' @export

fastSVD = function(X, k, q = 2){
  
  
  
  fastSVDCpp(X, k, q)
}

#' C++ implementation of a fast GRM function
#' 
#' @param X Matrix of the type 'MatrixXd' for which the fastGRM matrix will be calculated.
#' @param k Number of first k eigen vectors to return
#' @param q Auxiliary exponent
#' 
#' @return SVD matrix of size .
#' 
#' @importFrom Rdpack reprompt
#' @references https://doi.org/10.1016/j.ajhg.2010.11.011
#' @references https://doi.org/10.48550/arXiv.0909.4061
#' 
#' @examples
#' require(fastGRM)
#' Data = matrix(sample(0:1, 9000, TRUE, c(.9,.1)), 90)
#' X = fastGRM(Data, 2)
#' 
#' @export
#' 

fastGRM = function(obj, k, monoclonals = NULL, polyclonals = NULL, Pop = NULL, q = 2){
  
  gt = obj@gt
  loci = obj@loci_table
  metadata = obj@metadata
  
  
  X = handle_ploidy(gt = gt, w = 1, n = 1, monoclonals = monoclonals, polyclonals = polyclonals)
  
  
  X = matrix(as.numeric(X), ncol = ncol(X),
                                        nrow = nrow(X), 
                                        dimnames = list(
                                          rownames(X),
                                          colnames(X)
                                        ))
  
  X[is.na(X)] = 0
  
  ibs_matrix = grm(X)
  
  ibs_pca = princomp(ibs_matrix)
  
  evector = ibs_pca$scores
  evalues = ibs_pca$sdev
  
  contrib = 100*(evalues)^2/sum((evalues)^2)
  
  # evector = fastGRMCpp(X, k, q)
  # 
  # evalues = NULL
  # 
  # for(i in 1:k){
  #   evalues = c(evalues, unlist((ibs_matrix %*% evector[,i])/evector[,i])[1])
  # }
  # 
  # contrib = 100*(abs(evalues)/sum(abs(evalues)))
  
  #### Add metadata to the PCA----
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(X)),
                             order = 1:ncol(X)), metadata[,c('Sample_id', Pop)], by = 'Sample_id', all.x = T)
  
  Pop_col = Pop_col[order(Pop_col$order),]
  
  evector = data.frame(Pop_col, evector)
  names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))
  
  ibs_pca = list(eigenvectors = evector, eigenvalues = evalues, contrib = contrib)
  
  return(ibs_pca)
  
}

Pop = 'site_of_collection_world_region'

GRM_evectors = function(obj, k = 10, monoclonals = NULL, polyclonals = NULL, Pop = NULL, q = 2, cor = T, method = 'fastGRM'){

    gt = obj@gt
    loci = obj@loci_table
    metadata = obj@metadata


    X = handle_ploidy(gt = gt, w = 1, n = 1, monoclonals = monoclonals, polyclonals = polyclonals)


    X = matrix(as.numeric(X), ncol = ncol(X),
               nrow = nrow(X),
               dimnames = list(
                 rownames(X),
                 colnames(X)
               ))

    X[is.na(X)] = 0

    ibs_matrix = grm(X)


  #### Add metadata to the PCA----
  Pop_col = merge(data.frame(Sample_id = gsub('_C[1,2]$','',colnames(X)),
                             order = 1:ncol(X)), metadata[,c('Sample_id', Pop)], by = 'Sample_id', all.x = T)

  Pop_col = Pop_col[order(Pop_col$order),]
 

  if(method == 'fastGRM'){

    evector = fastSVDCpp(ibs_matrix, k, q)

    for(i in 1:k){
      evector[,i] = sign(evector[1,i])*evector[,i]
    }

    evalues = NULL

    for(i in 1:k){
      evalues = c(evalues, unlist((ibs_matrix %*% evector[,i])/evector[,i])[1])
    }

    contrib = 100*(evalues)^2/sum((evalues)^2)

    evector = data.frame(Pop_col, evector)
    names(evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))

    ibs_pca = list(eigenvector=evector, eigenvalues = evalues, contrib = contrib)

  }else if(method == 'princomp'){

    ## Using princomp from R

    ibs_pca = princomp(ibs_matrix, cor = cor)
    ibs_evector = ibs_pca$scores
    ibs_evalues = ibs_pca$sdev

    ibs_contrib = 100*(ibs_evalues)^2/sum((ibs_evalues)^2)

    ibs_evector = data.frame(Pop_col, ibs_evector)
    names(ibd_evector) = c(colnames(Pop_col), paste0(rep('PC', k), 1:k))

    ibs_pca = list(eigenvector = ibs_evector, eigenvalues = ibs_evalues, contrib = ibs_contrib)

  }

  return(ibs_pca)

}

# merge_rGenome----

merge_rGenome = function(obj1 = NULL, obj2 = NULL, na.rm = T){
  
  # Check if the column Alleles is not present in the first data set
  if(sum(grepl('Alleles', names(obj1@loci_table))) == 0){ # if not generate Alleles and Allele_Counts
    # count Alleles
    Locus_info_temp1 = get_AC(obj1, update_alleles = T)
    Locus_info_temp1 = cbind(obj1@loci_table, Locus_info_temp1)
  }else{
    Locus_info_temp1 = obj1@loci_table[, c('CHROM',
                                           'POS',
                                           'REF', 
                                           'ALT',
                                           'Alleles',
                                           'Allele_Counts')]
  }
  
  # Check if the column Alleles is not present in the second data set
  if(sum(grepl('Alleles', names(obj2@loci_table))) == 0){ # if not generate Alleles and Allele_Counts
    # count Alleles
    Locus_info_temp2 = get_AC(obj2, update_alleles = T)
    Locus_info_temp2 = cbind(obj2@loci_table, Locus_info_temp2)
  }else{
    Locus_info_temp2 = obj2@loci_table[, c('CHROM',
                                           'POS',
                                           'REF', 
                                           'ALT',
                                           'Alleles',
                                           'Allele_Counts')]
  }
  
  # Relabel columns of the first data set
  
  Locus_info_temp1 = data.frame(locus_id = rownames(Locus_info_temp1), Locus_info_temp1)
  names(Locus_info_temp1) = c('locus_id',paste0(names(Locus_info_temp1)[-1], '.temp1'))
  
  # Relabel columns of the second data set
  Locus_info_temp2 = data.frame(locis_id = rownames(Locus_info_temp2), Locus_info_temp2)
  names(Locus_info_temp2) = c('locus_id',paste0(names(Locus_info_temp2)[-1], '.temp2'))
  
  
  # Store the genotype tables of both data sets
  gt_temp1 = obj1@gt
  gt_temp2 = obj2@gt
  
  # Merge the loci_table of both data sets
  Locus_info_merged = merge(Locus_info_temp1, Locus_info_temp2, by = 'locus_id', all = T)
  
  # Check if the reference allele is the same for both data sets in each position
  Locus_info_merged %<>% mutate(
    REF = case_when(
      REF.temp1 == REF.temp2 ~ REF.temp1,
      !is.na(REF.temp1) & is.na(REF.temp2) ~ REF.temp1,
      is.na(REF.temp1) & !is.na(REF.temp2) ~ REF.temp2,
      REF.temp1 != REF.temp2 ~ 'ERROR'
    )
  )
  
  # Check if the alternative alleles  and their labels are the same for both data sets in each position
  Locus_info_merged$Alleles = apply(Locus_info_merged, 1, function(pos){
    if(!is.na(pos['Alleles.temp1'])&is.na(pos['Alleles.temp2'])){ # If the second data set is all missing data
      paste0(pos['Alleles.temp1'], ';Only First data set') 
    }else if(is.na(pos['Alleles.temp1'])&!is.na(pos['Alleles.temp2'])){ # If the first data set is all missing data
      paste0(pos['Alleles.temp2'], ';Only Second data set') 
    }else if(is.na(pos['Alleles.temp1'])&is.na(pos['Alleles.temp2'])){ # If the first and second data sets are all missing data
      NA
    }else if(pos['Alleles.temp1'] == pos['Alleles.temp2']){ # If the alternative alleles  and their labels are the same for both data sets
      pos['Alleles.temp1']
    }else if(sum(!(unlist(str_split(pos['Alleles.temp2'], ',', simplify = T)) %in% unlist(str_split(pos['Alleles.temp1'], ',', simplify = T)))) == 0){ # If the alternative alleles  and their labels of second data set are contained in the first data set
      pos['Alleles.temp1']
    }else if(sum(!(unlist(str_split(pos['Alleles.temp1'], ',', simplify = T)) %in% unlist(str_split(pos['Alleles.temp2'], ',', simplify = T)))) == 0){ # If the alternative alleles  and their labels of first data set are contained in the second data set
      pos['Alleles.temp2']
    }else{'ERROR'} # If at least one alternative allele  and its label does not coincide between both data sets
  })
  
  # Split the loci table between:
  
  ## loci which reference alleles do not coincide between data sets
  Locus_info_merged_w_REFdiscrepancies = Locus_info_merged %>% filter(REF == 'ERROR')
  
  ## loci which reference alleles coincide but the alternative alleles do not coincide between data sets
  Locus_info_merged_w_ALTdiscrepancies = Locus_info_merged %>% filter(Alleles == 'ERROR', REF != 'ERROR')
  
  ## loci which reference and alternative alleles coincide between data sets and do not have missing data
  Locus_info_merged_good = Locus_info_merged %>% filter(Alleles != 'ERROR' & REF != 'ERROR' & !grepl('Only',Alleles))
  
  ## loci with no amplification in one data set
  Locus_info_merged_w_missing = Locus_info_merged %>% filter(grepl('Only',Alleles) & REF != 'ERROR')
  
  # Add row names for each splited set of loci
  rownames(Locus_info_merged_w_REFdiscrepancies) = Locus_info_merged_w_REFdiscrepancies$locus_id
  rownames(Locus_info_merged_w_ALTdiscrepancies) = Locus_info_merged_w_ALTdiscrepancies$locus_id
  rownames(Locus_info_merged_good) = Locus_info_merged_good$locus_id
  rownames(Locus_info_merged_w_missing) = Locus_info_merged_w_missing$locus_id
  
  # Relabel each position in the set of loci that reference alleles coincide but the alternative alleles do not coincide between data sets
  
  if(length(rownames(Locus_info_merged_w_ALTdiscrepancies)) > 0){
    
    for(pos in rownames(Locus_info_merged_w_ALTdiscrepancies)){
      
      # Get Alleles, allele labels, and allele counts for the first data set
      Alleles.temp1 = as.character(unlist(str_split(gsub(':\\d+','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp1']]), ',', simplify = T)))
      Allele_labels.temp1 = unlist(str_split(gsub('([ATGC]+|\\*):','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp1']]), ',', simplify = T))
      Allele_counts.temp1 = as.integer(unlist(str_split(gsub('\\d+:','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Allele_Counts.temp1']]), ',', simplify = T)))
      names(Alleles.temp1) = Allele_labels.temp1
      names(Allele_counts.temp1) = Allele_labels.temp1
      
      # Get Alleles, allele labels, and allele counts for the first data set
      Alleles.temp2 = as.character(unlist(str_split(gsub(':\\d+','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp2']]), ',', simplify = T)))
      Allele_labels.temp2 = unlist(str_split(gsub('([ATGC]+|\\*):','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles.temp2']]), ',', simplify = T))
      Allele_counts.temp2 = as.integer(unlist(str_split(gsub('\\d+:','',Locus_info_merged_w_ALTdiscrepancies[pos,][['Allele_Counts.temp2']]), ',', simplify = T)))
      names(Alleles.temp2) = Allele_labels.temp2
      names(Allele_counts.temp2) = Allele_labels.temp2
      
      # Get the reference and unique alternative alleles observed in both data sets
      REF.Allele = unique(c(Alleles.temp1[names(Alleles.temp1)[names(Alleles.temp1) == '0']], Alleles.temp2[names(Alleles.temp2)[names(Alleles.temp2) == '0']]))
      ALT.Alleles = unique(c(Alleles.temp1[names(Alleles.temp1)[names(Alleles.temp1) != '0']], Alleles.temp2[names(Alleles.temp2)[names(Alleles.temp2) != '0']]))
      
      # If reference allele is not present in either set
      if(length(REF.Allele) == 0){
        REF.Allele = Locus_info_merged_w_ALTdiscrepancies[pos,][['REF']]
      }
      
      # Calculate the total number of samples that have each ealternative allele in both data sets
      ALT.Allele_counts = NULL
      
      for(allele in ALT.Alleles){ # For each alternative allele
        
        ALT.Allele_counts = c(ALT.Allele_counts,
                              
                              if(allele %in% Alleles.temp1 & allele %in% Alleles.temp2){ # if the allele is present in both data sets
                                Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == allele])] +  
                                  Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == allele])]
                              }else if(allele %in% Alleles.temp1 & !(allele %in% Alleles.temp2)){ # if allele is present only in the first data set
                                Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == allele])]
                              }else if(!(allele %in% Alleles.temp1) & allele %in% Alleles.temp2){ # if allele is present only in the second data set
                                Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == allele])]
                              }
        )
      }
      
      # name (index) the allele count with the allele
      names(ALT.Allele_counts) = ALT.Alleles
      
      # sort the alternative alleles based on their allele count
      ALT.Allele_counts = sort(ALT.Allele_counts, decreasing = T)
      
      # Calculate the total number of samples that have the reference allele in both data sets
      REF.Allele_count = if(REF.Allele %in% Alleles.temp1 & REF.Allele %in% Alleles.temp2){ # if the allele is present in both data sets
        Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == REF.Allele])] +  
          Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == REF.Allele])]
      }else if(REF.Allele %in% Alleles.temp1 & !(REF.Allele %in% Alleles.temp2)){ # if allele is present only in the first data set
        Allele_counts.temp1[names(Alleles.temp1[Alleles.temp1 == REF.Allele])]
      }else if(!(REF.Allele %in% Alleles.temp1) & REF.Allele %in% Alleles.temp2){ # if allele is present only in the second data set
        Allele_counts.temp2[names(Alleles.temp2[Alleles.temp2 == REF.Allele])]
      }
      
      if(!is.null(REF.Allele_count)){
        
        # name (index) the allele count with the allele
        names(REF.Allele_count) = REF.Allele
        # Combine the reference and the alternative alleles in one object
        Allele_counts = c(REF.Allele_count, ALT.Allele_counts)
        Alleles = names(Allele_counts)
        
      }else{
        
        Allele_counts = ALT.Allele_counts
        Alleles = names(Allele_counts)
        
      }
      
      # RELABEL the alleles for the combined data set
      
      if(!is.null(REF.Allele_count)){
        Allele_labels = as.character(0:(length(Alleles) - 1))
        names(Alleles) = Allele_labels
      }else{
        
        Allele_labels = as.character(1:length(Alleles))
        names(Alleles) = Allele_labels
      }
      
      # If the alternative alleles  and their labels of First data set are NOT contained in the merged data set
      if(sum(!(paste(Alleles.temp1, Allele_labels.temp1, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))) != 0){
        
        # Identify and select the alleles which lables have changed
        Changed.Alleles.temp1 = Alleles.temp1[!(paste(Alleles.temp1, Allele_labels.temp1, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))]
        
        # Identify the new labels of the alleles
        Relabeled.Alleles.temp1 = Alleles[Alleles %in% Changed.Alleles.temp1]
        
        # For each changed label, modify the label in the genotype table
        for(clabel in names(Changed.Alleles.temp1)){
          
          rlabel = names(Relabeled.Alleles.temp1[Relabeled.Alleles.temp1 == Changed.Alleles.temp1[clabel]])
          
          gt_temp1[pos,] = gsub(paste0(clabel,':'), paste0(rlabel,'R:'), gt_temp1[pos,]) # The R: denotes that that label is been modify and avoids overwriting
          
        }
        
        # Once screening of all alleles is been completed, delete the R: identifier
        gt_temp1[pos,] = gsub('R:', ':', gt_temp1[pos,])
        
      }
      
      # If the alternative alleles  and their labels of Second data set are NOT contained in the merged data set
      if(sum(!(paste(Alleles.temp2, Allele_labels.temp2, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))) != 0){
        
        # Identify and select the alleles which lables have changed
        Changed.Alleles.temp2 = Alleles.temp2[!(paste(Alleles.temp2, Allele_labels.temp2, sep = ':') %in% paste(Alleles, Allele_labels, sep = ':'))]
        
        # Identify the new labels of the alleles
        Relabeled.Alleles.temp2 = Alleles[Alleles %in% Changed.Alleles.temp2]
        
        # For each changed label, modify the label in the genotype table
        for(clabel in names(Changed.Alleles.temp2)){
          
          rlabel = names(Relabeled.Alleles.temp2[Relabeled.Alleles.temp2 == Changed.Alleles.temp2[clabel]])
          
          gt_temp2[pos,] = gsub(paste0(clabel,':'), paste0(rlabel,'R:'), gt_temp2[pos,]) # The R: denotes that that label is been modify and avoids overwriting
          
        }
        
        # Once screening of all alleles is been completed, delete the R: identifier
        gt_temp2[pos,] = gsub('R:', ':', gt_temp2[pos,])
        
      }
      
      # Update Alleles in the Locus_info_merged_w_ALTdiscrepancies data.frame
      Locus_info_merged_w_ALTdiscrepancies[pos,][['Alleles']] = paste(paste(Alleles, Allele_labels,sep = ':'), collapse = ',')
      
    }
    
  }
  
  
  # Filter positions that were able to merge
  
  if(na.rm){
    Locus_info_merged_final = rbind(Locus_info_merged_good, Locus_info_merged_w_ALTdiscrepancies)
  }else{
    
    # adding missing positions in obj2
    Locus_info_merged_w_missing_temp1 = Locus_info_merged_w_missing[grepl('First', Locus_info_merged_w_missing$Alleles),]
    
    if(nrow(Locus_info_merged_w_missing_temp1) > 0){
      
      gt_temp2 = rbind(gt_temp2, matrix(
        NA, nrow = nrow(Locus_info_merged_w_missing_temp1), ncol = ncol(gt_temp2), dimnames = list(rownames(Locus_info_merged_w_missing_temp1),
                                                                                                   colnames(gt_temp2))
      ))
      
    }
    
    # adding missing positions in obj1
    
    Locus_info_merged_w_missing_temp2 = Locus_info_merged_w_missing[grepl('Second', Locus_info_merged_w_missing$Alleles),]
    
    if(nrow(Locus_info_merged_w_missing_temp2) > 0){
      
      gt_temp1 = rbind(gt_temp1, matrix(
        NA, nrow = nrow(Locus_info_merged_w_missing_temp2), ncol = ncol(gt_temp1), dimnames = list(rownames(Locus_info_merged_w_missing_temp2),
                                                                                                   colnames(gt_temp1))
      ))
      
    }
    
    
    Locus_info_merged_w_missing$Alleles = gsub(';Only (First|Second) data set', '', Locus_info_merged_w_missing$Alleles)
    Locus_info_merged_final = rbind(Locus_info_merged_good, Locus_info_merged_w_ALTdiscrepancies, Locus_info_merged_w_missing)
  }
  
  
  # Update Alternative Alleles (ALT)
  Locus_info_merged_final$ALT  = apply(Locus_info_merged_final, 1, function(ALT){
    gsub(':\\d+', '', gsub(paste0('^',paste(ALT['REF'], '0', sep = ':'), ','), '', ALT['Alleles']))
  })
  
  Locus_info_merged_final %<>% mutate(ALT = case_when(
    REF ==  gsub(':.+', '',Alleles) ~ NA,
    REF !=  gsub(':.+', '',Alleles) ~ ALT
  ))
  
  # Sort positions by position and chromosome
  Locus_info_merged_final = Locus_info_merged_final[order(Locus_info_merged_final$POS.temp1),]
  Locus_info_merged_final = Locus_info_merged_final[order(Locus_info_merged_final$CHROM.temp1),]
  
  # Select columns that will be returned
  Locus_info_merged_final = Locus_info_merged_final[,c('CHROM.temp1',
                                                       'POS.temp1',
                                                       'REF',
                                                       'ALT',
                                                       'Alleles')]
  
  colnames(Locus_info_merged_final) = c('CHROM',
                                        'POS',
                                        'REF',
                                        'ALT',
                                        'Alleles')
  
  # Filter and merge the genotype tables
  gt_final = cbind(gt_temp1[rownames(Locus_info_merged_final),], gt_temp2[rownames(Locus_info_merged_final),])
  
  
  colnames(gt_final) = c(colnames(gt_temp1), colnames(gt_temp2))
  
  # Merge the metadata
  merged_metadata = merge(obj1@metadata, obj2@metadata, by = 'Sample_id', all = T)
  
  rownames(merged_metadata) = merged_metadata$Sample_id
  
  final_merged_metadata = data.frame(merged_metadata[colnames(gt_final),])
  
  colnames(final_merged_metadata) = colnames(merged_metadata)
  
  rownames(final_merged_metadata) = merged_metadata$Sample_id
  
  # Create the merged rGenome object
  merged_rGenome = rGenome(gt = gt_final,
                           loci_table = Locus_info_merged_final,
                           metadata = final_merged_metadata)
  
  return(merged_rGenome)
  
}

# filter_loci_table----

filter_loci_table = function(loci_table = NULL, # table with loci information
                             vcf_object = NULL, # table with loci information in vcf format
                            filter_table, # Table of filtered segments or positions
                            by = 'segments'){
  
  
  if(!is.null(loci_table)){
    
    positions = NULL
    
    options(scipen=999)
    
    for(pos in 1:nrow(filter_table)){
      positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                     filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
      
    }
    
    options(scipen=0)
    
    loci_table = loci_table[rownames(loci_table) %in% positions,]
  
  }else if(!is.null(vcf_table)){
    
    rownames(vcf_object) = paste(vcf_object[['CHROM']], vcf_object[['POS']], sep = '_')
    
    positions = NULL
    
    options(scipen=999)
    
    for(pos in 1:nrow(filter_table)){
      positions = c(positions, paste(filter_table[pos,][['CHROM']],
                                     filter_table[pos,][['start']]:filter_table[pos,][['end']], sep = '_'))
      
    }
    
    options(scipen=0)
    
    vcf_object = vcf_object[rownames(vcf_object) %in% positions,]
    
    loci_table = vcf_object[,c(1,2,4,5)]
    
    gt = NULL
    
    for(w in 1:10){
      gt = rbind(gt, get_GTAD_matrix(vcf = vcf_object, w = w, n = 10, threshold = 5))
    }
    
    loci_table = cbind(loci_table, get_AC(gt, loci_table, w = 1, n = 1))
    
    loci_table$Type_of_polymorphism = get_type_of_polymorphism(loci_table, w = 1, n = 1)
    
  }
  
  return(loci_table)
  
}


# get_MHAP_fragment ----

get_MHAP_fragments = function(MHAPs = NULL,
                             upstream_buffer = 50,
                             downstream_buffer = 50,
                             path_to_reference_gff = 'genes.gff',
                             path_to_reference_genome = 'PvP01.v1.fasta',
                             loci_table = Locus_info
                             ){
  library(Biostrings)
  
  reference_gff = ape::read.gff(path_to_reference_gff)
  reference_genome = Biostrings::readDNAStringSet(path_to_reference_genome, format = 'fasta')
  
  filtered_loci_table = NULL
  target_sequences = NULL
  
  for(MHAP in 1:nrow(MHAPs)){
    
    temp_filtered_loci_table = filter_loci_table(loci_table = loci_table, filter_table = MHAPs[MHAP,])
    
    
    upstream_filter_table = data.frame(CHROM = MHAPs[MHAP,][['CHROM']],
                                       start = MHAPs[MHAP,][['start']] - upstream_buffer,
                                       end = min(temp_filtered_loci_table$POS - 1))
    
    temp_filtered_upstream_loci_table = filter_loci_table(loci_table = loci_table, filter_table = upstream_filter_table)
    
    
    downstream_filter_table = data.frame(CHROM = MHAPs[MHAP,][['CHROM']],
                                         start = max(temp_filtered_loci_table$POS) + 1,
                                         end = MHAPs[MHAP,][['end']] + downstream_buffer)
    
    temp_filtered_downstream_loci_table = filter_loci_table(loci_table = loci_table, filter_table = downstream_filter_table)
    
    
    temp_filtered_loci_table$polymorphism_location = 'Target sequence'
    if(nrow(temp_filtered_upstream_loci_table) > 0){
      temp_filtered_upstream_loci_table$polymorphism_location = 'Upstream region'
    }
    if(nrow(temp_filtered_downstream_loci_table)>0){
      temp_filtered_downstream_loci_table$polymorphism_location = 'Downstream region'
    }
    
    temp_filtered_loci_table = rbind(
      temp_filtered_upstream_loci_table,
      temp_filtered_loci_table,
      temp_filtered_downstream_loci_table
    )
    
    temp_filtered_loci_table = temp_filtered_loci_table[,c('CHROM', 'POS', 'Alleles', 'Allele_Counts','Type_of_polymorphism', 'polymorphism_location')]
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = as.character(subseq(reference_genome[grep(MHAPs[MHAP,][['CHROM']], names(reference_genome))],
                                                                      start = MHAPs[MHAP,][['start']] - upstream_buffer,
                                                                      end = MHAPs[MHAP,][['end']] + downstream_buffer))
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = unlist(str_split(target_sequences[[MHAPs[MHAP,][['Locus']]]], ''))
    
    for(polymorphism in 1:nrow(temp_filtered_loci_table)){
      
      if(temp_filtered_loci_table[polymorphism,][['Type_of_polymorphism']] == 'SNP'){
        
        target_sequences[[MHAPs[MHAP,][['Locus']]]][temp_filtered_loci_table[polymorphism,][['POS']] - (MHAPs[MHAP,][['start']] - upstream_buffer) + 1] = mergeIUPACLetters(
          paste(unlist(str_split(gsub(':\\d+','',temp_filtered_loci_table[polymorphism,][['Alleles']]),',')), collapse = ''))
        
        
      }else{
        target_sequences[[MHAPs[MHAP,][['Locus']]]][temp_filtered_loci_table[polymorphism,][['POS']] - (MHAPs[MHAP,][['start']] - upstream_buffer) + 1] = paste(
          rep('-',max(nchar(unlist(str_split(gsub(':\\d+','',temp_filtered_loci_table[polymorphism,][['Alleles']]),','))))), collapse = '')
        
      }
      
    }
    
    
    upstream_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][1:(min(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                       (MHAPs[MHAP,][['start']] - upstream_buffer))]
    
    targeted_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][(min(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                     (MHAPs[MHAP,][['start']] - upstream_buffer) + 1):(max(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                                                                         (MHAPs[MHAP,][['start']] - upstream_buffer) + 1)]
    
    downstream_region = target_sequences[[MHAPs[MHAP,][['Locus']]]][(max(temp_filtered_loci_table[temp_filtered_loci_table$polymorphism_location == 'Target sequence',][['POS']]) -
                                                                       (MHAPs[MHAP,][['start']] - upstream_buffer) + 2):length(target_sequences[[MHAPs[MHAP,][['Locus']]]])]
    
    
    target_sequences[[MHAPs[MHAP,][['Locus']]]] = paste(paste(upstream_region, collapse = ''),
                                                        '[',
                                                        paste(targeted_region, collapse = ''),
                                                        ']',
                                                        paste(downstream_region, collapse = ''), sep = '')
    
    
    temp_filtered_loci_table = data.frame(MHAPs[MHAP,c('Locus', 'start', 'end')], temp_filtered_loci_table)
    
    temp_filtered_loci_table$start = temp_filtered_loci_table$start - upstream_buffer
    temp_filtered_loci_table$end = temp_filtered_loci_table$end + downstream_buffer
    
    filtered_loci_table = rbind(filtered_loci_table,
                                temp_filtered_loci_table)
    
  }
  
  return(target_sequences)
  
}


# get_MHAPs_for_geneTarget ----

get_MHAPs_for_geneTarget = function(gene_id,
                                    upstream_buffer = 50,
                                    downstream_buffer = 50,
                                    path_to_reference_gff = 'genes.gff',
                                    path_to_reference_genome = 'PvP01.v1.fasta',
                                    loci_table = Locus_info,
                                    max_length){
  
  library(Biostrings)
  
  reference_gff = ape::read.gff(path_to_reference_gff)
  reference_genome = Biostrings::readDNAStringSet(path_to_reference_genome, format = 'fasta')
  
  selected_CDS = reference_gff[grepl(gene_id,reference_gff$attributes)&reference_gff$type == 'CDS',]
  
  selected_CDS = selected_CDS[,c('seqid', 'start', 'end')]
  names(selected_CDS) = c('CHROM', 'start', 'end')
  
  selelcted_CDS_loci_table = filter_loci_table(loci_table = loci_table,
                                               filter_table = selected_CDS,
                                               by = 'segments')
  
  selelcted_CDS_loci_table$dist = c(selelcted_CDS_loci_table$POS[-1] -
                                      selelcted_CDS_loci_table$POS[-length(selelcted_CDS_loci_table$POS)],
                                    Inf)
  
  MHAPs = NULL
  pos = 1
  start_pos = 1
  end_pos = 1
  
  temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                         start = selelcted_CDS_loci_table$POS[pos],
                         end = selelcted_CDS_loci_table$POS[pos])
  
  while(pos <= length(selelcted_CDS_loci_table$POS)){
    
    if(temp_MHAP$end - temp_MHAP$start + selelcted_CDS_loci_table$dist[pos] <= max_length){
      pos = pos + 1
      end_pos = pos
      
      temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                             start = selelcted_CDS_loci_table$POS[start_pos],
                             end = selelcted_CDS_loci_table$POS[end_pos])
      
    }else{
      pos = pos + 1
      start_pos = pos
      end_pos = pos
      
      MHAPs = rbind(MHAPs, temp_MHAP)
      
      temp_MHAP = data.frame(CHROM = selelcted_CDS_loci_table$CHROM[pos],
                             start = selelcted_CDS_loci_table$POS[start_pos],
                             end = selelcted_CDS_loci_table$POS[end_pos])
    }
    
  }
  
  MHAPs = data.frame(Locus = paste(MHAPs$CHROM, MHAPs$start, MHAPs$end, sep = '_'), MHAPs)
  
  target_sequences = get_MHAP_fragments(MHAPs = MHAPs)
  
  return(target_sequences)
  
}


# find_sequence----

find_DNAsequence = function(sequences = NULL,
                         reference_genome = NULL 
                         ){
  
  sequence_location = NULL
  
  for(DNAsequence in sequences){
    
    sequence_length = nchar(DNAsequence)
    
    print(paste0('searching DNA sequence ', DNAsequence))
    
    #grep(DNAsequence, as.character(reference_genome)) 
    
    chrom = grep(DNAsequence, as.character(reference_genome)) 
    
    if(length(chrom) == 0){
      DNAsequence = as.character(reverseComplement(DNAString(DNAsequence)))
      chrom = grep(DNAsequence, as.character(reference_genome)) 
    }
    
    chrom = chrom[1] # CORRECT THIS LINE
    
    print(paste0('searching DNA sequence in chrom ', chrom))
    
    temp_chrom = unlist(str_split(as.character(reference_genome[[chrom]]),''))
    
    w = 100000
    
    low = 1
    high = nchar(as.character(reference_genome[[chrom]]))
    
    while(w > sequence_length){
      
      s = c(seq(low, high, w), high + 1)
      
      temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = s[-length(s)],
                               end = s[-1] -1)
      
      seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
      
      i = grep(DNAsequence, seqs)
      
      i = i[length(i)] # CORRECT THIS LINE
      
      while(length(i) == 0){
        s = s + 1
        
        temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = s[-length(s)],
                                 end = s[-1] -1)
        
        seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
        
        i = grep(DNAsequence, seqs)
        
        
      }
      
      low = s[i]
      high = s[i+1]
      
      w = w/10
      
    }
    
    temp_search = data.frame(CHROM =names(reference_genome)[chrom], start = low:(high - (sequence_length - 1)),
                             end = (low + sequence_length -1):high)
    
    seqs = apply(temp_search, 1, function(x){paste(temp_chrom[seq(x['start'], x['end'],1)], collapse = '')})
    
    sequence_location = rbind(sequence_location, data.frame(DNAsequence = DNAsequence, temp_search[grep(DNAsequence, seqs)[1]
                                                                                                   ,]))
    
  }
  
  rownames(sequence_location) = 1:nrow(sequence_location)
  
  return(sequence_location)
  
}

# get_haplotypes_respect_to_reference----

get_haplotypes_respect_to_reference = function(obj,
                                               gene_ids = NULL,
                                               gene_labels = NULL,
                                               gff_file = NULL,
                                               fasta_file = NULL,
                                               monoclonals = NULL,
                                               polyclonals = NULL,
                                               variables = NULL,
                                               mon_poly_ratio = 1,
                                               aa_format = 'compact', #
                                               na_pos_rm = TRUE){
  library(ape)
  library(Biostrings)
  library(msa)
  
  # Call reference genome and its corresponding anotation in the gff file
  reference_gff = read.gff(gff_file)
  reference_genome = readDNAStringSet(fasta_file)
  
  # Get gene ids, names, and descriptions 
  gene_names = data.frame(gene_id = gsub(';.+$','',gsub('ID=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())),
                          gene_name = gsub(';.+$','',gsub('^.+Name=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())),
                          gene_description = gsub(';.+$','',gsub('^.+description=','',reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes),  grepl('gene', type)) %>% select(attributes) %>% unlist())))
  
  # Split gene coordinates by coding sequences (CDS)
  genes_gff = reference_gff %>% filter(grepl(paste(gene_ids, collapse = '|'), attributes) & type == 'CDS')
  
  # Get gene ids
  genes_gff %<>% mutate(gene_id = gsub(';.+$|\\..+$','',gsub('^.+gene_id=|^ID=','',attributes)))
  
  # Merge CDS coordinates with gene names and descriptions by gene_id
  genes_gff = merge(genes_gff, gene_names, by = 'gene_id', all.x = T)
  
  
  
  Sample_id = obj@metadata$Sample_id
  
  monoclonals_ids = Sample_id[Sample_id %in% monoclonals]
  polyclonals_ids = Sample_id[Sample_id %in% polyclonals]
  
  if(mon_poly_ratio > 1 & length(polyclonals_ids) > 0){
    
    haplotypes_ids = c(paste0(rep(monoclonals_ids, mon_poly_ratio), '_C', rep(1:mon_poly_ratio, each = length(monoclonals_ids))),
                       paste0(polyclonals_ids, '_C1'),
                       paste0(polyclonals_ids, '_C2'))
    
  }else if(mon_poly_ratio == 1 & length(polyclonals_ids) > 0){
    
    haplotypes_ids = c(monoclonals_ids,
                       paste0(polyclonals_ids, '_C1'),
                       paste0(polyclonals_ids, '_C2'))
  }else if(length(polyclonals_ids) == 0){
    
    haplotypes_ids = monoclonals_ids
      
    }
  
  haplotypes = matrix(NA, ncol = length(gene_ids), nrow = length(haplotypes_ids),
                      dimnames = list(haplotypes_ids, gene_ids))
  
  
  gDNA_haplotypes = matrix(NA, ncol = length(gene_ids), nrow = length(haplotypes_ids),
                           dimnames = list(haplotypes_ids, gene_ids))
  
  
  # For each gene:
  
  for(gene in unique(genes_gff$gene_id)){
    
    # Filter polymorphic positions located in the gene of interest (GOI)
    ## Get the cds coordinates of the goi
    cds_gff = genes_gff[genes_gff$gene_id == gene,]
    ## Get the chromosome sequence where the goi is located
    ref_seq = reference_genome[[grep(unique(cds_gff$seqid), names(reference_genome))]]
    
    # Generate a vector with all DNA coordinates in the chromosome
    dna_regions = NULL
    gene_seq = NULL
    
    for(cds in 1:nrow(cds_gff)){
      
      # Get nucleotide coordinates
      dna_regions = c(dna_regions,
                      paste(cds_gff[cds,][['seqid']],
                            cds_gff[cds,][['start']]:cds_gff[cds,][['end']],
                            sep = '_'))
      
      # Get DNA sequence
      gene_seq = paste0(gene_seq, as.character(subseq(ref_seq, start = cds_gff[cds,][['start']], end = cds_gff[cds,][['end']])))
      
    }
    
    names(gene_seq) = 'reference'
    
    if(sum(rownames(obj@loci_table) %in% dna_regions) > 0){
      
      # Filter the rGenome object based on the gene coordinates
      
      gene_obj = filter_loci(obj, v = rownames(obj@loci_table) %in% dna_regions)
      
      ## Remove read abundance
      
      gt = handle_ploidy(gt = gene_obj@gt, monoclonals = monoclonals, polyclonals = polyclonals)
      
      
      if(mon_poly_ratio > 1){
        
        gt_monoclonals = NULL
        
        for(i in 1:as.integer(mon_poly_ratio)){
          
          gt_monoclonals_temp = matrix(gt[, monoclonals], nrow = nrow(gt), dimnames = list(rownames(gt), paste0(monoclonals, '_C', i)))
          
          gt_monoclonals = cbind(gt_monoclonals, gt_monoclonals_temp)
          
        }
        
        gt = cbind(gt_monoclonals, matrix(gt[, grepl(paste(polyclonals, collapse = '|'), colnames(gt))],
                                          nrow = nrow(gt),
                                          dimnames = list(rownames(gt),
                                                          c(paste0(polyclonals_ids, '_C1'),
                                                            paste0(polyclonals_ids, '_C2')))))
        
      }
      
      # Tarnsform haplotype codes to nucleotides
      nuc_gt = NULL
      for(pos in 1:nrow(gt)){
        temp_nuc_gt = gt[pos,]
        allele_codes = unique(unlist(strsplit(gt[pos,], '/')))
        allele_codes = min(as.integer(allele_codes), na.rm = T):max(as.integer(allele_codes), na.rm = T)
        alleles = unlist(strsplit(gsub(':\\d+','',gene_obj@loci_table[pos,"Alleles"]), ','))
        
        for(allele in allele_codes){
          temp_nuc_gt = gsub(allele_codes[which(allele_codes == allele)],
                             alleles[which(allele_codes == allele)],
                             temp_nuc_gt)
        }
        
        temp_nuc_gt = gsub('\\*', '',temp_nuc_gt)
        nuc_gt = rbind(nuc_gt, temp_nuc_gt)
      }
      rownames(nuc_gt) = rownames(gt)
      
      #nuc_gt[grepl('NA', nuc_gt)] = NA # remove this line later
      
      # Get DNA haplotypes
      
      gDNA_haplotypes[,gene] = apply(nuc_gt, 2, function(sample){
        
        sample[is.na(sample)] = '?'
        
        pos = gene_obj@loci_table$POS
        ref = gene_obj@loci_table$REF
        
        changes = (ref != sample)
        
        ifelse(sum(changes) == 0, 'g.(=)',paste(paste0('g.', pos[changes], ref[changes], '>', sample[changes]), collapse = ' '))
        
        })
      
      
      
      ref_polymorphism_length = nchar(gene_obj@loci_table$REF)
      
      # Define the coordinates of the polymorphims in the gene sequence
      polymorphic_positions = which(dna_regions %in% rownames(nuc_gt))
      
      # For each sample for each codon define the aminoacid changes in the gene
      
      sample_seqs = NULL
      
      for(sample in 1:ncol(nuc_gt)){
        
        sample_seq = gene_seq
        
        for(pos in 1:length(polymorphic_positions)){
          
          nucleotide = nuc_gt[pos,sample]
          
          if(!is.na(nucleotide)){
            
            nucleotide_to_substitute = substr(sample_seq,
                   polymorphic_positions[pos],
                   polymorphic_positions[pos] + ref_polymorphism_length[pos] - 1)
            
            if(nucleotide_to_substitute != '-'){
              
              substr(sample_seq,
                     polymorphic_positions[pos],
                     polymorphic_positions[pos] + ref_polymorphism_length[pos] - 1) = nuc_gt[pos,sample]
              
            }
            
          }
          
        }
        
        sample_seq = gsub('(-+|\\*+)', '', sample_seq)
        
        sample_seqs = c(sample_seqs, sample_seq)
        
      }
      
      names(sample_seqs) = colnames(nuc_gt)
      
      if(cds_gff$strand[1] == '+'){
        
        aa_seqs = AAStringSet(c(as.character(Biostrings::translate(Biostrings::DNAString(gene_seq))),
                                as.character(Biostrings::translate(Biostrings::DNAStringSet(sample_seqs)))))
        
      }else if(cds_gff$strand[1] == '-'){
        
        aa_seqs = AAStringSet(c(as.character(Biostrings::translate(Biostrings::reverseComplement(Biostrings::DNAString(gene_seq)))),
                                as.character(Biostrings::translate(Biostrings::reverseComplement(Biostrings::DNAStringSet(sample_seqs))))))
        
      }
      
      names(aa_seqs) = c('reference', names(sample_seqs))
      
      aa_alignment = msa(aa_seqs, method = 'ClustalOmega')

      
      aa_alignment_matrix = matrix(unlist(strsplit(as.character(aa_alignment@unmasked), '')), ncol = length(aa_alignment@unmasked),
                                   nrow = nchar(as.character(aa_alignment@unmasked[['reference']])),
                                   dimnames = list(
                                     1:nchar(as.character(aa_alignment@unmasked[['reference']])),
                                     names(aa_alignment@unmasked)
                                   )
      )
      
      
      
      if(cds_gff$strand[1] == '+'){
        
        aa_expected_polymorphic_positions = unique(ceiling(polymorphic_positions/3))
        
      }else if(cds_gff$strand[1] == '-'){
        
        aa_expected_polymorphic_positions = sort(unique(ceiling((nchar(gene_seq) - polymorphic_positions + 1)/3)))
        
      }
      
      
      aa_haplotypes = NULL
      
      for(sample in names(sample_seqs)){
        
        ref_aa_aligned = aa_alignment_matrix[,'reference']
        
        samp_aa_aligned = aa_alignment_matrix[,sample]
        
        haplotype = paste0(ref_aa_aligned[aa_expected_polymorphic_positions], aa_expected_polymorphic_positions, samp_aa_aligned[aa_expected_polymorphic_positions])
        
        
        for(pos in 1:length(polymorphic_positions)){
          
          if(cds_gff$strand[1] == '+'){
            
            codon = ceiling(polymorphic_positions[pos]/3)
            
          }else{
            
            codon = ceiling((nchar(gene_seq) - polymorphic_positions[pos] + 1)/3)
            
          }
          
          if(sum(grepl(codon, haplotype)) > 0){
            
            if(is.na(nuc_gt[pos, sample])){
              
              haplotype[grepl(codon, haplotype)] = paste0(gsub('\\w$', '', haplotype[grepl(codon, haplotype)]), '?')
              
            }
            
          }
          
        }
        
        
        if(aa_format == 'compact'){
          
          sample_polymorphic_positions = sapply(haplotype, function(pos){
            gsub('\\d+\\w+$', '', pos) != gsub('^\\w\\d+', '', pos)})
          
          if(sum(sample_polymorphic_positions) == 0){
            
            haplotype = 'p.(=)'
            
          }else if(sum(sample_polymorphic_positions) > 0 & 
                   sum(sample_polymorphic_positions) < length(sample_polymorphic_positions)){
            
            haplotype = haplotype[sample_polymorphic_positions]
            
          }
          
          sample_missing_positions = sapply(haplotype, function(pos){
            grepl('\\?', pos)})
          
          if(sum(sample_missing_positions) == length(sample_polymorphic_positions)){
            
            haplotype = 'p.(?)'
            
          }else if(sum(sample_missing_positions) > 0 & 
                   sum(sample_missing_positions) < length(haplotype) &
                   na_pos_rm){
            
            haplotype = haplotype[!sample_missing_positions]
            
            haplotype = c(haplotype, 'Partial')
            
          }else if(sum(sample_missing_positions) > 0 & 
                   sum(sample_missing_positions) == length(haplotype) &
                   sum(sample_missing_positions) < length(sample_polymorphic_positions)){
            
            haplotype = 'p.(=?)'
            
          }
          
        }
        
        
        haplotype = paste(haplotype, collapse = ' ')
        
        names(haplotype) = sample
        
        aa_haplotypes = c(aa_haplotypes, haplotype)
        
      }
      
      haplotypes[rownames(haplotypes) %in% names(aa_haplotypes), gene] = aa_haplotypes
      
    }else{
      
      haplotypes[, gene] = 'p.(=)'
      
    }
    
  }
  
  haplo_freqs = data.frame(Sample_id = rownames(haplotypes), haplotypes)
  
  metadata = obj@metadata[,c('Sample_id', variables)]
  
  if(length(variables) == 1){
    colnames(metadata) = c('Sample_id', 'Var1')
  }else if(length(variables) == 2){
    colnames(metadata) = c('Sample_id', 'Var1', 'Var2')
  }

  if(!is.null(polyclonals)){
    
    metadata_monoclonals = metadata[metadata$Sample_id %in% monoclonals,]
    metadata_polyclonals = metadata[metadata$Sample_id %in% polyclonals,]
    
    metadata_polyclonals1 = metadata_polyclonals
    metadata_polyclonals1$Sample_id = paste0(metadata_polyclonals1$Sample_id, '_C1')
    metadata_polyclonals2 = metadata_polyclonals
    metadata_polyclonals2$Sample_id = paste0(metadata_polyclonals2$Sample_id, '_C2')
    
    
    if(mon_poly_ratio > 1){
      
      metadata = NULL
      
      for(i in 1:mon_poly_ratio){
        
        metadata_monoclonals_temp = metadata_monoclonals
        metadata_monoclonals_temp$Sample_id = paste0(metadata_monoclonals_temp$Sample_id, '_C', i)
        
        metadata = rbind(metadata, metadata_monoclonals_temp)
      }
      
    }else{
      
      metadata = metadata_monoclonals
      
    }
    
    metadata = rbind(metadata, metadata_polyclonals1, metadata_polyclonals2)
    
  }
  
  
  rownames(metadata) = metadata$Sample_id
  
  haplo_freqs = merge(haplo_freqs, metadata, by = 'Sample_id')
  
  haplo_freqs %<>% pivot_longer(cols = colnames(haplotypes), names_to = 'Gene_id', values_to = 'Haplotype')
  
  gene_labels = data.frame(gene_ids, gene_labels)
  
  haplo_freqs$gene_label = NA
  
  for(gene in gene_labels$gene_ids){
    
    haplo_freqs[haplo_freqs$Gene_id == gene,][['gene_label']] = gene_labels[gene_labels$gene_ids == gene,][['gene_labels']]
    
  }
  
  if(length(variables) == 1){
    haplo_freqs %<>% 
      mutate(Haplotype = paste(gene_label, Haplotype, sep = ':'))%>%
      dplyr::summarise(Sample_id = Sample_id,
                Haplotype = Haplotype,
                Sample_size = n(),
                .by = c(Var1, gene_label)
      ) %>%
      dplyr::summarise(Count = n(),
                       Freq = binconf(n(), unique(Sample_size))[1],
                       lower = binconf(n(), unique(Sample_size))[2],
                       upper = binconf(n(), unique(Sample_size))[3],
                .by= c(Var1, gene_label, Haplotype))
    
    haplo_freq_plot = haplo_freqs%>%
      ggplot(aes(y = Freq, x = Var1, color = Haplotype, group = Haplotype))+
      geom_point()+
      geom_errorbar(aes(ymin = lower, ymax = upper), alpha = .5, width = .2)+
      facet_grid(.~gene_label)+
      labs(x = variables[1])+
      theme(legend.position = 'bottom')
  }else if(length(variables) == 2){
    haplo_freqs %<>% 
      mutate(Haplotype = paste(gene_label, Haplotype, sep = ':'))%>%
      dplyr::summarise(Sample_id = Sample_id,
                Haplotype = Haplotype,
                Sample_size = n(),
                .by = c(Var1, Var2, gene_label)
      ) %>%
      dplyr::summarise(Count = n(),
                       Freq = binconf(n(), unique(Sample_size))[1],
                       lower = binconf(n(), unique(Sample_size))[2],
                       upper = binconf(n(), unique(Sample_size))[3],
                .by= c(Var1, Var2, gene_label, Haplotype))
    
    haplo_freq_plot = haplo_freqs %>%
      ggplot(aes(y = Freq, x = Var2, color = Haplotype, group = Haplotype))+
      geom_point()+
      geom_errorbar(aes(ymin = lower, ymax = upper), alpha = .5, width = .2)+
      geom_line()+
      facet_grid(Var1~gene_label)+
      labs(x = variables[2])+
      theme(legend.position = 'bottom')
  }
  
  # Collapse polyclonal gDNA haplotypes per sample
  
  gDNA_haplotypes_collapsed = matrix(NA, 
                                     nrow = length(unique(gsub('_C\\d+$', '', haplotypes_ids))),
                                     ncol = length(colnames(gDNA_haplotypes)),
                                     dimnames = list(
                                       unique(gsub('_C\\d+$', '', haplotypes_ids)),
                                              colnames(gDNA_haplotypes)
                                              
                                       )
                                     )
  
  for(gene in colnames(gDNA_haplotypes)){
    
    for(sample in unique(gsub('_C\\d+$', '', haplotypes_ids))){
      
      genotype = unique(gDNA_haplotypes[grepl(sample, rownames(gDNA_haplotypes)), gene])
      
      if(length(genotype) > 1){
        
        genotype = strsplit(genotype, ' ')
        
        positions = sort(unique(gsub('>\\w+$', '', unique(unlist(genotype)))))
        
        for(position in positions){
          
          genotype_temp = NULL
          
          for(haplo in 1:length(genotype)){
            
            haplo_temp = genotype[[haplo]]
            
            haplo_temp = haplo_temp[grepl(position, haplo_temp)]
            
            if(length(haplo_temp) > 0){
              
              if(is.null(genotype_temp)){
                
                genotype_temp = haplo_temp
                
              }else{
                
                if(!grepl(gsub('g.\\d+\\w+>', '', haplo_temp), gsub('g.\\d+\\w+>', '', genotype_temp))){
                  
                  genotype_temp = paste(genotype_temp, gsub('g.\\d+\\w+>', '', haplo_temp), sep = '|')
                  
                }
                
              }
              
              
            }else{
              
              haplo_temp = paste(position, gsub('^g.\\d+|>\\w+$', '',position), sep =  '>')
              
              if(is.null(genotype_temp)){
                
                genotype_temp = haplo_temp
                
              }else{
                
                if(!grepl(gsub('g.\\d+\\w+>', '', haplo_temp), gsub('g.\\d+\\w+>', '', genotype_temp))){
                  
                  genotype_temp = paste(genotype_temp, gsub('g.\\d+\\w+>', '', haplo_temp), sep = '|')
                  
                }
                
              }
              
              
            }
            
            
          }
          
          if(!is.na(gDNA_haplotypes_collapsed[sample,gene])){
            
            gDNA_haplotypes_collapsed[sample,gene] = paste(gDNA_haplotypes_collapsed[sample,gene], genotype_temp, sep = ' ')
            
          }else{
            
            gDNA_haplotypes_collapsed[sample,gene] = genotype_temp
            
          }
          
          
          
        }
        
        
          
        }else{
        
        gDNA_haplotypes_collapsed[sample,gene] = genotype
        
      }
      
    }
    
  }
  
  gDNA_haplotypes_collapsed = data.frame(Sample_id = rownames(gDNA_haplotypes_collapsed),
                                         gDNA_haplotypes_collapsed)
  
  colnames(gDNA_haplotypes_collapsed) = c('Sample_id', paste(gene_labels$gene_labels,
                                                           gene_labels$gene_ids,
                                                           sep = ': '))
  
  gDNA_haplotypes_collapsed %<>% pivot_longer(cols = all_of(colnames(gDNA_haplotypes_collapsed)[-1]),
                                            names_to = 'Gene',
                                            values_to = 'Haplotype')
  
  
  
  # Collapse polyclonal aa haplotypes per sample
  aa_haplotypes_collapsed = matrix(NA, 
                                     nrow = length(unique(gsub('_C\\d+$', '', haplotypes_ids))),
                                     ncol = length(colnames(haplotypes)),
                                     dimnames = list(
                                       unique(gsub('_C\\d+$', '', haplotypes_ids)),
                                       colnames(haplotypes)
                                       
                                     )
  )
  
  for(gene in colnames(haplotypes)){
    
    for(sample in unique(gsub('_C\\d+$', '', haplotypes_ids))){
      
      genotype = unique(haplotypes[grepl(sample, rownames(haplotypes)), gene])
      
      if(length(genotype) > 1){
        
        genotype = strsplit(genotype, ' ')
        
        positions = unique(gsub('([A-Z]|\\*)+$', '', unique(unlist(genotype))))
        
        for(position in positions){
          
          genotype_temp = NULL
          
          for(haplo in 1:length(genotype)){
            
            haplo_temp = genotype[[haplo]]
            
            haplo_temp = haplo_temp[grepl(position, haplo_temp)]
            
            if(length(haplo_temp) > 0){
              
              if(is.null(genotype_temp)){
                
                genotype_temp = haplo_temp
                
              }else{
                
                if(!grepl(gsub('^([A-Z]+|\\*)\\d+', '', haplo_temp), gsub('^([A-Z]+|\\*)\\d+', '', genotype_temp))){
                  
                  genotype_temp = paste(genotype_temp, gsub('^([A-Z]+|\\*)\\d+', '', haplo_temp), sep = '|')
                  
                }
                
              }
              
              
            }else{
              
              haplo_temp = paste(position, gsub('\\d+([A-Z]+\\*)$', '',position), sep =  '>')
              
              if(is.null(genotype_temp)){
                
                genotype_temp = haplo_temp
                
              }else{
                
                if(!grepl(gsub('^([A-Z]+|\\*)\\d+', '', haplo_temp), gsub('^([A-Z]+|\\*)\\d+', '', genotype_temp))){
                  
                  genotype_temp = paste(genotype_temp, gsub('^([A-Z]+|\\*)\\d+', '', haplo_temp), sep = '|')
                  
                }
                
              }
              
              
            }
            
            
          }
          
          if(!is.na(aa_haplotypes_collapsed[sample,gene])){
            
            aa_haplotypes_collapsed[sample,gene] = paste(aa_haplotypes_collapsed[sample,gene], genotype_temp, sep = ' ')
            
          }else{
            
            aa_haplotypes_collapsed[sample,gene] = genotype_temp
            
          }
          
          
          
        }
        
        
        
      }else{
        
        aa_haplotypes_collapsed[sample,gene] = genotype
        
      }
      
    }
    
  }
  
  
  aa_haplotypes_collapsed = data.frame(Sample_id = rownames(aa_haplotypes_collapsed),
                                       aa_haplotypes_collapsed)
  
  colnames(aa_haplotypes_collapsed) = c('Sample_id', paste(gene_labels$gene_labels,
                                            gene_labels$gene_ids,
                                            sep = ': '))
  
  aa_haplotypes_collapsed %<>% pivot_longer(cols = all_of(colnames(aa_haplotypes_collapsed)[-1]),
                                            names_to = 'Gene',
                                            values_to = 'Haplotype')
  
  return(list(aa_haplotypes = aa_haplotypes_collapsed,
              gDNA_haplotypes = gDNA_haplotypes_collapsed,
              haplo_freqs = haplo_freqs,
              haplo_freq_plot = haplo_freq_plot))
  
}


write_rGenome = function(rGenome_object, format = c('excel', 'tsv', 'json'), 
                         name = 'wb.xlsx',
                         sep = '\t'){
  
  if(format == 'excel'){
    
    if(file.exists(name)){
      system(paste0('rm ', name))
    }
    
    excel_wb = loadWorkbook(name, create = T)
    
    for(temp_slot in c('gt', 
                       'metadata', 
                       'loci_table')){
      
      if(temp_slot %in% c('gt', 'loci_table')){
        
        temp_sheet = data.frame(Position_id = rownames(slot(rGenome_object, temp_slot)),
                                as.data.frame(slot(rGenome_object, temp_slot)))
        
      }else{
        
        if(!is.null(slot(rGenome_object, temp_slot))){
          temp_sheet = as.data.frame(slot(rGenome_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        createSheet(excel_wb, name = temp_slot)
        
        writeWorksheet(excel_wb,
                       temp_sheet,
                       sheet = temp_slot,
                       header = T)
      }
      
      
    }
    
    saveWorkbook(excel_wb)
    
  }else if(format == 'tsv'){
    
    if(file.exists(name)){
      system(paste0('rm -r ', name))
    }
    
    system(paste0('mkdir ', name))
    
    for(temp_slot in c('gt', 
                       'loci_table',
                       'metadata')){
      
      
      
      if(temp_slot %in% c('gt', 'loci_table')){
        
        temp_sheet = data.frame(Position_id = rownames(slot(rGenome_object, temp_slot)),
                                as.data.frame(slot(rGenome_object, temp_slot)))
        
      }else{
        
        if(!is.null(slot(rGenome_object, temp_slot))){
          temp_sheet = as.data.frame(slot(rGenome_object, temp_slot))
        }else{
          temp_sheet = NULL
        }
        
        
      }
      
      if(!is.null(temp_sheet)){
        
        write.table(temp_sheet, paste0(file.path(name, temp_slot), '.tsv'), 
                    quote = F, 
                    row.names = F,
                    sep = sep)
        
      }
      
      
    }
    
    
  }else if(format == 'json'){
    # In development
    
  }
  
}

## read_rGenome ----
read_rGenome = function(file = NULL, format = 'excel', sep = '\t'){
  
  rGenome_object = rGenome()
  
  if(format == 'excel'){
    
    temp_wb = loadWorkbook(file)
    
    for(sheet in getSheets(temp_wb)){
      if(sheet == 'gt'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table'){
        
        temp_sheet = readWorksheet(temp_wb, sheet = sheet)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, sheet, check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'tsv'){
    
    for(sheet in list.files(file)){
      if(sheet == 'gt.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = as.matrix(temp_sheet[,-1])
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'metadata.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }else if(sheet == 'loci_table.tsv'){
        
        temp_sheet = read.table(file.path(file, sheet), header = T, sep = sep)
        temp_sheet_rownames = temp_sheet[,1]
        temp_sheet = temp_sheet[,-1]
        rownames(temp_sheet) = temp_sheet_rownames
        
        slot(rGenome_object, gsub('.tsv','',sheet), check = TRUE) = temp_sheet
        
      }
    }
    
  }else if(format == 'json'){
    # In development
  }
  
  return(rGenome_object)
  
}




