#!/bin/bash
source /broad/software/scripts/useuse
use R-4.1

function json_extract() {
  # Citation: https://stackoverflow.com/questions/1955505/parsing-json-with-unix-tools
  local key=$1
  local json=$2
  local string_regex='"([^"\]|\.)*"'
  local string_logial='([A-Za-z]+)*'
  local number_regex='-?(0|[1-9][0-9]*)(\.[0-9]+)?([eE][+-]?[0-9]+)?'
  local time_regex='[0-9][0-9]+:(([0-5][0-9])|60):(([0-5][0-9])|60)'
  local value_regex="${string_regex}|${number_regex}|${time_regex}|${string_logial}"
  local pair_regex="\"${key}\"[[:space:]]*:[[:space:]]*(${value_regex})"

  if [[ ${json} =~ ${pair_regex} ]]
  then
	echo $(sed 's/^"\|"$//g' <<< "${BASH_REMATCH[1]}")
  else
	return 1
  fi
	}
	
# Extract variables
json=$1

wd=$(json_extract wd "$(cat ${json})")
wd=${wd%\"}
wd=${wd#\"}

echo "wd: "${wd}

nTasks=$(json_extract nTasks "$(cat ${json})")

echo "nTasks: "${nTasks}

if [[ ${nTasks} > 1 ]]
  then
  merge_rgenome='true';
fi

echo "merge_rgenome: "${merge_rgenome}

fd=$(json_extract fd "$(cat ${json})")
fd=${fd%\"}
fd=${fd#\"}

echo "fd: "${fd}

output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

echo "output: "${output}

remove_filtered_all=$(json_extract remove_filtered_all "$(cat ${json})")
include_from_gff=$(json_extract include_from_gff "$(cat ${json})")
include_from_gff=${include_from_gff%\"}
include_from_gff=${include_from_gff#\"}
exclude_from_gff=$(json_extract exclude_from_gff "$(cat ${json})")
exclude_from_gff=${exclude_from_gff%\"}
exclude_from_gff=${exclude_from_gff#\"}

bed=$(json_extract bed "$(cat ${json})")
bed=${bed%\"}
bed=${bed#\"}

exclude_bed=$(json_extract exclude_bed "$(cat ${json})")
exclude_bed=${exclude_bed%\"}
exclude_bed=${exclude_bed#\"}

keep_regexp=$(json_extract keep_regexp "$(cat ${json})")
non_ref_ac_any=$(json_extract non_ref_ac_any "$(cat ${json})")

# Check if the prefiltering was requested
vcf_file=${output}

if [[ ${non_ref_ac_any} != NaN ]]
then
  vcf_file+="_NonREFacAny"${non_ref_ac_any};
fi

if [[ ${remove_filtered_all} == true ]]
  then
    vcf_file+='_filtered';
fi

if [[ ${include_from_gff} != NaN ]]
then
  vcf_file+="_PATTERN"${include_from_gff}"SELECTED";
fi

if [[ ${exclude_from_gff} != NaN ]]
then
  vcf_file+="_PATTERN"${exclude_from_gff}"REMOVED";
fi

if [[ ${bed} != NaN ]]
then
  vcf_file+="_SelectedFromBed";
fi

if [[ ${exclude_bed} != NaN ]]
then
  vcf_file+="_RemovedFromBed";
fi

if [[ ${keep_regexp} != NaN ]]
then
  vcf_file+="_SelectedSamplesOnly";
fi

if [[ ${remove_filtered_all} == true || ${include_from_gff} != NaN  || ${exclude_from_gff} != NaN || ${bed} != NaN  || ${exclude_bed} != NaN || ${keep_regexp} != NaN || ${non_ref_ac_any} != NaN ]]
then
  vcf_file+=".recode.vcf";
else
  vcf_file=$(json_extract vcf_file "$(cat ${json})");
  vcf_file=${vcf_file%\"};
  vcf_file=${vcf_file#\"};
fi

echo ${vcf_file}

ref_gff_file=$(json_extract ref_gff_file "$(cat ${json})")  
ref_gff_file=${ref_gff_file%\"}
ref_gff_file=${ref_gff_file#\"}

echo "ref_gff_file: "${ref_gff_file}

do_vcf2rGenome='false'

echo "do_vcf2rGenome: "${do_vcf2rGenome}

nchunks=$(json_extract nchunks "$(cat ${json})")

echo "nchunks: "${nchunks}

pop=$(json_extract pop "$(cat ${json})")
pop=${pop%\"}
pop=${pop#\"}

echo "pop: "${pop}

post_filtering=$(json_extract post_filtering "$(cat ${json})")

echo "post_filtering: "${post_filtering}

type_of_polymorphism_to_remove=$(json_extract type_of_polymorphism_to_remove "$(cat ${json})")
type_of_polymorphism_to_remove=${type_of_polymorphism_to_remove%\"}
type_of_polymorphism_to_remove=${type_of_polymorphism_to_remove#\"}

echo "type_of_polymorphism_to_remove: "${type_of_polymorphism_to_remove}

sample_ampl_rate=$(json_extract sample_ampl_rate "$(cat ${json})")

echo "sample_ampl_rate: "${sample_ampl_rate}

locus_ampl_rate=$(json_extract sample_ampl_rate "$(cat ${json})")

echo "locus_ampl_rate: "${sample_ampl_rate}

ObsHet_quantile=$(json_extract ObsHet_quantile "$(cat ${json})")

echo "ObsHet_quantile: "${ObsHet_quantile}

SNP_density_quantile=$(json_extract SNP_density_quantile "$(cat ${json})")

echo "SNP_density_quantile: "${SNP_density_quantile}


mask_formula=$(json_extract mask_formula "$(cat ${json})")
echo "mask_formula: "${mask_formula}

alignment_filter=$(json_extract alignment_filter "$(cat ${json})")
alignment_filter=${alignment_filter%\"}
alignment_filter=${alignment_filter#\"}

echo "alignment_filter: "${alignment_filter}


cd ${wd}

# Run post filtering
Rscript ${fd}/vcf2rGenome_filtering.R \
  -wd ${wd} \
  -fd ${fd} \
  -o ${output} \
  -vcf ${vcf_file} \
  -gff ${ref_gff_file} \
  -t 1 \
  -tid 1 \
  -pop ${pop} \
  -merge_rgenome ${merge_rgenome} \
  -post_filtering ${post_filtering} \
  -n ${nchunks} \
  -samprate ${sample_ampl_rate} \
  -lamprate ${locus_ampl_rate} \
  -type_of_polymorphism_to_remove ${type_of_polymorphism_to_remove} \
  -ohetq ${ObsHet_quantile} \
  -snpdq ${SNP_density_quantile} \
  -alignment_filter ${alignment_filter}
  -mask_formula ${mask_formula}
