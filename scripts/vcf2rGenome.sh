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

nTasks=$(json_extract nTasks "$(cat ${json})")
Task_id=$(json_extract Task_id "$(cat ${json})")

wd=$(json_extract wd "$(cat ${json})")
wd=${wd%\"}
wd=${wd#\"}

fd=$(json_extract fd "$(cat ${json})")
fd=${fd%\"}
fd=${fd#\"}

output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

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

# Extract variables for vcf2rGenome
do_vcf2rGenome=$(json_extract do_vcf2rGenome "$(cat ${json})")
coverage_data=$(json_extract coverage_data "$(cat ${json})")
merge_rgenome='false'

nchunks=$(json_extract nchunks "$(cat ${json})")
ReadDepthThreshold=$(json_extract ReadDepthThreshold "$(cat ${json})")

metadata=$(json_extract metadata "$(cat ${json})")
metadata=${metadata%\"}
metadata=${metadata#\"}

join_by=$(json_extract join_by "$(cat ${json})")
join_by=${join_by%\"}
join_by=${join_by#\"}

pop=$(json_extract pop "$(cat ${json})")
pop=${pop%\"}
pop=${pop#\"}

# Change working directory
cd ${wd}

echo ${vcf_file}

# Run vcf2rGenome
Rscript ${fd}/vcf2rGenome_filtering.R \
  -wd ${wd} \
  -fd ${fd} \
  -vcf ${vcf_file} \
  -o ${output} \
  -t ${nTasks} \
  -tid ${SGE_TASK_ID} \
  -do_vcf2rGenome ${do_vcf2rGenome} \
  -coverage_data ${coverage_data} \
  -metadata ${metadata} \
  -join_by ${join_by} \
  -pop ${pop} \
  -merge_rgenome ${merge_rgenome} \
  -n ${nchunks} \
  -RDthres ${ReadDepthThreshold}

  