#!/bin/bash
#source /broad/software/scripts/useuse
#use R-4.1

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

fd=$(json_extract fd "$(cat ${json})")
fd=${fd%\"}
fd=${fd#\"}

vcf_file=$(json_extract vcf_file "$(cat ${json})")
vcf_file=${vcf_file%\"}
vcf_file=${vcf_file#\"}

output=$(json_extract output "$(cat ${json})")
output=${output%\"}
output=${output#\"}

HardFilteringReport=$(json_extract HardFilteringReport "$(cat ${json})")
remove_filtered_all=$(json_extract remove_filtered_all "$(cat ${json})")

include_from_gff=$(json_extract include_from_gff "$(cat ${json})")
include_from_gff=${include_from_gff%\"}
include_from_gff=${include_from_gff#\"}

exclude_from_gff=$(json_extract exclude_from_gff "$(cat ${json})")
exclude_from_gff=${exclude_from_gff%\"}
exclude_from_gff=${exclude_from_gff#\"}

ref_gff_file=$(json_extract ref_gff_file "$(cat ${json})")
ref_gff_file=${ref_gff_file%\"}
ref_gff_file=${ref_gff_file#\"}

bed=$(json_extract bed "$(cat ${json})")
bed=${bed%\"}
bed=${bed#\"}

exclude_bed=$(json_extract exclude_bed "$(cat ${json})")
exclude_bed=${exclude_bed%\"}
exclude_bed=${exclude_bed#\"}

keep_regexp=$(json_extract keep_regexp "$(cat ${json})")
keep_regexp=${keep_regexp%\"}
keep_regexp=${keep_regexp#\"}

positions=$(json_extract positions "$(cat ${json})")
positions=${positions%\"}
positions=${positions#\"}

non_ref_ac_any=$(json_extract non_ref_ac_any "$(cat ${json})") 

# Run pre-filtering

Rscript ${fd}/vcf2rGenome_filtering.R \
  -wd ${wd} \
  -fd ${fd} \
  -vcf ${vcf_file} \
  -o ${output} \
  -tid 1 \
  -HardFilteringReport ${HardFilteringReport} \
  -remove_filtered_all ${remove_filtered_all} \
  -include_from_gff ${include_from_gff} \
  -exclude_from_gff ${exclude_from_gff} \
  -gff ${ref_gff_file} \
  -bed ${bed} \
  -ebed ${exclude_bed} \
  -rkeep ${keep_regexp} \
  -positions ${positions} \
  -non_ref_ac_any ${non_ref_ac_any}
  