#!/usr/bin/env sh



${1}/fileutilities.py T ${2} --dir point.stats\$ | ${1}/fileutilities.py P --loop ${1}/postprocess_mpileup_summaries.R {abs} $3 $4 $5
