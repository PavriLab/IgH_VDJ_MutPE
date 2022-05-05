#!/bin/bash
set -e
TAG=$1

python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[CGcg][CTct][Cc]' no 'SYC' > reference_and_annotation/"$TAG"_SYC.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[Gg][AGag][CGcg]' no 'GRS' > reference_and_annotation/"$TAG"_GRS.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[CGcg][CTct][Cc]' yes 'SYC/GRS' > reference_and_annotation/"$TAG"_SYC_GRS.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[ATat][AGag][Cc]' no 'WRC' > reference_and_annotation/"$TAG"_WRC.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[Gg][CTct][ATat]' no 'GYW' > reference_and_annotation/"$TAG"_GYW.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[ATat][AGag][Cc]' yes 'WRC/GYW' > reference_and_annotation/"$TAG"_WRC_GYW.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[AGTagt][Gg][CTct][ATat]' no 'DGYW' > reference_and_annotation/"$TAG"_DGYW.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[AGTagt][Gg][CTct][ATat]' yes 'WRCH/DGYW' > reference_and_annotation/"$TAG"_WRCH_DGYW.bed
python sequtilities.py T reference_and_annotation/"$TAG".fa --regex2bed '[ATat][AGag][Cc][ACTact]' no 'WRCH' > reference_and_annotation/"$TAG"_WRCH.bed
