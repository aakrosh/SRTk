#!/usr/bin/env bash

indx=$1
lastz=$2
reference=$3
reads=$4
sd=$5
threads=$6
tmpdir=$7

lzscores="--match=1,1 --gap=3,1"
lzparams="--ambiguous=iupac --seed=match11 --xdrop=2 --hspthresh=20 --ydrop=5"
lzfmt="--format=softsam"
lzflt="--limitperquery=2011"

  ${lastz} \
  ${reference}[unmask,multi] \
  ${reads}[nameparse=full,subsample=${indx}/${threads}] \
  ${lzscores} ${lzparams} ${lzfmt} ${lzflt} \
| ${sd}/find_best_split -s 1 -i ${reference} \
> ${tmpdir}/alignments.${indx}.lz
