#!/usr/bin/env bash

for i in {1..22}; do
  bash create_namefile.sh hg002.mat_chr$i.fasta hg002.pat_chr$i.fasta ./raw_data $i pathfile.chr$i.txt;
done;

for i in {1..22}; do
  docker run --rm -it -v $(pwd):/data quay.io/comparative-genomics-toolkit/cactus:latest cactus-pangenome /data/work.chr$i /data/pathfile.chr$i.txt --outDir /data/out_chr$i --outName pan --reference hg002_mat;
done;

for i in {1..22}; do
  tantan -f4 hg002.mat_chr$i.fasta > ref.tantan.chr$i.bed;
done;

for i in {1..22}; do
  zcat out_chr$i/pan.gfa.gz | python3 processing.py --repeats ref.tantan.chr$i.bed --log log.chr$i.txt --chr $i --ref hg002_mat;
done
