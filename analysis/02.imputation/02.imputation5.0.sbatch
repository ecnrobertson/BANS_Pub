#!/bin/bash
#
#SBATCH --job-name=impute
#SBATCH --output=impute.%j.out
#SBATCH --error=impute.%j.err
#SBATCH -t 6:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 
#SBATCH --mem=40G
#SBATCH --mail-type=END,FAIL
#SBATCH  --mail-user=ericacnr@colostate.edu

#echo commands to stdout

set -x
source ~/.bashrc

BEAGLE_JAR="/projects/ericacnr@colostate.edu/mambaforge/envs/bioinf/java/beagle.27Feb25.75f.jar"
vcf="/scratch/alpine/ericacnr@colostate.edu/BANS/01.fastq_processing/mega-non-model-wgs-snakeflow/results/bqsr-round-0/vcf_final/BANS.all.ds6x.pass-maf-0.05.bialSNP.filtered.ind0.3miss.0.8miss.recode.vcf"

java -Xss5m -Xmx90g -jar $BEAGLE_JAR \
  gt=$vcf \
  out=BANS.all.ds6x.pass-maf-0.05.bialSNP.filtered.ind0.3miss.0.8miss.recode.vcf_imputed5.0 \
  nthreads=48
