#!/bin/bash
#
#SBATCH --job-name=impute
#SBATCH --output=impute.%j.out
#SBATCH --error=impute.%j.err
#SBATCH -t 6:00:00
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node 11
#SBATCH --mem=40G
#SBATCH --mail-type=END,FAIL
#SBATCH  --mail-user=ericacnr@colostate.edu

#echo commands to stdout

set -x
source ~/.bashrc

BEAGLE_JAR="/projects/ericacnr@colostate.edu/mambaforge/envs/bioinf/java/beagle.27Jan18.7e1.jar"
vcf="BANS.all.ds6x.pass-maf-0.05.bialSNP.filtered.ind03miss.0.8miss.rm_scaf_lth10.recode.vcf"

java -Xss5m -Xmx90g -jar $BEAGLE_JAR \
  gtgl=$vcf \
  out=BANS.all.ds6x.pass-maf-0.05.bialSNP.filtered.ind03miss.0.8miss.rm_scaf_lth10.recode.vcf \
  nthreads=48
