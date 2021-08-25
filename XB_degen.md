# XB degen
* Goal: evaluate evidence for degen ancestral Y and explore whether there is a turnover in the west
* Use genomic seqs mapped to XL to estimate polymorphism and dNdS in XB M and F
* Could explore sex-linkage in West Kenya samples using existing RADseq or with new GBS
* If new data, could add samples from SE Kenya near Kilimanjaro
* Drawback: best case scenario would include analysis of a family from west Kenya.
* Could use RADseq data to test whether west is fixed for W

# mapped RADseq data
evanslab:
```
/scratch/evanslab/2019_RADseq_KenyaXBXL_GhanaEastfamily/plate1/genotyped/mpileup_raw_justBorealis_allChrs.vcf.gz
```

## Mapping WGS data to XL

Working on this on graham: `/home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/`

* trimmed data:
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing a directory argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


module load StdEnv/2020
module load trimmomatic/0.39

#v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*_1.fq.gz; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
    echo java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-8}_1.fq.gz ${file::-8}_2.fq.gz ${file::-8}_trim.R1
.fq.gz ${file::-8}_trim.R1_single.fq.gz ${file::-8}_trim.R2.fq.gz ${file::-8}_trim.R2_single.fq.gz ILLUMINACLIP:TruSeq2_a
nd_3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-8}_1.fq.gz ${file::-8}_2.fq.gz ${file::-8}_trim.R1.
fq.gz ${file::-8}_trim.R1_single.fq.gz ${file::-8}_trim.R2.fq.gz ${file::-8}_trim.R2_single.fq.gz ILLUMINACLIP:TruSeq2_an
d_3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi
done 
```
Map to XL
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=96:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh ../../2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa.gz pathtofqfilez

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*_trim.R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	echo bwa mem ${1} ${file::-14}_trim.R1.fq.gz ${file::-14}_trim.R2.fq.gz -t 16 | samtools view -Shu - | samtools s
ort - -o ${file::-14}_sorted.bam
	bwa mem ${1} ${file::-14}_trim.R1.fq.gz ${file::-14}_trim.R1.fq.gz -t 16 | samtools view -Shu - | samtools sort -
 -o ${file::-8}_sorted.bam
	samtools index ${file::-14}_sorted.bam
  fi
done

```
Haplotypecaller:
```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=30gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh ref dir_of_bam chr
# sbatch 2021_HaplotypeCaller.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa /home/ben/proj
ects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/ 

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_sorted.bam_rg.bam
do
    gatk --java-options -Xmx24G HaplotypeCaller  -I ${file} -R ${1} -L ${3} -O ${file}_${3}.g.vcf -ERC GVCF
done
```
CombineGVCFs
```
#!/bin/sh
#SBATCH --job-name=CombineGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=12gb
#SBATCH --output=CombineGVCFs.%J.out
#SBATCH --error=CombineGVCFs.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch 2021_CombineGVCFs.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa 2021_new_borealis
_male_genome chr1L

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx10G CombineGVCFs -R ${1}"
for file in ${2}*{$3}*g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${3} -O ${2}MandF_${3}.g.vcf.gz"

${commandline}
```
