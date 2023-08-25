# Read Me: Genome Assembly and Search for Sex-Linked Contigs in Wikstroemia spp.

A 30 credit, 20 weeks long project that aims fill knowledge gaps of sex chromosomes by understanding how they evolve using *Wikstoemia oahuensis* as species of study and hopefully detect some sex-lined contigs. The scripts and commands run on UPPMAX and are adaptable for HPCs.

# Table of Contents

1. Navigation through SLURM
    1. Account Sign-In
    2. Basic Commands
    3. Sbatch Script Components
2. Navigating through the files and Quality Check
    1. FastQC Analysis
    2. Trimming 
    3. Post Trimming FastQC Analysis
3. Genome Assembly
    1. Genome Assembly using HiFiasm
    2. Converting the gfa files to Fasta
    3. Genome Quality Assessment
        1. BUSCO
        2. QUAST
4. Mapping 
    1. Trimming the Illumina Reads
    2. Mapping of Illumina Reads to the Generated Assembly
    3. Coverage Analysis
    4. Plotting the coverage analysis
5. Identification of sex-linked contigs: FindZX
    1. Run_1: All sample run
    2. Run_2: Sample selection run
    3. Run_3: Selected sample run

# 0) Navigation through SLURM

---

What all Simple Linux Utility for Resource Management (SLURM) does?

- job submission using sbatch command
- resource allocation
- job scheduling
- job monitoring

## 0.1) Accounts Signed-In

- [SNIC 2022/5-71](https://supr.snic.se/project/21759/)
- [SNIC 2021/2-12](https://supr.snic.se/project/18967/)

in Rakham Cluster

## 0.2) SLURM Basic Commands

1. sbatch : run a bash script
    
    ```bash
    sbatch my_job_script_file.sh
    ```
    
2. interactive : run an interactive session 
    
    ```bash
    interactive -A <project_name> -t <DD-hh:mm:ss>
    ```
    
3. jobinfo : to get status and other information about job
    
    ```bash
    jobinfo -u <user_name>
    #OR
    jobinfo -A account
    
    ```
    
4. scancel : to cancel a job 
    
    ```bash
    scancel -n <job_name>
    ```
    

## 0.3) Sbatch Script Components

| #!/bin/bash -l | shebang to make file executable  |
| --- | --- |
| #SBATCH -A snic2022-5-71 | project name  |
| #SBATCH -p node | nodes |
| #SBATCH -n 20 | number of nodes |
| #SBATCH -t 2-00:00:00 | tine for the run  |
| #SBATCH -J hifiasm_wikstroemia_hmc20s55.job | job name; it will help to search the project status  |
| #SBATCH -o hifiasm_wikstroemia_hmc20s55.out | output file name  |
| #SBATCH -e hifiasm_wikstroemia_hmc20s55.err | error file name |
| #SBATCH --mail-user=ayushipathakofficial@gmail.com | email for the job updates like queued, running, ending and cancelled |
| #SBATCH --mail-type=ALL | what kind of mailing listyou prefer |

# 1) Navigating through the files and Quality Check

Make new folders in the designated location 

```bash
cd /proj/snic2020-2-25/nobackup/ayushi
mkdir 1_HiFiasm 2_gfa_files 3_fasta_files 4_quast_output 5_busco_output 6_findZX_out 7_trimmed_reads 8_map_female_and_male_reseq_to_assembly 9_get_con_reseq
```

Locations of the of the PacBio reads and Illumina reads

```bash
#PacBio Reads
/proj/snic2020-2-25/nobackup/ayushi/hifi_reads 
#Illumina Reads
/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/*Woah*
```

# 2) Genome Assembly

---

Software: HiFiasm

HiFiasm version 0.16.0 

GitHub: [GitHub - chhylp123/hifiasm: Hifiasm: a haplotype-resolved assembler for accurate Hifi reads](https://github.com/chhylp123/hifiasm)

Usage: hifiasm <input file> -o <output file name> 

Installation 

```bash
cd bin 
git clone https://github.com/chhylp123/hifiasm.git
```

Exploring the flags and checking if installation was complete

```bash
/home/ayuship/bin/hifiasm-0.16.1/hifiasm -h
```

## 2.1) Genome Assembly using HiFiasm

### 2.1.1) Firstrun

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_wikstroemia_firstrun.job
#SBATCH -o hifiasm_wikstroemia_firstrun.out
#SBATCH -e hifiasm_wikstroemia_firstrun.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#creating directory
outdir="/proj/snic2020-2-25/nobackup/ayushi"
readdir="/proj/snic2020-2-25/nobackup/ayushi/hifi_reads"

cd $readdir

#copying seq to temporary folder
cp *.fastq.gz $SNIC_TMP

cd $SNIC_TMP

#making needed directories
mkdir hifiasm_wikstroemia_firstrun

#running the code in specified folder
cd hifiasm_wikstroemia_firstrun

/home/ayuship/bin/hifiasm-0.16.1/hifiasm   -o hifiasm_wikstroemia_firstrun -t20 ../*.fastq.gz

cd ..

cp -r hifiasm_wikstroemia_firstrun $outdir
```

### 2.1.2) k60

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_wikstroemia_k60.job
#SBATCH -o hifiasm_wikstroemia_k60.out
#SBATCH -e hifiasm_wikstroemia_k60.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

outdir="/proj/snic2020-2-25/nobackup/ayushi"
readdir="/proj/snic2020-2-25/nobackup/ayushi/hifi_reads"

cd $readdir

cp *.fastq.gz $SNIC_TMP

cd $SNIC_TMP

mkdir hifiasm_wikstroemia_k60

cd hifiasm_wikstroemia_k60

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  -k 60 -o hifiasm_wikstroemia_k60 -t20 ../*.fastq.gz

cd ..

cp -r hifiasm_wikstroemia_k60 $outdir
```

As the assembly size was larger than expected so --hom-cov and -s flags were introduced with combination of various settings

Flag specification :

| --hom-cov | adjustment of homologous coverage |
| --- | --- |
| -s | used to balance out the length of assembly due to differences in the paternal and maternal haplotypes |

2.1.3) Hmc20s45, Hmc22s45 and Hmc24s45

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_wikstroemia_hmc20s45.job
#SBATCH -o hifiasm_wikstroemia_hmc20s45.out
#SBATCH -e hifiasm_wikstroemia_hmc20s45.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

outdir="/proj/snic2020-2-25/nobackup/ayushi"
readdir="/proj/snic2020-2-25/nobackup/ayushi/hifi_reads"

cd $readdir

cp *.fastq.gz $SNIC_TMP

cd $SNIC_TMP

mkdir hifiasm_wikstroemia_hmc20s45 hifiasm_wikstroemia_hmc22s45 hifiasm_wikstroemia_hmc24s45

cd hifiasm_wikstroemia_hmc20s45

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 20 -s 0.45 -o hifiasm_wikstroemia_hmc20s45 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc22s45

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 22 -s 0.45 -o hifiasm_wikstroemia_hmc22s45 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc24s45

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 24 -s 0.45 -o hifiasm_wikstroemia_hmc24s45 -t20 ../*.fastq.gz

cd ..

cp -r hifiasm_wikstroemia_* $outdir
```

### 2.1.4) Hmc20s50, Hmc22s50 and Hmc24s50

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_wikstroemia_hmc20s50.job
#SBATCH -o hifiasm_wikstroemia_hmc20s50.out
#SBATCH -e hifiasm_wikstroemia_hmc20s50.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

outdir="/proj/snic2020-2-25/nobackup/ayushi"
readdir="/proj/snic2020-2-25/nobackup/ayushi/hifi_reads"

cd $readdir

cp *.fastq.gz $SNIC_TMP

cd $SNIC_TMP

mkdir hifiasm_wikstroemia_hmc20s50 hifiasm_wikstroemia_hmc22s50 hifiasm_wikstroemia_hmc24s50

cd hifiasm_wikstroemia_hmc20s50

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 20 -s 0.50 -o hifiasm_wikstroemia_hmc20s50 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc22s50

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 22 -s 0.50 -o hifiasm_wikstroemia_hmc22s50 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc24s50

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 24 -s 0.50 -o hifiasm_wikstroemia_hmc24s50 -t20 ../*.fastq.gz

cd ..

cp -r hifiasm_wikstroemia_* $outdir
```

### 2.1.5) Hmc20s55, Hmc22s55 and Hmc24s55

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 2-00:00:00
#SBATCH -J hifiasm_wikstroemia_hmc20s55.job
#SBATCH -o hifiasm_wikstroemia_hmc20s55.out
#SBATCH -e hifiasm_wikstroemia_hmc20s55.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

outdir="/proj/snic2020-2-25/nobackup/ayushi"
readdir="/proj/snic2020-2-25/nobackup/ayushi/hifi_reads"

cd $readdir

cp *.fastq.gz $SNIC_TMP

cd $SNIC_TMP

mkdir hifiasm_wikstroemia_hmc20s55 hifiasm_wikstroemia_hmc22s55 hifiasm_wikstroemia_hmc24s55

cd hifiasm_wikstroemia_hmc20s55

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 20 -s 0.55 -o hifiasm_wikstroemia_hmc20s55 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc22s55

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 22 -s 0.55 -o hifiasm_wikstroemia_hmc22s55 -t20 ../*.fastq.gz

cd ..

cd hifiasm_wikstroemia_hmc24s55

/home/ayuship/bin/hifiasm-0.16.1/hifiasm  --hom-cov 24 -s 0.55 -o hifiasm_wikstroemia_hmc24s55 -t20 ../*.fastq.gz

cd ..

cp -r hifiasm_wikstroemia_* $outdir
```

## 2.2) Converting the gfa files to Fasta

The output from hifiasm is in gfa format that needs to be converted to fasta format for further quality assessment.

```bash
for file in /proj/snic2020-2-25/nobackup/ayushi/hifiasm_wikstroemia_*/*.bp.p_ctg.gfa;do name=$(echo $file |cut -d '/' -f 7 | cut -d '.' -f 1-3); awk '/^S/{print ">"$2;print $3}' $file > /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/$name'.fa'; done
```

code breakdown:

| for file in /proj/snic2020-2-25/nobackup/ayushi/hifiasm_wikstroemia_*/*.bp.p_ctg.gfa; | creating the loop for the needed file format that is bp.p_ctg.gfa |
| --- | --- |
| do name=$(echo $file |cut -d '/' -f 7 | cut -d '.' -f 1-3); | creating the variable ‘name’ without file extensions |
|  awk '/^S/{print ">"$2;print $3}' $file | awk command for extracting useful information form file  |
| > /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/$name'.fa'; done | saving file with proper file foramat at specified location and closing the loop |

## 2.3) Genome Quality Assessment

Quast reports summary statistics such as N50 and overall size of the assembly to give an idea how contiguous the assembly is. BUSCO checks for the presence of lineage-specific single copy genes. Better assemblies are expected to have higher numbers of complete genes and lower number of duplicated genes.

### 2.3.1) BUSCO

Software: BUSCO

Version: 5.3.1

GitHub: [GitHub - metashot/busco: Assessing the quality of genomes using busco](https://github.com/metashot/busco)

Usage: run_BUSCO.py -i <SEQUENCE_FILE> -l <LINEAGE> -o <OUTPUT_NAME> -m <MODE> <OTHER OPTIONS>

Flags: 

| --in | input file name  |
| --- | --- |
| --out | output file name |
| --lineage | lineage closest to the species in use |
| --mode | for genomic data use mode as genome |

Troubleshooting BUSCO error

```bash
module load bioinfo-tools
module load BUSCO/5.3.1
#Error: Cannot write to Augustus directory, please make sure you have write permissions to <directory> then you need to create a local copy of the config directory from the augustus module with 'source $AUGUSTUS_CONFIG_COPY'.  See 'module help augustus/3.4.0'
module help augustus/3.4.0
#Lmod Warning:  Failed to find the following module(s): "augustus/3.4.0" in your MODULEPATH Try: $ module spider augustus/3.4.0 to see if the module(s) are available across all compilers and MPI implementations.

'''
**If you see this error:

          Error: Cannot write to Augustus directory, please make sure you have write permissions to <directory>

      then you likely need to create a local copy of the config directory from
      the augustus module.  With this augustus module loaded, you can do this with:

          source $AUGUSTUS_CONFIG_COPY

      This will create a directory named augustus_config/ in the current directory
      and the AUGUSTUS_CONFIG_PATH environment variable will be updated to reflect
      this new location.  If you continue to see the error from BUSCO, the script
      was not run correctly.  This local directory will contain additional training
      sets created by BUSCO during its run.

      If you already have your own directory for augustus config files, then you may
      wish to adjust AUGUSTUS_CONFIG_PATH to name this directory.
'''**

```

Sourcing the Augustus config file 

```bash
**source $AUGUSTUS_CONFIG_COPY
# creates a augustus_config in current folder**
```

Checking if installation was successful

```bash
run_BUSCO.py --help
```

Running the BUSCO script

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 6-00:00:00
#SBATCH -J BUSCO_run.job
#SBATCH -o BUSCO_run.out
#SBATCH -e BUSCO_run.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#loading the modules 
moldule load bioinfo-tools
module load BUSCO/5.3.1

#Specify input file directory
input_dir='/proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/'

#Specify output file directory
out_dir-'/proj/snic2020-2-25/nobackup/ayushi/5_busco_output'

#copying input files to temporary folder
cd $input_dir
cp *.fa $SNIC_TMP
cd $SNIC_TMP

for file in $(ls | grep fa$); do file_name=$( echo $file | cut -d '_' -f 3| cut -d '.' -f 1) ;run_BUSCO.py --in $file --lineage $BUSCO_LINEAGE_SETS/eudicots_odb10 --out $file_name'_BUSCO_out' --mode genome  ;cp -r $file_name'_BUSCO_out' $out_dir ; done
```

### 2.3.2) QUAST

Software: QUAST 

Version: 0.39 

GitHub: [GitHub - ablab/quast: Genome assembly evaluation tool](https://github.com/ablab/quast)

Usage: 

```
./quast.py test_data/contigs_1.fasta \
           test_data/contigs_2.fasta \
        -r test_data/reference.fasta.gz \
        -g test_data/genes.txt \
        -1 test_data/reads1.fastq.gz -2 test_data/reads2.fastq.gz \
        -o quast_test_output
```

```bash
quast.py -h
```

Quast analysis is not too heavy so it can run a loop 

```bash
for file in /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/*.fa;do file_name=$(echo $file | cut -d '/' -f 7 | cut -d '.' -f 1); echo $file_name;output=('quast_output_'$file_name); python /home/ayuship/bin/quast/quast.py $file -o /proj/snic2020-2-25/nobackup/ayushi/4_quast_output/$output  ; done
```

Quast script

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 04:00:00
#SBATCH -J QUAST_run.job
#SBATCH -o QUAST_run.out
#SBATCH -e QUAST_run.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#loading the modules/tools
moldule load bioinfo-tools
module load quast/5.0.2

#providing the output directory
out_dir='/proj/snic2020-2-25/nobackup/ayushi/4_quast_output'

#specifiy input directory 
in_dir='/proj/snic2020-2-25/nobackup/ayushi/3_fasta_files'

#copy fast files
cd $in_dir
cp *fa $SNIC_TMP
cd $SNIC_TMP

prefixes=$(ls *.fa  | cut -d "_" -f3 | cut -d '.' -f1)
for prefix in $prefixes
do
python /sw/bioinfo/quast/5.0.2/rackham/bin/quast.py hifiasm_wikstroemia_${prefix}.bp.p_ctg.fa -o $out_dir/quast_out_${prefix}
done
```

# 3) Mapping

---

The mapping step consists of :

- trimming
- indexing of assembled genome
- mapping of reads
- conversion of generated .sam file to .bam file
- sorting the bam file v) removing duplicates
- indexing of sorted deduplicated file.

## 3.1) Trimming the Illumina Reads

Software: Trimmomatic

Version: 0.39 

GitHub: [GitHub - timflutre/trimmomatic: Read trimming tool for Illumina NGS data.](https://github.com/timflutre/trimmomatic)

Usage: java -jar $TRIMMOMATIC_ROOT/trimmomatic.jar

Flags:

| PE | Paired ends |
| --- | --- |
| -basein  | input file format |
| -baseout | output file format |
| ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE-2.fa:2:30:10  | adapter : /TruSeq3-PE-2.fa |
|  LEADING:15  | remove leading (in the starting) low quality bases (below quality 15) |
| TRAILING:20  | remove trailing (at the end) low quality bases (below quality 15) |
| SLIDINGWINDOW:4:20 | reads are scanned with a 4 base wide sliding window , cutting when quality (average) drops below 15  |
|  MINLEN:90 | remove the reads below 90 bp length |

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 8-00:00:00
#SBATCH -J wikstroemia_trimmomatic.job
#SBATCH -o wikstroemia_trimmomatic.out
#SBATCH -e wikstroemia_trimmomatic.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

#loading the modules
module load bioinfo-tools trimmomatic

#specify the directories
raw_read_dir="/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX"
out_dir="/proj/snic2020-2-25/nobackup/ayushi/7_trimmed_reads"

cd $raw_read_dir

#select specific files form folders
files=$(find | grep Woah | grep gz)

for file in $files
do
cp $file $SNIC_TMP
done

cd $SNIC_TMP

prefixes=$(ls *.fastq.gz  | cut -d "_" -f1-3)
for prefix in $prefixes
do
trimmomatic PE -threads 20 -basein ${prefix}_R1_001.fastq.gz -baseout ${prefix}_trimmed.fastq.gz ILLUMINACLIP:$TRIMMOMATIC_ROOT/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:15 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:90
cp ${prefix}_trimmed_[12]P.fastq.gz $out_dir
done
```

## 3.2) Mapping of Illumina Reads to the Generated Assembly

Software: BWA (Burrows-Wheeler Aligner)

Version: 0.7.17

GitHub: [GitHub - lh3/bwa: Burrow-Wheeler Aligner for short-read alignment (see minimap2 for long-read alignment)](https://github.com/lh3/bwa)

Usage: bwa mem <reference_genome> <read_1.fq> <read_2.fq> > <output_file>

Flags:

| mem | used  for Illumina reads ~70 bp length |
| --- | --- |
| -t | number of threads |
| -M  | to map short split hits as secondary  |
| -R | to create group header in this case : "@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:Illumina" |

Software: samtools 

Version: 1.9

GitHub: [Releases · samtools/samtools (github.com)](https://github.com/samtools/samtools/releases/)

Usage: samtools <commands> [additional_options]

Flags:

| view | sam to bam conversion (in this case) |
| --- | --- |
| sort  | sort alignment from leftmost coordinate |
| @ | number of sorting and compression threads  |
| index  | index the alignment |

Software: picard

Version: 2.27.5

GitHub: [GitHub - broadinstitute/picard: A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF.](https://github.com/broadinstitute/picard)

Usage: java -jar $PICARD_ROOT/picard.jar MarkDuplicates <input_file > <out_file> <metrics_file>

Flags:

| MarkDuplicates | command to identify dulicate reads from file |
| --- | --- |
| REMOVE_DUPLICATES=true | removing the duplicates from original file  |
| I | input file  |
| O | output file |
| M | metrics file  |

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p node
#SBATCH -n 20
#SBATCH -t 5-00:00:00
#SBATCH -J wik_map_female_and_male_reseq_to_assembly.job
#SBATCH -o wik_map_female_and_male_reseq_to_assembly.out
#SBATCH -e wik_map_female_and_male_reseq_to_assembly.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools samtools/1.9 bwa/0.7.17 picard

#Provide trimmed read dir
trimmed_read_dir="/proj/snic2020-2-25/nobackup/ayushi/7_trimmed_reads/"

#Provide the name of the genome fasta file
genome_dir="/proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/"

#Provide the name of the genome fasta file
genome="hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa"

#Specify an output directory - need to create it beforehand
output_dir="/proj/snic2020-2-25/nobackup/ayushi/8_map_female_and_male_reseq_to_assembly"

#Copy trimmed reads to the working directory
cd $trimmed_read_dir
cp *.gz $SNIC_TMP

#Copy the genome assembly to the working directory
cd $genome_dir
cp $genome $SNIC_TMP

#Go to the working directory
cd $SNIC_TMP

#Index the genome before mapping
bwa index $genome

#Loop through the samples and map them separately

samples=$(ls *.fastq.gz | cut -f1-2 -d "_" | sort | uniq)

for sample in $samples
do
bwa mem -t19 -M -R "@RG\tID:$sample\tLB:$sample\tSM:$sample\tPL:Illumina" $genome ${sample}_L001_trimmed_[12]P.fastq.gz  | samtools view -bS - > $sample.bamsamtools sort -@ 20 $sample.bam > $sample.sorted.bam
java -jar -Xmx110G $PICARD_TOOLS_DIR/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$sample.sorted.bam O=$sample.sorted.dedup.bam M=$sample.duplicatedata.txt
samtools index -@ 20 $sample.sorted.dedup.bam
cp $sample.sorted.dedup.bam* $output_dir
cp $sample.duplicatedata.txt $output_dir
done
```

## 3.3) Coverage Analysis

Software: BEDTools

Version: 2.29.2

GitHub: [arq5x/bedtools2: bedtools - the swiss army knife for genome arithmetic (github.com)](https://github.com/arq5x/bedtools2)

Usage: 

Flags:

| makewindows | used to create windows dor analysis |
| --- | --- |
| -w | window size |
| -g | genome file  |
| multicov | creating multi coverage plots/ tables that is counting coverage from multiple files  |
| -q | quality mapping  |
| -p | for paired end reads  |
| -bed | follow bed format  |

```bash
#!/bin/bash -l
#SBATCH -A snic2022-5-71
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 1-16:00:00
#SBATCH -J wik_generic_get_cov_female_and_male_reference.job
#SBATCH -o wik_generic_get_cov_female_and_male_reference.out
#SBATCH -e wik_generic_get_cov_female_and_male_reference.err
#SBATCH --mail-user=ayushipathakofficial@gmail.com
#SBATCH --mail-type=ALL

module load bioinfo-tools BEDTools/2.29.2 samtools

#Specify output directory
output_dir="/proj/snic2020-2-25/nobackup/ayushi/9_get_con_reseq"

#Specify genome directory
genome_dir="/proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/"

#Specify genome assembly fasta file
genome="hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa"

#Specify directory where the mapped reads (bam files) are
map_dir="/proj/snic2020-2-25/nobackup/ayushi/8_map_female_and_male_reseq_to_assembly/"

#Copy the genome assembly to the working directory
cd $genome_dir
cp $genome $SNIC_TMP

#Copy mapped reads to the working directory
cd $map_dir
cp *.bam* $SNIC_TMP

#Go to the working directory
cd $SNIC_TMP

#Get 1 kb window intervals from the genome
window_size=1000
samtools faidx $genome
cat $genome.fai | cut -f1-2 | sort -k1,1 > $genome.contig_sizes.txt
bedtools makewindows -w $window_size -g $genome.contig_sizes.txt > $genome.1kb.bed

#Get read cov for each of the samples in 1 kb windows. -q1: mapping quality of at least one; -p paired end reads

bedtools multicov -q 1 -p -bams SJ-2333-Woah-000_S4.sorted.dedup.bam SJ-2333-Woah-003_S2.sorted.dedup.bam SJ-2333-Woah-006_S10.sorted.dedup.bam SJ-2333-Woah-996_S9.sorted.dedup.bam SJ-2333-Woah-002_S12.sorted.dedup.bam SJ-2333-Woah-005_S6.sorted.dedup.bam SJ-2333-Woah-007_S11.sorted.dedup.bam SJ-2333-Woah-998_S8.sorted.dedup.bam -bed $genome.1kb.bed > $genome.female_male_cov.$window_size.windows.bed

cp $genome.female_male_cov.$window_size.windows.bed $output_dir
```

## 3.4) Plotting the coverage analysis

```bash
#Go to the directory that stores the coverage bed file

cd /proj/snic2020-2-25/nobackup/ayushi/9_get_con_reseq
```

```bash
#FCreate a directory that we will contain the plots

mkdir plots_male_female_assembly
```

R script to visualize the plots

```r
#Read window coverage file and specify plot directory

cov_data=read.delim("hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa.female_male_cov.1000.windows.bed",header=F)

plot_dir="plots_male_female_assembly"

#Get median genome-wide coverage for each of the samples

n_samples=8

med_cov=rep(0,n_samples)

for(i in c(1:n_samples)){
 med_cov[i]=median(cov_data[,i+3])
}

#Calculate normalized coverage in each window
#by, for each sample, divide values with the median coverage for that sample
norm_cov_data=cov_data

for(i in c(1:n_samples)){
 norm_cov_data[,i+3]=cov_data[,i+3]/med_cov[i]
}

#Get a mean coverage for males and females, respectively
#Here it is very important to correctly specify the male and female intevals

#SJ-2333-Woah-000_S4 (1:male) SJ-2333-Woah-003_S2 (2:female) SJ-2333-Woah-006_S10 (3:male) SJ-2333-Woah-996_S9 (4:female) SJ-2333-Woah-002_S12 (5:female) SJ-2333-Woah-005_S6 (6:female) SJ-2333-Woah-007_S11 (7:male) SJ-2333-Woah-998_S8 (8:male)
#male samples: 1,3,7,8
#female samples: 2,4,5,6
#convert these into column indices by adding 3 (first three column have positional information)

male_mean=rowMeans(norm_cov_data[,c(4,6,10,11)])
female_mean=rowMeans(norm_cov_data[,c(5,7,8,9)])

norm_cov_data=cbind(norm_cov_data,male_mean,female_mean)

#Plot for each contig

contigs=unique(norm_cov_data[,1])

for(j in c(1:length(contigs))){

	contig=contigs[j]
	contig_data=norm_cov_data[which(norm_cov_data[,1]==contig),]

	filename=gsub(" ","",paste(plot_dir,"/",contig,".png"))

	png(file=filename,width=960,height=480)

	#par(mfrow=c(2,1))

	plot(rowMeans(contig_data[,2:3]),contig_data$male_mean,ylim=c(0,3),col="blue",type="l",xlab="position",ylab="norm cov",main=contig)
	points(rowMeans(contig_data[,2:3]),contig_data$female_mean,col="green",type="l")

	#plot(rowMeans(scaff_data[,2:3]),log2(scaff_data$male_mean/scaff_data$female_mean),ylim=c(-2,2),type="l",xlab="position",ylab="norm cov",main=scaffold)

	dev.off()

	}
```

# 4) Identification of sex-linked contigs

---

Software/ Pipeline: FindZX

GitHub: [GitHub - hsigeman/findZX](https://github.com/hsigeman/findZX)

Usage: snakemake -s workflow/findZX 

Installation 

```bash
git clone https://github.com/hsigeman/findZX.git
cd findZX # Go to directory
```

Create conda environment 

```bash
conda create -n findZX -c conda-forge -c bioconda python=3.9.4 snakemake-wrapper-utils=0.2.0 snakemake=6.4.0 mamba=0.15.3
```

| -n | name of conda environment |
| --- | --- |
| -c  | needed channels |
|  python=3.9.4 | python version  |
| snakemake | installs snakemake version 6.4.0 |
| mamba | installa mamba version 0.15.3 |

Activating the environment

```bash
conda activate findZX
```

After activating the environment check the installation

```bash
python -V # Should give Python 3.9.4
```

Running basic command

```bash
snakemake -s workflow/findZX --configfile .test/config.yml --cores 1 -R all --use-conda -k
snakemake -s workflow/findZX-synteny --configfile .test/config.yml --cores 1 -R all -k --use-conda
```

basic command to run on SLURM

```bash
snakemake -s workflow/findZX -j 15 -R all --configfile test_config.yml --cluster-config .test/config.yml --cluster " sbatch -A {cluster.account} -t {cluster.time} -n {cluster.n} " -k --use-conda
```

|  --configfile | config file with path and settings  |
| --- | --- |
| --cluster | the cluster specification |
| -A | cluster account |
| -t | cluster time  |
| -n | cluster nodes |
| -R | generates all the output files  |
| -k | all the others jobs continue when one fails |
| --use-conda | using conda in different environments |

Inside FindZX folder you find various subfolders that are used .

```bash
#inside FindZX folder the following subfolders are present
README.md   cluster  environment.yml  reports  workflow  Rplots.pdf  config   figures  results
```

All the config files are stored in **config folder**. This folder also contain the units files in tsv format that contains the sample, group, fq1 and fq2. These are the included samples in run.

NOTE : All the config files format is taken from the github repository.

The basic command in FindZX remains almost same. The main edits are made in the config files. These files open us to many editing and visualizing options.

## 4.1 **Run_1**: All sample run

```bash
snakemake -s workflow/findZX -j 15 -R all --configfile config/Run_1.yml --cluster-config cluster/cluster_Wikstroemia.yml --cluster " sbatch -A snic2022-5-71 -t 8:00:00 -n 3 " -k --use-conda
```

config file for Run_1

```bash
## ================================= ##
## findZX config file (test dataset) ##
## ================================= ##

# Variables marked with "[findZX]" or "[findZX-synteny]" are only used when deploying 
# the pipeline with either snakefile. Other variables are used for all analyses. 

threads_max: 18
mem_max: 16000
  # Specify maximum number of cores [threads_max] and memory [mem_max] allocation

# ============================ #
# Analysis name and input data #

run_name: Wikstroemia/scaf_assembly
  # Select an analysis name. Output files will be stored under "results/[run_name]"

units: config/Wikstroemia/scaf_assembly/Wikstroemia_units.tsv
  # Path to sample information file

ref_genome: /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa 
  # Path to study-species reference genome (not .gz format)

synteny_ref: 
  # [findZX-synteny] Path to synteny-species reference genome (not .gz format)

synteny_name: 
  # [findZX-synteny] Synteny-species name (can be any string, will be used for file and directory names)

# ================= #
# Plotting settings #

window_sizes: [50000, 100000, 1000000]
  # Choose genome window sizes for plotting (as many as you want)
  # Optimal sizes depend on reference genome fragmentation and size of the sex-linked region
  # Recommended sizes to start are: [50000, 100000, 1000000] (i.e. 50 kb, 100 kb, 1Mb)

chr_file: None
  # [findZX] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

chr_highlight:  
  # [findZX] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

synteny_chr_file: None
  # [findZX-synteny] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

synteny_chr_highlight:
  # [findZX-synteny] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

# ================================== #
# Trimming and subsampling of reads  #

## These three variables control trimming and subsampling of reads
## Set all to "false" to disable trimming and subsampling
## Only one variable is allowed to be "true"

trim_reads: true 
  # Set to true for trimming of reads

trim_and_subsample: false
  # Set to true for trimming and subsampling of reads

subsample_only: false
  # Set to true for subsampling of reads (but not trimming)

subsample_basepairs: 1888226
  # Specify the total number of basepairs to extract from both fastq files
  # Will be used if [trim_and_subsample] or [subsample_basepairs] is set to "true"

  # Use this script to calculate expected coverage:
  # ./code/subsampling_cov_calv.sh <REF.fasta> <WANTED_COV> 

# ========================== #
# findZX-specific parameters # 
# ===== (edit if needed) ===== #
 
mismatch_settings: [0.0, 0.2]
  # Genome coverage results will be generated from the original BAM files ("unfiltered"),
  # and two other (modifiable) mismatches settings.
  # "0.0" = 0 mismatches allowed
  # "0.2" = <=2 mismatches allowed

minSizeScaffold: "7000"
  # The mimimum size of scaffolds in the reference genome to be included in the results

# ============================ #
# External software parameters # 
# ===== (edit if needed) ===== #

params:
  trimmomatic:
  # Control Trimmomatic settings here
    pe:
      trimmer:
        - "LEADING:20"
        - "TRAILING:20"
        - "SLIDINGWINDOW:4:20"
        - "MINLEN:120"
        - "ILLUMINACLIP:workflow/meta/adapters/TruSeq3-PE.fa:2:30:10"
```

Units file with path of samples

File name : Wikstroemia_units.tsv

```
sample	group	fq1	fq2
Woah-000	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-000/SJ-2333-Woah-000_S4_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-000/SJ-2333-Woah-000_S4_L001_R2_001.fastq.gz
Woah-002	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-002/SJ-2333-Woah-002_S12_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-002/SJ-2333-Woah-002_S12_L001_R2_001.fastq.gz
Woah-003	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-003/SJ-2333-Woah-003_S2_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-003/SJ-2333-Woah-003_S2_L001_R2_001.fastq.gz
Woah-005	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-005/SJ-2333-Woah-005_S6_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-005/SJ-2333-Woah-005_S6_L001_R2_001.fastq.gz
Woah-006	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-006/SJ-2333-Woah-006_S10_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-006/SJ-2333-Woah-006_S10_L001_R2_001.fastq.gz
Woah-007	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-007/SJ-2333-Woah-007_S11_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-007/SJ-2333-Woah-007_S11_L001_R2_001.fastq.gz
Woah-996	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-996/SJ-2333-Woah-996_S9_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-996/SJ-2333-Woah-996_S9_L001_R2_001.fastq.gz
Woah-998	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-998/SJ-2333-Woah-998_S8_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-998/SJ-2333-Woah-998_S8_L001_R2_001.fastq.gz
```

## 4.1 **Run_2**: Sample selection run

```
snakemake -s workflow/findZX -j 15 -R all --configfile config/Run_2.yml --cluster-config cluster/cluster_Wikstroemia.yml --cluster " sbatch -A snic2022-5-71 -t 8:00:00 -n 3 " -k --use-conda
```

config file for Run_2

```bash
## WITHOUT WOAH-000 ###########################################################################################################################################
## ================================= ##
## findZX config file (test dataset) ##
## ================================= ##

# Variables marked with "[findZX]" or "[findZX-synteny]" are only used when deploying 
# the pipeline with either snakefile. Other variables are used for all analyses. 

threads_max: 18
mem_max: 16000
  # Specify maximum number of cores [threads_max] and memory [mem_max] allocation

# ============================ #
# Analysis name and input data #

run_name: Wikstroemia/scaf_assembly
  # Select an analysis name. Output files will be stored under "results/[run_name]"

units: config/Wikstroemia/scaf_assembly/Wikstroemia_units_no_Woah000.tsv
  # Path to sample information file

ref_genome: /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa 
  # Path to study-species reference genome (not .gz format)

synteny_ref: 
  # [findZX-synteny] Path to synteny-species reference genome (not .gz format)

synteny_name: 
  # [findZX-synteny] Synteny-species name (can be any string, will be used for file and directory names)

# ================= #
# Plotting settings #

window_sizes: [50000, 100000, 1000000]
  # Choose genome window sizes for plotting (as many as you want)
  # Optimal sizes depend on reference genome fragmentation and size of the sex-linked region
  # Recommended sizes to start are: [50000, 100000, 1000000] (i.e. 50 kb, 100 kb, 1Mb)

chr_file: None
  # [findZX] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

chr_highlight:  
  # [findZX] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

synteny_chr_file: None
  # [findZX-synteny] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

synteny_chr_highlight:
  # [findZX-synteny] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

# ================================== #
# Trimming and subsampling of reads  #

## These three variables control trimming and subsampling of reads
## Set all to "false" to disable trimming and subsampling
## Only one variable is allowed to be "true"

trim_reads: true 
  # Set to true for trimming of reads

trim_and_subsample: false
  # Set to true for trimming and subsampling of reads

subsample_only: false
  # Set to true for subsampling of reads (but not trimming)

subsample_basepairs: 1888226
  # Specify the total number of basepairs to extract from both fastq files
  # Will be used if [trim_and_subsample] or [subsample_basepairs] is set to "true"

  # Use this script to calculate expected coverage:
  # ./code/subsampling_cov_calv.sh <REF.fasta> <WANTED_COV> 

# ========================== #
# findZX-specific parameters # 
# ===== (edit if needed) ===== #
 
mismatch_settings: [0.0, 0.2]
  # Genome coverage results will be generated from the original BAM files ("unfiltered"),
  # and two other (modifiable) mismatches settings.
  # "0.0" = 0 mismatches allowed
  # "0.2" = <=2 mismatches allowed

minSizeScaffold: "7000"
  # The mimimum size of scaffolds in the reference genome to be included in the results

# ============================ #
# External software parameters # 
# ===== (edit if needed) ===== #

params:
  trimmomatic:
  # Control Trimmomatic settings here
    pe:
      trimmer:
        - "LEADING:20"
        - "TRAILING:20"
        - "SLIDINGWINDOW:4:20"
        - "MINLEN:120"
        - "ILLUMINACLIP:workflow/meta/adapters/TruSeq3-PE.fa:2:30:10"
```

Units file with path of samples

File name : Wikstroemia_units_no_Woah000.tsv

```
sample	group	fq1	fq2
Woah-002	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-002/SJ-2333-Woah-002_S12_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-002/SJ-2333-Woah-002_S12_L001_R2_001.fastq.gz
Woah-003	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-003/SJ-2333-Woah-003_S2_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-003/SJ-2333-Woah-003_S2_L001_R2_001.fastq.gz
Woah-005	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-005/SJ-2333-Woah-005_S6_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-005/SJ-2333-Woah-005_S6_L001_R2_001.fastq.gz
Woah-006	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-006/SJ-2333-Woah-006_S10_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-006/SJ-2333-Woah-006_S10_L001_R2_001.fastq.gz
Woah-007	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-007/SJ-2333-Woah-007_S11_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-007/SJ-2333-Woah-007_S11_L001_R2_001.fastq.gz
Woah-996	homogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-996/SJ-2333-Woah-996_S9_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-996/SJ-2333-Woah-996_S9_L001_R2_001.fastq.gz
Woah-998	heterogametic	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-998/SJ-2333-Woah-998_S8_L001_R1_001.fastq.gz	/proj/snic2020-2-25/backup/raw_data/delivery02942/INBOX/SJ-2333/200110_A00181_0136_AHTWTCDSXX/Sample_SJ-2333-Woah-998/SJ-2333-Woah-998_S8_L001_R2_001.fastq.gz
```

## 4.1 **Run_3**: Selected sample run

```
snakemake -s workflow/findZX -j 15 -R all --configfile congig/Run_3.yml --cluster-config cluster/cluster_Wikstroemia.yml --cluster " sbatch -A snic2022-5-71 -t 8:00:00 -n 3 " -k --use-conda
```

config file for Run_3

The additional parameter ‘chr_highlight’ contains list of the contigs that are highlighted in the run.

```bash
## WITH WOAH-000 -- NOT INCLUDED IN THE FINAL RESULTS ###########################################################################################################################################
## ================================= ##
## findZX config file (test dataset) ##
## ================================= ##

# Variables marked with "[findZX]" or "[findZX-synteny]" are only used when deploying 
# the pipeline with either snakefile. Other variables are used for all analyses. 

threads_max: 18
mem_max: 16000
  # Specify maximum number of cores [threads_max] and memory [mem_max] allocation

# ============================ #
# Analysis name and input data #

run_name: Wikstroemia/scaf_assembly
  # Select an analysis name. Output files will be stored under "results/[run_name]"

units: config/Wikstroemia/scaf_assembly/Wikstroemia_units.tsv
  # Path to sample information file

ref_genome: /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa 
  # Path to study-species reference genome (not .gz format)

synteny_ref: 
  # [findZX-synteny] Path to synteny-species reference genome (not .gz format)

synteny_name: 
  # [findZX-synteny] Synteny-species name (can be any string, will be used for file and directory names)

# ================= #
# Plotting settings #

window_sizes: [50000, 100000, 1000000]
  # Choose genome window sizes for plotting (as many as you want)
  # Optimal sizes depend on reference genome fragmentation and size of the sex-linked region
  # Recommended sizes to start are: [50000, 100000, 1000000] (i.e. 50 kb, 100 kb, 1Mb)

chr_file: config/Wikstroemia/scaf_assembly/list_COI.txt
  # [findZX] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

chr_highlight: ["ptg000010l", "ptg000024l", "ptg000034l", "ptg000045l", "ptg000066l", "ptg000074l", "ptg000117l"]
  # [findZX] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

synteny_chr_file: None
  # [findZX-synteny] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

synteny_chr_highlight:
  # [findZX-synteny] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

# ================================== #
# Trimming and subsampling of reads  #

## These three variables control trimming and subsampling of reads
## Set all to "false" to disable trimming and subsampling
## Only one variable is allowed to be "true"

trim_reads: true 
  # Set to true for trimming of reads

trim_and_subsample: false
  # Set to true for trimming and subsampling of reads

subsample_only: false
  # Set to true for subsampling of reads (but not trimming)

subsample_basepairs: 1888226
  # Specify the total number of basepairs to extract from both fastq files
  # Will be used if [trim_and_subsample] or [subsample_basepairs] is set to "true"

  # Use this script to calculate expected coverage:
  # ./code/subsampling_cov_calv.sh <REF.fasta> <WANTED_COV> 

# ========================== #
# findZX-specific parameters # 
# ===== (edit if needed) ===== #
 
mismatch_settings: [0.0, 0.2]
  # Genome coverage results will be generated from the original BAM files ("unfiltered"),
  # and two other (modifiable) mismatches settings.
  # "0.0" = 0 mismatches allowed
  # "0.2" = <=2 mismatches allowed

minSizeScaffold: "7000"
  # The mimimum size of scaffolds in the reference genome to be included in the results

# ============================ #
# External software parameters # 
# ===== (edit if needed) ===== #

params:
  trimmomatic:
  # Control Trimmomatic settings here
    pe:
      trimmer:
        - "LEADING:20"
        - "TRAILING:20"
        - "SLIDINGWINDOW:4:20"
        - "MINLEN:120"
        - "ILLUMINACLIP:workflow/meta/adapters/TruSeq3-PE.fa:2:30:10"
```

```bash
##WITHOUT WOAH-000 ###########################################################################################################################################
## ================================= ##
## findZX config file (test dataset) ##
## ================================= ##

# Variables marked with "[findZX]" or "[findZX-synteny]" are only used when deploying 
# the pipeline with either snakefile. Other variables are used for all analyses. 

threads_max: 18
mem_max: 16000
  # Specify maximum number of cores [threads_max] and memory [mem_max] allocation

# ============================ #
# Analysis name and input data #

run_name: Wikstroemia/scaf_assembly
  # Select an analysis name. Output files will be stored under "results/[run_name]"

units: config/Wikstroemia/scaf_assembly/Wikstroemia_units_no_Woah000.tsv
  # Path to sample information file

ref_genome: /proj/snic2020-2-25/nobackup/ayushi/3_fasta_files/hifiasm_wikstroemia_hmc24s50.bp.p_ctg.fa 
  # Path to study-species reference genome (not .gz format)

synteny_ref: 
  # [findZX-synteny] Path to synteny-species reference genome (not .gz format)

synteny_name: 
  # [findZX-synteny] Synteny-species name (can be any string, will be used for file and directory names)

# ================= #
# Plotting settings #

window_sizes: [50000, 100000, 1000000]
  # Choose genome window sizes for plotting (as many as you want)
  # Optimal sizes depend on reference genome fragmentation and size of the sex-linked region
  # Recommended sizes to start are: [50000, 100000, 1000000] (i.e. 50 kb, 100 kb, 1Mb)

chr_file: config/Wikstroemia/scaf_assembly/list_COI.txt
  # [findZX] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

chr_highlight: ["ptg000010l", "ptg000024l", "ptg000034l", "ptg000045l", "ptg000066l", "ptg000074l", "ptg000117l"]
  # [findZX] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

synteny_chr_file: None
  # [findZX-synteny] Specify a file with list of chromosomes to only plot these (otherwise leave as "None")

synteny_chr_highlight:
  # [findZX-synteny] Specify chromosomes/scaffolds to highlight in plot type 4, or leave empty

# ================================== #
# Trimming and subsampling of reads  #

## These three variables control trimming and subsampling of reads
## Set all to "false" to disable trimming and subsampling
## Only one variable is allowed to be "true"

trim_reads: true 
  # Set to true for trimming of reads

trim_and_subsample: false
  # Set to true for trimming and subsampling of reads

subsample_only: false
  # Set to true for subsampling of reads (but not trimming)

subsample_basepairs: 1888226
  # Specify the total number of basepairs to extract from both fastq files
  # Will be used if [trim_and_subsample] or [subsample_basepairs] is set to "true"

  # Use this script to calculate expected coverage:
  # ./code/subsampling_cov_calv.sh <REF.fasta> <WANTED_COV> 

# ========================== #
# findZX-specific parameters # 
# ===== (edit if needed) ===== #
 
mismatch_settings: [0.0, 0.2]
  # Genome coverage results will be generated from the original BAM files ("unfiltered"),
  # and two other (modifiable) mismatches settings.
  # "0.0" = 0 mismatches allowed
  # "0.2" = <=2 mismatches allowed

minSizeScaffold: "7000"
  # The mimimum size of scaffolds in the reference genome to be included in the results

# ============================ #
# External software parameters # 
# ===== (edit if needed) ===== #

params:
  trimmomatic:
  # Control Trimmomatic settings here
    pe:
      trimmer:
        - "LEADING:20"
        - "TRAILING:20"
        - "SLIDINGWINDOW:4:20"
        - "MINLEN:120"
        - "ILLUMINACLIP:workflow/meta/adapters/TruSeq3-PE.fa:2:30:10"
```

There are 5 types of plots generated from findZX pipeline. 1) shows sex differences 2) shows values for every sex separately 3) shows heterozygosity and genome coverage values for every contig/chromosome/scaffold 4) shows heterozygosity and genome coverage values for each contig/chromosome/scaffold 5) shows heterozygosity and genome coverage profiles for every focused individual. 

All the plots are generated in the result folder with the outlier table.
