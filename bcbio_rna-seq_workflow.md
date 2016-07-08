# Bcbio workflow
**Documentation for bcbio:** [bcbio-nextgen readthedocs](http://bcbio-nextgen.readthedocs.org/en/latest/contents/pipelines.html#rna-seq)

## Set-up
1. One time only follow set-up instructions for [Rory's bcbio.rnaseq](https://github.com/roryk/bcbio.rnaseq): 
	- Install [lein](https://github.com/technomancy/leiningen) - I installed lein in ~/bin
	- Add lein location to path in `~/.bashrc`:
		- `export PATH=~/bin:$PATH`
	- I could not get pandoc installed
2. Make directory structure 
    - `cd path-to-consult-folder`
    - `mkdir analysis meta config data`
    
3. Download fastq files from facility to data folder
	
	- Download fastq files from a non-password protected url
		- `wget --mirror url` (for each file of sample in each lane)
   	 	- Rory's code to concatenate files for the same samples on multiple lanes: 
    
    			barcodes="BC1 BC2 BC3 BC4"
    			for barcode in $barcodes
    			do
    			find folder -name $barcode_*R1.fastq.gz -exec cat {} \; > data/${barcode}_R1.fastq.gz
    			find folder -name $barcode_*R2.fastq.gz -exec cat {} \; > data/${barcode}_R2.fastq.gz
    			done

   	- Download fastq files from BioPolymers: 
		- `sftp username@bpfngs.med.harvard.edu`
		- `cd` to correct folder
		- `mget *.tab`
		- `mget *.bz2`



4. Settings for bcbio- make sure you have following settings in `~/.bashrc` file:
 
 ```
 unset PYTHONHOME
 unset PYTHONPATH
 module load dev/java/jdk1.7
 module load stats/R/3.2.1
 module load dev/perl/5.18.1
 export PATH=/opt/bcbio/centos/bin:$PATH
 ```
    
5. Within the `meta` folder, add your comma-separated metadata file (`projectname_rnaseq.csv`)
	- first column is `samplename` and is the names of the fastq files as they appear in the directory
	- second column is `description` and is unique names to call samples (can be the file name without the extension (.fastq or R#.fastq for paired-end reads))
	- column entitled `samplegroup` is your sample groups
	- **FOR CHIP-SEQ** need additional columns:
		- `phenotype`: `chip` or `input` for each sample
		- `batch`: batch1, batch2, batch3, ... for grouping each input with it's appropriate chip(s)
	- additional specifics regarding the metadata file: [http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration](http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration) 
        
6. Within the `config` folder, add your custom Illumina template
    - Example template for human RNA-seq using Illumina prepared samples (genome_build for mouse = mm10):

	```
        details:
          - analysis: RNA-seq
            genome_build: hg19
            algorithm:
              aligner: star
              quality_format: Standard
              trim_reads: read_through
              strandedness: firststrand 
        upload:
          dir: ../results
        star-illumina-rnaseq.yaml 
```
	- Additional parameters can be found: [http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration](http://bcbio-nextgen.readthedocs.org/en/latest/contents/configuration.html#automated-sample-configuration) 
	- Best practice templates can be found: [https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates](https://github.com/chapmanb/bcbio-nextgen/tree/master/config/templates)

 
7. Within the `data` folder, add all your fastq files to analyze.

##Analysis

1. Go to analysis folder and create the full Illumina instructions using the Illumina template created in Set-up: step #6.
    - `bsub -Is -n 4 -q interactive bash` start interactive job
    - `cd path-to-folder/*_RNAseq/analysis` change directories to analysis folder
    - `bcbio_nextgen.py -w template ../config/star-illumina-rnaseq.yaml ../meta/*-rnaseq.csv ../data/*fastq.gz` run command to create the full yam file

2. Create script for running the job (in analysis folder)

	```
	#!/bin/sh
	#BSUB -q priority
	#BSUB -J *-rnaseq
	#BSUB -o *-rnaseq.out
	#BSUB -N
	#BSUB -u "email@hsph.harvard.edu"
	#BSUB -n 1
	#BSUB -R "rusage[mem=8024]"
	#BSUB -W 50:00
	#
	date

	bcbio_nextgen.py ../config/*-rnaseq.yaml -n 64 -t ipython -s lsf -q parallel '-rW=90:00' -r mincores=2 -rminconcores=2 --retries 3 --timeout 580

	date
	```

3. Go to work folder and start the job - make sure in an interactive session 

	```
cd path-to-folder/*-rnaseq/analysis/*-rnaseq/work
bsub < ../../runJob-*-rnaseq.sh
```

### Exploration of region of interest

1. The bam files will be located here: `path-to-folder/*-rnaseq/analysis/*-rnaseq/work/align/SAMPLENAME/NAME_*-rnaseq_star/`

2. Extracting interesting region (example)
	- `samtools view -h -b  sample1.bam "chr2:176927474-177089906" > sample1_hox.bam`

	- `samtools index sample1_hox.bam`

## Report generation
1. Report creation and creating project_summary.csv

```
source ~/.bashrc
cd ~/bcbio.rnaseq
bsub -Is -q interactive bash
lein run summarize path-to-project-summary.yaml -f "~batch+panel"
```
2. Copy to local computer the results/*-rnaseq/ folder and the results/*-rnaseq/summary/qc-summary.Rmd
    - `scp -r username@orchestra.med.harvard.edu:path-to-folder/*-rnaseq/analysis/*-rnaseq/results/
date_*-rnaseq/ .`

3. Within R Studio:
	- load library(knitrBootstrap)
	- three dashes at top and bottom of knitrBootstrap specifics
	- Copy over header info for knitrBootstrap
	- Alter paths to files
