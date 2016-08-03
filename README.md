# Snakemake Workflow for 4C-ker

*4C-Workflow* is a snakemake implementation for the *4C-ker* library.
This generalized workflow aims to automate the 4C data analysis process allowing
you to get to your results faster. 

This workflow depends heavily upon *Miniconda*, *R.4Cker*, and *Snakemake* so if you
want know what they are checkout the links below. Otherwise just follow the **setup**
process and you should be able to get the *4C-Workflow* up and running.

* [Miniconda](http://conda.pydata.org/miniconda.html)
* [R.4Cker Github](https://github.com/rr1859/R.4Cker)
* [4C-ker Paper](http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004780)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)


## Setup

1. Clone the *4C-Workflow* directory into your local directory
2. Install [Miniconda](http://conda.pydata.org/miniconda.html)
3. Once installed use the terminal to create a conda while inside the *4C-Workflow* folder
   - `conda create --name 4C-Workflow --file environment.txt`
4. Inside terminal activate the environment, start R and install *R.4Cker* using the following commands:
   - `source activate 4C-Workflow`
   - `R`
   - `library(devtools)`
   - `install_github("rr1859/R.4Cker")`

You are now ready to use *4C-Workflow* to analyze your raw 4C data.

## Analysis

1. *Snakemake* uses the *config.yaml* file to understand the experimental parameters, personalize it accordingly
   - The file is divided into three general blocks:
	 - **comparisons**
	   - Raw bedGraphs for each of the experimental comparisons are placed here
	   - Make sure to follow the exact formatting displayed in the example
		 - **Note:** Each of the conditions in the example have two replicates
		 - *bedGraph/* is the folder that *Snakemake* will place the files into
	 - **samples**
	   - *Snakemake* will use the path name to find all the raw *.fastq* files
	   - Make sure to insert the **full path** to your *.fastq* files
	 - **other**
	   - *bait_chr*: short form of the chromosome name
	   - *bait_coord*: numerical start site for your primer, usually after restriction enzyme site
	   - *bait_name*: the shortform bait name
	   - *primary_enz*: sequence of the primary enzyme
	   - *species*: shortform species name
	   - *reduced_genome*: shortform name of reference genome
	   - *fragment_len*: barcode + primer + primary restriction enzyme, eg. 37
	   - *primer*: genome sequence for the primer
2. You may also notice the following files:
   * mock_vs_siAGO2.tsv
   * mock_vs_dsLamin.tsv
   * siAGO2_vs_dsLamin.tsv
   - A *.tsv* file is created for each of the experimental conditions specified in the **comparisons** section from above
4. Lastly add both a *.fasta* reference genome and primary enzyme sequence to the *reduced_genome* folder
   - **Make Sure** the names for both the reference genome and enzyme match the names provided for *primary_enz_name* and *reduced_genome* from **other**
	 - eg. dm6 == dm6.fasta & hindiii == hindiii.fasta

**Note:** Make sure that your *.tsv* file names have *_vs_* between the two conditions you are comparing

## RUN!!!

1. Finally open the terminal and while inside the *4C-Workflow/* folder run the following command:
   - `sh runscript`
2. Wait for the analysis to finish and find your results in *Output/* folder


## View Results

Either IGV or UCSC Genome Browser can be used to view the .bedGraph and .bed files.

