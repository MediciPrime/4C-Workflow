# Snakemake Workflow for 4C-ker

*4C-Workflow* is a snakemake implementation for the *4C-ker* library.
This generalized workflow aims to automate the 4C data analysis process allowing
you to get to your results faster. 


* [R.4Cker Github](https://github.com/rr1859/R.4Cker)
* [4C-ker Paper](http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004780)
* [Snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home)


## Setup

1. Clone the *4C-Workflow* directory into your local directory
2. Install [Miniconda](http://conda.pydata.org/miniconda.html)
3. Once installed use conda and the *environment.txt* file to create the conda environment
4. Activate that environment and start R and load `library(devtools)`
5. Once inside R enter `install_github("rr1859/R.4Cker")`

You are now ready to use *4C-Workflow* to perform your analysis.

## Analysis

1. Personalize the *config.yaml* file according to your experimental parameters
2. You may also notice the following files:
   * mock_vs_siAGO2.tsv
   * mock_vs_dsLamin.tsv
   * siAGO2_vs_dsLamin.tsv
3. You will need to create a *.tsv* for each of the experimental conditions you want to compare
4. Lastly add both a *.fasta* reference genome and primary enzyme sequence to the *reduced_genome* folder

**Note:** Make sure that your *.tsv* file names have *_vs_* between the two conditions you are comparing

## RUN!!!

1. Run the following command: `sh runscript`
2. Wait for the analysis to finish and find your results in <em>Output</em> folder


## View Results

Either IGV or UCSC Genome Browser can be used to view the .bedGraph and .bed files.

