# Snakemake Workflow for 4C-ker

*4C-Workflow* is a snakemake implementation for the *4C-ker* library.
This generalized workflow aims to automate the 4C data analysis process allowing
you to get to your results faster. 

<ul>
<li><a href="https://github.com/rr1859/R.4Cker" target=_blank>R.4Cker Github</a></li>
<li><a href="http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004780" target=_blank>4C-ker Paper</a></li>
<li><a href="https://bitbucket.org/snakemake/snakemake/wiki/Home" target=_blank>Snakemake</a></li>
</ul>

## Setup
<ol>
<li>Clone the <em>4C-Workflow</em> directory to your local folder</li>
<li>Install <a href="http://conda.pydata.org/miniconda.html" target=_blank>Miniconda</a></li>
<li>Once installed use conda and the <em>environment.txt</em> file to create the conda environment</li>
<li>Activate that environment and start R and load <code>library(devtools)</code></li>
<li>Once inside R enter <code>install_github("rr1859/R.4Cker")</code></li>
</ol>

You are now ready to use <em>4C-Workflow</em> to perform your analysis.

## Analysis
<ol>
<li>Personalize the <em>config.yaml</em> file according to your experimental parameters</li>
<li>You may also notice the following files:</li>
<ul>
<li><em>mock_vs_siAGO2.tsv</em></li>
<li><em>mock_vs_dsLamin.tsv</em></li>
<li><em>siAGO2_vs_dsLamin.tsv</em></li>
</ul>
<li>You will need to create a <em>.tsv</em> for each of the experimental conditions you want to compare</li>
<li>Lastly add both a <em>.fasta</em> reference genome and primary enzyme sequence to the <em>reduced_genome</em> folder</li>
</ol>

<b>Note:</b> Make sure that your <em>.tsv</em> file names have <em>"_vs_"</em> between the two conditions you are comparing

## RUN!!!
<ul>
<li>Run the following command: <code>sh runscript</code></li>
<li>Wait for the analysis to finish and find your results in <em>Output</em> folder</li>
</ul>

## View Results
<p>
Either IGV or UCSC Genome Browser can be used to view the .bedGraph and .bed files.
</p>
