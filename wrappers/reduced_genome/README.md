# Wrapper for reduced_Genome

Reduced_Genome is a general purpose wrapper that will take any reference genome
and extract only the user defined relevant segments. It is specifically made for
the 4C-ker algorithm since the algorithm requires a reduced genome. Three params
*sequence fragment length*, *restriction enzyme (RE)* and *reference genome* are
required to get the algorithm to work. 

Note: *sequence fragment length* include the sizes of the barcode, primer and RE

[Link to homepage](https://github.com/rr1859/R.4Cker)

## Input
* `enzyme`: fasta file containing the sequence of the restriction enzyme
* `genome`: fasta file of the reference genome used in the experiment
* `fragLen`: although not an actual file, fragLen should represent the correct
length of the working fragment (see **params**)

## Output
* `reducedGenome`: a total of 13 reduced genome files are created for downstream
analysis

## Threads
Threads not supported.

## Params
* `fragLen`: *important* to change based on the size of the barcode, primer and
restriction enzyme length
* `enzyme`: name of the enzyme being used ex. hindiii
* `genome`: name of reference genome being use ex. dm6

## Example

```python
rule reduced_genome:
    input:
		enzyme = "reduced_genome/hindiii.fa",
        genome = "reduced_genome/dm6.fa"
    output:
        "reduced_genome/dm6.fa.fai"
    params: 
		extra = "",
		fragLen = 37,
		enzyme = 'hindiii',
		genome = 'dm6'
    wrapper:
        "/path/to/wrapper/location"
```
