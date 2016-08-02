# Wrapper for 4C_purify

*4C_purify* is a general purpose wrapper that will take a fasta, bed and primer
sequence and identify the bait location in the bed file. This is done because
the bait region contains a disproportionate number or reads thus amplifying the
bait signal. 

Note: *primer* should only contain the primer sequence

[Link to homepage](https://github.com/rr1859/R.4Cker)

## Input
* `fasta`: fasta file containing the flaking sequence for the restriction enzyme
* `bed`: bed file conversion of the fasta from file from above

## Output
* `temp.bed`: temp bed file containing the sequence for the bait region and the
regions adjacent to it

## Threads
Threads not supported.

## Params
* `primer`: primer sequence w/out barcode

## Example

```python
rule purify_helper:
    input:
		fasta='reduced_genome/dm6_hindiii_flanking_sequences_37_unique_2.fa',
		bed='reduced_genome/dm6_hindiii_flanking_sites_37_unique_2.bed'
    output:
        temp('temp.bed')
    params: 
		primer='ACCGTTTTTGACAACAGCAGCTGTA'
    wrapper:
        "/path/to/wrapper/location"
```
