# Wrapper for sam_bedGraph

*Sam_bedGraph* is a general purpose wrapper that will take any aligned sam file
and convert it into a bedGraph.  It's used to create the bedGraph files need to
identify the location containing 4C peeks. 

[Link to homepage](https://github.com/rr1859/R.4Cker)

## Input
* `sam`: aligned sam file containing the reads

## Output
* `bedGraph`: bedGraph file that contains all the called peaks

## Threads
Threads not supported.

## Params
None

## Example

```python
rule bedGraph_Counts:
    input:
		sam="sam_files/{sample}_aligned.sam"
    output: "sam_files/{sample}_aligned.bedGraph"
    wrapper:
        "/path/to/wrapper/location"
```
