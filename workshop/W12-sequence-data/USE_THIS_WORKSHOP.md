---
title: Workshop in markdown
author: Bin He, Leslie Speight
date: 2020-03-24
---

# Working with Sequence Data Tuesday Workshop 
Leslie Speight and Anna Ward

## Get data
Go to fastx or your laptop home directory (or wherever you prefer to store 
the workshop data, and clone the following repository

```bash
$ cd <path/to/workshop> # replace the second part with your directory name
$ git clone https://github.com/hezhaobin/2020-Data-Skills.git
$ cd 2020-Data-Skills/workshop/seqeunce-data/ # hint: Tab don’t type!
$ ls
```

Next open Rstudio and install required packages

```r
> if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
> BiocManager::install('qrqc')
# This is a Bioconductor package and helps us to visualize quality distribution across bases in reads
# Now you can minimize the RStudio window and get back to the terminal

```

## Trimming low-quality bases with sickle and seqtk

- Now, go back to your terminal and check to make sure that you are in the 'sequence-data'folder. Use 'ls' to see if the untreated1_chr4.fq is in the folder
- **TAB, DON'T TYPE**

```bash 
$ ./sickle se -f untreated1_chr4.fq -t sanger -o untreated1_chr4_sickle.fq
# Note you must have the “./” before “sickle” to make this work
# you should see the following output
#   Total FastQ records: 204355
#   FastQ records kept: 203121
#   FastQ records discarded: 1234
# sickle takes an input file through -f, a quality type through -t, and trimmed output with -o

# trimming low quality bases with seqtk
$ ./seqtk trimfq untreated1_chr4.fq > untreated1_chr4.trimfq.fq
# this takes a single argument and outputs trimmed sequence through standard out
```

## Comparing these results in R

Now, go back to RStudio. Be sure to check that your working directory is pointed to `sequence-data`. You can check the current working directory by `getwd()` and set the working directory by either `setwd(path/to/working/directory)` or use the file pane to navigate to the folder and use the "more" icon -> set as working directory

```r
> library(qrqc) 

# if you get an error, most likely you didn’t install the package successfully
                # go back to the beginning of the workshop and install `qrqc`
                
# Specify the file names using a character vector “c()”

> fqfiles <- c(none= "untreated1_chr4.fq", sickle= "untreated1_chr4_sickle.fq", trimfq= "untreated1_chr4.trimfq.fq")

# Below we use the `lapply()` function to apply a custom function, `readSeqFile()` to read the content of each of the three files
# We are loading each file in, and using qrqc’s readSeqFile

> seq_info <- lapply( fqfiles, function(file) { readSeqFile(file, hash=FALSE, kmer=FALSE) } )

# We only need qualities, se we turn off some of readSeqFiles’s other features
# note, if you type all 5 lines of code on the same line, you need a semi-colon “;” after each of line 2 and 3.

> quals <- mapply(function(sfq, name){
                    qs <- getQual(sfq)
                    qs$trimmer <- name
                    qs
        }, seq_info, names(fqfiles), SIMPLIFY=FALSE)   
#here we are extracting qualities as dataframe, and appending a column of which trimmer is used
# this is used in later plots 

> d <- do.call(rbind,quals)
# this combines separate dataframes in a list into single dataframes

> p1 <- ggplot(d)+geom_line(aes(x=position,y=mean,linetype=trimmer))
> p1 <- p1 + ylab("mean quality (sanger)")+ theme_bw()
> print(p1)
# visualize qualities 

# what do you see? 

> p2 <- qualPlot(seq_info, quartile.color=NULL, mean.color=NULL) + theme_bw()
> p2 <- p2 + ylab("quality (sanger)")
> print(p2)
# uses qrqc’s qualPlot with list produces panel plots
# only shows 10% to 90% quantiles and lowes curve
```

## Counting Nucleotides
- We will be working in the terminal for this section
- We will use the Python implementation of `readfq`, `readfq.py`
    - `readfq()` takes a file object (a filename that has been opened) and will generate FASTA/FASTQ entried until it reaches the enc of the file or stream.
    - Each FASTA/FASTQ entry is returned by `readfq()` containing the entry's description, sequence, and quality.

```bash
$ cd Desktop/Seq_data/2020-Data-Skills/workshop/seqeunce-data/
$ ls  
# we are now gonna look at and use a simple program that counts the number of each IUPAC nucleotide in a file
$ vim nuccount.py
```

```python
#!/usr/bin/env python
# nuccount.py -- tally nucleotides in a file
import sys
from collections import Counter
from readfq import readfq
# Use Counter in the collections module to count nucleotides. Works much like Python dict

IUPAC_BASES = "ACGTRYSWKMBDHVN-."
# Defines all IUPAC nucleotides (standardized and unstandardized) , which we will use later to print bases in a consistent order. 

# initialize counter
counts = Counter()
# create a new Counter object

for name, seq, qual in readfq(sys.stdin):
    # uses the readfq function in the readfq module to read FASTQ/FASTA entries from the file handle argument into the name, seq, qual variables
    # for each sequence entry, add all its bases to the counter
    counts.update(seq.upper())
    # takes any iterable object and adds them to the counter. We convert all characters to uppercase with the upper( ) method, so that       lowercase soft-masked bases are also counted. 

# print the results
for base in IUPAC_BASES:
    print (base + "\t" + str(counts[base]))
# we iterate over all IUPAC bases and print their counts
```

- Now run the script
    ```bash
    $ cat contam.fastq | ./nuccount.py
    ```

## Indexed FASTA Files
- Make sure you are working within /sequence-data/
- Extracting numerous random subsequences from a FASTA file can be quite computationally costly.
- A common computational strategy that allows for easy and fast random access is indexing the file. It allows us to easily extract subsequences

```bash
$ gunzip Mus_musculus.GRCm38.75.dna.chromosome.8.fa.gz
$ samtools faidx Mus_musculus.GRCm38.75.dna.chromosome.8.fa 8:123407082-123410744
```

- What do you see?
- `samtools faidx` allows for multiple regions at once, so we could do:

    ```bash
    $ samtools faidx Mus_musculus.GRCm38.75.dna.chromosome.8.fa \ 8:123407082-123410744 8:123518835-123536649
    ```
