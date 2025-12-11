# GenomicRanges & IRanges – R Practice

This repository contains two R practice scripts exploring the **IRanges**, **GenomicRanges**, and **Biostrings** packages from the Bioconductor ecosystem.  
These exercises focus on genomic interval manipulation, biological sequence handling, and fundamental operations commonly used in computational genomics and bioinformatics.

The goal of this mini-project is to build confidence with core Bioconductor data structures and functions.

---

## 📄 Files Overview

### 1. `IRanges_introduction.R`

This script introduces the **IRanges** package and demonstrates how to work with genomic intervals.

Main topics covered:

- Creating `IRanges` objects (`start`, `end`, `width`)
- Accessing and modifying interval attributes
- Naming and subsetting ranges
- Combining multiple IRanges objects
- Interval normalization functions:
  - `reduce()`
  - `disjoin()`
  - `gaps()`
- Custom visualization of IRanges
- Inter-range vs. intra-range transformations:
  - `shift()`
  - `narrow()`
  - `resize()`
  - `flank()`
  - `restrict()`

---

### 2. `BioStrings_GenomicRanges.R`

This script introduces the **Biostrings** and **GenomicRanges** packages.

#### Biostrings section

- Working with DNA, RNA, and protein sequences:
  - `DNAString`, `RNAString`, `AAString`, `StringSet`
- Using IUPAC ambiguity codes
- Sequence transformations:
  - `reverse()`
  - `complement()`
  - `reverseComplement()`
  - `translate()`
- Frequency analyses:
  - `letterFrequency()`
  - `dinucleotideFrequency()`
  - `trinucleotideFrequency()`
  - `oligonucleotideFrequency()`
- Sliding-window frequency analysis
- Generating consensus matrices from aligned sequences

#### GenomicRanges section

- Creating `GRanges` objects
- Exploring genomic interval metadata:
  - `seqnames`
  - `strand`
  - `ranges`
  - `seqlengths()`
  - `genome()`
  - `isCircular()`
- Introduction to RLE (Run-Length Encoding) concepts

---

## 🎯 Skills Developed

### Genomic Data Structures
- Understanding and manipulating `IRanges` and `GRanges`
- Performing inter- and intra-range transformations
- Working with genomic coordinates and strand information

### Biological Sequence Handling
- Representing DNA, RNA, and amino acid sequences
- Computing sequence frequencies and performing translations
- Handling ambiguous nucleotide symbols

### Data Transformation & Visualization
- Plotting genomic intervals with a custom visualization function
- Identifying overlaps, reduced ranges, and gaps

### Bioconductor Ecosystem Familiarity
- Installing and using `IRanges`, `GenomicRanges`, and `Biostrings`
- Understanding their role in genomic data analysis workflows

---

## ▶️ How to Run the Scripts

### Install dependencies

```r
install.packages("BiocManager")
BiocManager::install(c("IRanges", "GenomicRanges", "Biostrings"), dependencies = TRUE)
library(IRanges)
library(GenomicRanges)
library(Biostrings)
```

### Purpose of the Repository

This repository serves as:
- Personal practice notebook
- Personal reference 
