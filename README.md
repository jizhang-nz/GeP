# Genome Profiler (GeP)
Genome Profiler (GeP) is a program to perform whole-genome multilocus sequence typing (wgMLST) analysis for bacterial isolates using gene-by-gene allele-calling strategy.

This program is designed to run in Unix-like operating systems such as Mac OS X or Linux.

GeP uses conserved gene neighborhoods (CGN) to resolve gene paralogy. This approach allows the differentiation of orthologs from recently duplicated paralogs, which are often indistinguishable by pairwise sequence alignment.

GeP is best for resolving relationship of closely related isolates, such as isolates of the same ST (e.g. in an outbreak investigation). For distantly related isolates, one could consider using  [Fast-GeP](https://github.com/jizhang-nz/fast-GeP), as the running time of GeP would increase dramatically due to frequent calling of BLASTX to find alleles of greater genetic variation.

## Prerequisites
Before start, you need to make sure the following three programs were fully functional in your system:
   * [BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
   * [MAFFT](https://mafft.cbrc.jp/alignment/software/)
   
## Usage
Let's assume you have put the `GeP.pl` file in your PATH. If not, or if you prefer, you could always put the `GeP.pl` file along with your other input files, and use commands like:

    perl GeP.pl -g list.fas.txt -r reference.gbk

or something like:

    GeP.pl -g list.fas.txt -r reference.gbk

For full list of switches and options, please use commond like:

    GeP.pl -h


## Citation
Please see our publication: 
   * [Refinement of Whole-Genome Multilocus Sequence Typing Analysis by Addressing Gene Paralogy](http://jcm.asm.org/content/53/5/1765.abstract)
