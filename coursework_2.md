---
title: "Coursework 2"
author: Jennifer J. Stiens
output:
  html_document:
    keep_md: yes
    pdf_document: default
  header-includes:
  - \usepackage{color}
---
## Jennifer Stiens
j.j.stiens@gmail.com

[github.com/jenjane118/ngs_tutorials/coursework_1.md](https://github.com/jenjane118/ngs_tutorials)

##Question 1
$\color{red}{\text{Compile a set of virulence-related genes in the AFPN02.1 strain and other E. coli strains and compare them}}$

First I went to the VFPB website (http://www.mgc.ac.cn/VFs/main.htm) and looked up 'virulence factors' associated with E. coli strains. The most common of these are: adherence, toxin, invasion, Type III translocated protein, and immune evasion. I can use these terms to annotate genes associated with virulence. 

[VFDB website, Escherichia Coli](http://www.mgc.ac.cn/cgi-bin/VFs/genus.cgi?Genus=Escherichia)
\
\
$\color{red}{\text{a) Assemble a set of virulence-associated genes from public source or publication}}$
 \
 \
 
From VFDB I was able to download a database of DNA sequences of core set of experimentally-verified virulence-associated genes. This will be used to extract genes (and relevant sequences) for the associated adhesion genes for searching and annotating the AFPN02.1 genome. (VFDB_setA_nt.fas)
(http://www.mgc.ac.cn/VFs/download.htm) (Jin, *et.al.*, 2007). 

Using the information on the VFDB website (http://www.mgc.ac.cn/cgi-bin/VFs/genus.cgi?Genus=Escherichia), Brzuszkiewicz, *et. al.* (2011), and the spreadsheet of virulence genes from EHEC and EAEC strains from Cheung, *et.al*, (2011), I compiled a list of 18 pertinent virulence-related genes. These are mostly related to adhesion and toxin production. I used these to select virulence-related genes from the VFDB.
![virulence spreadsheet](spreadsheet.png){width=75%}
[Brzuszkiewicz, et al, 2011]

\
\


```bash
## help for this code came from (https://infoplatter.wordpress.com/2013/10/15/extracting-specific-fasta-records-from-a-multi-fasta-file/)
## grabs all ecoli files and sequences
awk 'BEGIN {RS=">"}/Escherichia/{print">"$0}' VFDB_setA_nt.fas > ecoli_vfdb_genes.fas
## grabs all files with these genes
cat ecoli_vfdb_genes.fas | awk 'BEGIN {RS=">"}/aggR|agg3B|aat|astA|hlyA|sepA|aggA|setA|espP|setC|aat|fimH|fimA|cpxA|cpxR|csgA|elfA|stx/{print">"$0}' > virulence_vfdb_genes.fa
```

These have really elaborate headings, so I needed to make them much more simple:

```bash
## removes part of heading before gene name
cat virulence_vfdb_genes.fa | sed 's/^>.*(.*)\s(//g' > temp.fa
## removes rest of heading after gene name
cat temp.fa | sed 's/)\s.*$//g' > virulence_vfdb_editedgenes.fa
mv virulence_vfdb_editedgenes.fa vfdb_virulence_genes.fa
```

\
\

$\color{red}{\text{b) Build a separate set of virulence-associated genes in the annotation file created for AFPN02.1}}$

The next step is to work with the annotation.gff file to retrieve annotations based on 'virulence', 'adherence', 'toxin', and 'Type III' (short for Type III aggregative adherence fimbriae) and convert into .bed format.



```bash
cat annotation.gff | grep -E 'virulence|adherence|toxin|invasion|Type III' | awk 'BEGIN {FS="\t"}  split($9, captured, /[(=);]/) >=10  {print "sequence1" "\t" $4 "\t" $5 "\t" captured[10] "\t" captured[4] "\t" $7}' > present_in_AFPN02_virulence_genes.bed
```

Then we have to retrieve the DNA sequence for these genes and save in fasta format:


```bash
/s/software/bedtools/v2.27.1/bin/bedtools getfasta -name -s -fi ${st_path}/results_GC/annotation/genome.fna -bed present_in_AFPN02_virulence_genes.bed -fo present_in_AFPN02_virulence_genes.fasta
```

A bit more formatting to get rid of extra elements in the heading:


```bash
cat present_in_AFPN02_virulence_genes.fasta | sed 's/(.*//g' > present_in_AFRN02_virulence_genes.fasta
```


Then we have to combine the two files:

```bash
cat present_in_AFPN02_virulence_genes.fasta vfdb_virulence_genes.fa > final_comparison_virulence.fasta
cp present_in_AFPN02_virulence_genes.fasta vfdb_virulence_genes.fa final_comparison_virulence.fasta ${st_path}/results_GC/wholeGenomeExamples

```
\
\

$\color{red}{\text{c) Use BRIG to visualise which of the virulence genes are present/absent in E. coli strains}}$
\
\

To run brig, use the following command:

```bash
/s/software/brig/BRIG-0.95-dist/brig.sh
```
<p>
![Brig analysis](final_comparison_virulence.fasta.png)
</p>
Most of the virulence genes were present in AFPN02.1. However, hylA was only identified with a short sequence, aggA was completely missing, esp was mostly missing, though some parts were  present at lower identity, fimH was completely missing, stxA was more weakly identified, and agg3B was more weakly identified. The other shiga-toxin genes, stx2A/B were present in the AFPN02.1 strain, as well as in the sakai strain, which is a EHEC (enterohemmoragic) strain, but isn't present at all in the other strains, except for short regions of identity with the other EHEC strain. The EHEC strain seemed the most different from the other strains, with only short regions of identity spread out over the genome, however, it seemed there were hits in all of the genes. 


## Question 2


$\color{red}{\text{Select ONE of the virulence genes (or if you prefer one operon) present in AFPN02.1
and study this}}$ $\color{red}{\text{gene/operon in more detail including its biological action/mechanism phylogeny.}}$




## References

Brzuszkiewicz, E., Thürmer, A., Schuldes, J., Leimbach, A., Liesegang, H., Meyer, F.-D., … Daniel, R. (2011). Genome sequence analyses of two isolates from the recent Escherichia coli outbreak in Germany reveal the emergence of a new pathotype: Entero-Aggregative-Haemorrhagic Escherichia coli (EAHEC). Archives of Microbiology, 193(12), 883–891. https://doi.org/10.1007/s00203-011-0725-6

Cheung, M., Li, L., Nong, W., & Kwan, H. (2011). 2011 German Escherichia coli O104:H4 outbreak: Whole-genome phylogeny without alignment. BMC Research Notes, 4(December). https://doi.org/10.1186/1756-0500-4-533

Jin, Q., Sun, L., Yang, J., Chen, L., & Yu, J. (2007). VFDB 2008 release: an enhanced web-based resource for comparative pathogenomics. Nucleic Acids Research, 36(Database), D539–D542. https://doi.org/10.1093/nar/gkm951