# GEMAP_NCLDV
**GE**nome **M**apping **A**nalysis **P**ipeline for giant viruses (but can also be used for other types of genomes as well. This was designed for use with timescale data (eg. mapping reads to a genome over multiple days/months).

## Overview of Use

1) Perform mapping to your desired genome using your mapping tool of choice but make sure you use CoverM for the coverage information. (a mapping script uing bbmap has been provided)
2) Run post_map_process.py to organize all of your read outputs into a single file.
3) Do genome annotations on the genes of your genome (GVOG HMM database and hmmscan was used in this example).
    *If using Hmmscan, you will have to do some additional formatting that will be described later. 
4) Run genome_annotation_many.py to map the hits from your hmmsearch to the functions from the annotted database.
5) Run genome_clustermaps.py to get figures of your data as well as summary csv files. 

## 1. Mapping pipeline

### Dependencies
1) python>= 3
2) bbmap
3) samtools
4) CoverM


The provided mappping pipeline script maps using bbmap with a minID of 97%. If you would like to change this parameter, simply go into the script and edit line 20, chinging 0.97 to your desired number. 

### Inputs
1) input directory of reads in either fastq or fastq.gz format.
2) Output directory
3) Reference genome you are mapping to

Example
`python Mapping_Pipeline_onerefgenome_97minID.py reads Outdir ecoli_genome.fna`

### Outputs
1) Sorted Bam files for all reads (these can be quite large if you have a long time period so I reccommend making a script to delete them)
2) Coverm output files for each read showing RPKM and count values for each gene. 

## 2. Post map processing

### Dependencies
1) python
2) numpy
3) pandas

This script will take inputs of CoverM files and rename them, sort them by date, and average days if you have multiple samples per date. It will combine all reads into one file. 

**IMPORTANT** Doing this will remove the first column (the column containing gene names) from your spreadsheet. You can easily just paste it back from one of the coverM files as the order doesn't change. 

### Inputs
-c: folder with all CoverM files in it
-r: an empty folder called 'renamed'
-ref: reference file for renaming your dates from read name to date name (eg. SRR77893 -> 11/19/2020 ). **NOTE** the read names must be in a column called SRR_run and the date must be in a column called Date
-com: an empty folder called 'combined'
-o: Output filename (can be anything really)

### example
`python Post_mapping_processing.py -c Coverm -r Renamed -ref reffile.csv -com Combined -o Genome_out.csv`

### Outputs
 1) a few useless folders containing renamed files and combined files. (These can be deleted)
 2) A file with all of your reads in order of date and averaged if there were multiples for each day. 

## 3. Genome annotations

### Dependencies
1) hmmer

### Steps
1) Download the hmm database of your choice (I used [GVDB](https://zenodo.org/record/4728209/files/GVOGs.tar.gz?download=1)).
2) Follow instructions from [HMMER](https://www.mankier.com/package/hmmer) to build an hmm and then scan your genomes against it.
3) Parse the output file to only include best hits with this line of code
`awk '!x[$3]++' MYOUTFILE.pfam  > MYBESTHITS.pfam`
4) You will still have to do some editing of the resulting spreadsheets in order to get them ready. The final product should be a 2 column spreadsheet wth one  column being the 'query' and the other being the 'Hit'. 

## 4. Using Genome_annotation_many.py

Inputs
1) Directory with all your formatted spreadsheets (space separated) from above step (4).
2) An annotation tsv file (this comes with most hmms that you download and it has functional information). 
For the GVDB, this file is called gvog.complete.annot.tsv.

Usage
`python genome_annotation_many.py input_directory gvog.complete.annot.tsv`

Output
1) This will map the functions to your annotation files and save them as .tsv files

## 5. Generating figures

### Dependencies
1) Seaborn
2) Matplotlib

This script will generate all sorts of valuable output files showing expression of different genes during the duration of your study.

### Inputs
1)-i: Your combined and renamed reads file from part 2. (CSV)
2) -a:Your annotation file from part 4. (TSV)
3) -c:COG categories that you want to graph. [List of COG categories](http://clovr.org/docs/clusters-of-orthologous-groups-cogs/)

### Usage
`python genome_clustermaps.py -i reads.csv -a annotated_file.tsv -c J,K,L,M,A,C`

### Outputs
1) annotated.csv : a csv that combines your counts and annotation data into one sheet
2) .png files of heatmaps for each COG functional category
![Heatmap!](/images/61D.42_rpkm_avgcoenzyme.png)

4) Heatmap of overall expression
5) nog_categories.csv: Number of genes in each COG category.
6) top50.csv: list of the top 50 most expressed genes and their average RPKM value
7) sums.png: A graph of expression over the days of the total genome.
![Sums!](/images/61D.42_rpkm_avg_sums.png)
