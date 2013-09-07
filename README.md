AGEM: Across-GEnome Matching
****************************

AGEM is a python package that 
1) pipe sequence alignment results ( or filtered blast results, or any other data that represent homologous relationships) into network analysis; 
2) implement several maximum bipartite matching algorithms
3) plot diagnostic curves to choose appropriate matching

Tutorial
********

1. Input Data
==========
Assuming you want to match across genome A and genome B. The first and foremost step is to align nucleotide/protein sequences of the two genomes against each other. The second step is to prune sequence alignment results into data that represent homologous relationships of genes/proteins between genome A and genome B. If you don't know how to do the abovementioned steps, please read "1.2 Obtain filtered blast results" below.



1.1 Input data format
----------------------

Now we assume you've already had filtered sequence alignment results between genome A and genome B. What worths noting is aligning genome A against genome B is different from aligning genome B against genome A. Reciprocal blast results are always recommendated, even unilateral blast results can also be run by AGEM. For example, the reciprocal blast results between genome A and genome B are in two files: A2B.dat (aligning genome A against genome B) and B2A.dat (aligning genome B against genome A). The two input files should look like:

format example of A2B.dat

genome.A_gene.i [tab-delimited] genome.B_gene.j [tab-delimited] homologys_score_between gene i and j

...                   ...                    ...

format example of B2A.dat

genome.B_gene.h [tab-delimited] genome.A_gene.k [tab-delimited] homologys_score_between gene h and k

...                   ...                    ...


In summary, the format has three columns: genes in one genome, genes in the other genome, homology score between two genes. If you don't know how to generate such files and file format, please follow "1.2 Obtain filtered blast results" below.

1.2 Obtain filtered blast results
-----------------------------
There could be many ways, but here the typical alignment tool NCBI Blast is taken as an example. 

1. Install standalone NCBI Blast so you can command run every NCBI Blast function on a server or a local machine. The most useful functions include building blast database from sequence data you generated or downloaded, run all kinds of sequence alignment and select the output file format. To install and use standalone NCBI Blast, please read NCBI help book: http://www.ncbi.nlm.nih.gov/books/NBK52640/
2. Build blast database of genome A and B; choose appropriate command to run blast according to your sequence type; when you run blast, do specify the output format in output format option. For more details on output format option, type blast command with help option(e.g. blastp -help).  Output format option 6 is recommendated. 
3. If you choose output format option 6, the results you get should be table-like. Some columns of the result table are properties such as alignment length and percent identity with which you can prune the alignment results. 
4. From the filtered result table, you can grep column 1(gene ids of genome A), column 2(gene ids of genome B) and the last column (bit score) to format the input file. 

If you use linux/unix or Mac terminal, you can implement step 3 and 4 by typing some awk commands as shown below.

awk commands example:

awk '{if ($3>=50 && $4>=200 && $12>=200) print $0}' rat2mouse.aln > rat2mouse.aln.flt

comment: this command filters blast output (format option 6) in the file "rat2mouse.aln" by column 3, 4, and 12 and directs the filtered output into the new file "rat2mouse.aln.flt". 

awk "{print $1, $2, $14}" rat2mouse.aln.flt >rat2mouse.dat

comment: this command greps column 1, 2, and 14 from filtered blast results in the file "rat2mouse.aln.flt" into a new file "rat2mouse.dat".



2. RUN AGEM
=========

2.1 Combine reciprocal blast results
-------------------------------------
Assuming now you've already had "A2B.dat" and "B2A.dat" (please read 1.1), the next step would be getting a representative homology score for each pair of homologous genes. For example, gene i in genome A and gene j in genome B are homologous genes,  and there are two conditions: 1) a reciprocally homologous relationship --gene i and gene j has homology score in both A2B.dat and B2A.dat; 2) a unilaterally homologous relationship -- gene i and gene j has homology score in either A2B.dat or B2A.dat. AGEM provides the user two options to get teh representative homology score for each pair of homologous genes. 

option 1: combine.py

combine.py gets a representative homology score for both (reciprocal and unilateral ) types of homologous gene pairs. To run this script, type:

python combine.py A2B.dat B2A.dat> AB.dat

option 2: combine_reciprocal.py

combine.py gets a representative homology score for only reciprocally homologous gene pairs, and the unilaterally homologous gene pairs are discarded. To run this script, type:

python combine_reciprocal.py A2B.dat B2A.dat> AB.dat

comment: both commands combine reciprocal blast results in the two files into a new file "AB.dat". In this new file, the homology score is the representative one. 

2.2 run maximum cardinality bipartite matching
----------------------------------------------
The algorithm of Hopcraft-Karp Algorithm is implemented to provide the user a fast way to estimate the maximal number of genes that can be matched betweeen two genomes and generate the matched gene pairs. Assuming you've already had "AB.dat" (read 2.1), to run this, type:

python hopcraft_karp.py AB.dat output_file_name

comment: hopcraft_karp.py generates the output file looks like:

The number of maximal matching is: number

genome.A_gene.i[tab-delimited]genome.B_gene.h(comment: gene.h and gene.i are the one-to-one match)
...                    ...                   ...

2.3 maximum cardinality vs. maximum weighted bipartite matching diagnostics
---------------------------------------------------------------------------

This step takes a much longer time than 2.2. We suggest you to run this on server or any powerful computer other than your personal laptop unless the genomes you want to match are really small.

weighted_diagnostics.py generates the data for you to diagnose the relationships between maximum cardinality and maximum weighted bipartite matching for the genomes you want to match. Specifically, it generates  a file that has two columns (the first column is the cardinality ranging from 0 to maximum, and the second column is the global weight of maximum weighted bipartite matching at a specific cardinality), as well as a plot to visualize the output data. To run this step, type:

python weighted_diagnostics.py AB.dat output_plot_filename >output_data_filename


2.4 maximum weighted bipartite matching at a specific cardinality
-----------------------------------------------------------------

Assuming (after 2.3 step diagnostics or not) you've already had an idea on which cardianlity you want to select as the one for a maximum weighted bipartite matching, now you want to print the matched pairs into an output file. To run, type:

python weighted.py AB.dat k >output_filename

comment: k in the above command is the cardinality or the number of matched pairs you want to specify for a maximum weighted bipartite matching.


2.5  more than one diagnostics curves on the same plot for better comparisons
-----------------------------------------------------------------------------

In some situations, you may want to plot more than one diagnostics curves on the same graph. For example, you want to match genomes in four contrasts such as mouse to rat, mouse to human, mouse to zebrafish etc., and you want to compare these curves. AGEM provides you this function with a maximum of curves being 7. These curves are in seven colors corresponding to input order. 

    input file 1 -black
    input file 2 -blue
    input file 3 -green
    input file 4 -red
    input file 5 -yellow
    input file 6 -cyan
    input file 7 -magenta
    
To run this, you can use weighted_diagnostics.py to generate the diagnostics data and then graph them into the same plot by typing:

python plot_all.py input_1 input_2 ... input_x output_plot_filename


