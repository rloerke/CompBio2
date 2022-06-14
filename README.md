# Computational Biology Assignment 2
## Written by Ray Loerke

Implementation of a primitive gene-finder and predictor program.
This program will scan a bacterial geneome for genes and predict their products.
This program accepts a FASTA (.fasta) file containing the sequence to be scanned.
Three output FASTA files are generated containing the gene sequence (.g.fasta), RNA sequences transcribed from the genes (.r.fasta) and the protein sequences translated from the genes (.p.fasts).

See the example files for the formatting of these file types. 
Generally, a leading '<' means a sequence follows and a leading '#' means a comment follows.
