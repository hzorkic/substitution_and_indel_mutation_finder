# Substitution, Insertion, and Deletion Finder

A python script used to parse through SARS-CoV-2 sequences from the GISAID database. 
AFter HOURS of searching for a pre-made python package that found both SUBSTITUTIONS and INSERTION and DELETIONS, 
I opted to make my own. This script is not compelete! But I will continue to update it as I make changes. :chart_with_upwards_trend:

## Input:
a .fasta file of sequences, and a reference sequence

## Output:
a csv file that lists subs and indels in the following format:

[sub/indel], | [RefAA+Position+SeqAA]

Insertion    | B789B

Deletion     | F45-

Substitution | F9T



