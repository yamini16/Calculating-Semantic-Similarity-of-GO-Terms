Final Project for 520 Fall 2015

Submitted to : Dr. Rasiah Loganantharaj


Gene ontology (GO) describes the attributes of genes and gene products (either RNA or protein, resulting from expression of a gene) using a structured and controlled vocabulary. GO consists of three ontologies: biological process (BP), cellular component (CC) and molecular function (MF), each of which is modeled as a directed acyclic graph[1].

In this project, we have implemented a method to measure the semantic similarity of GO terms for human genes using Python. For calculating the semantic similarity between two GO terms, we have used the formula from[1].

To execute the file:
run script geneSim.py with the 2 GENE_IDs as the command line parameters

On the first run, a 'cache' will be created to used for subsequent similarity queries.

Submitted by:

Prakash Subedi

Matinsadat Hosseini

Yamini Joshi  


References:
[1] Measure the Semantic Similarity of GO Terms Using Aggregate Information Content.
Xuebo Song, Lin Li, Pradip K. Srimani, Philip S. Yu, and James Z. Wang
