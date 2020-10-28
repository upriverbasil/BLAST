TEAM MEMBERS - Bhavay Aggarwal, Chandan Gupta, Gavish Gupta, Saad Ahmad
COURSE - Algorithms In BioInformatics
INTRODUCTION

This project contains an implementation of BLAST(Basic Local Alignment Search Tool). BLAST helps to align the query sequence with the sequence database and gives an optimal result in reasonable time. The steps used in the implementation of BLAST are:
	a) Creating K-mers of length = word_length from the database sequence.

                              ||
                              ||
                              \/

	b) Defining a suitable hash function to faster the search of sequences.

                              ||
                              ||
                              \/

	c) Creating a Binary Search Tree. Each node of BST contains a hashcode (for a particular K-mer) and its position in the database sequence.
                              ||
                              ||
                              \/

	d) Creating K-mers of length = word_length from the query sequence.

                              ||
                              ||
                              \/

	e) For each K-mer from the query sequence find all the possible mutated strings within the HSSP threshold. Form a list for each K-mer separately and assign a hash code to every member of the list using the same hash function.
                              ||
                              ||
                              \/

	f) Search the BST for exact hashcode matches and if found a match store the K-mer starting index in query sequence and corresponding matched K-mer's position(s) in the database.

                              ||
                              ||
                              \/

	g) Seed and Extend all the hits found between the K-mers using Smith–Waterman algorithm and traceback to find optimal alignments.
                              ||
                              ||
                              \/

	h) Generate Statistics(bit-scores, pvalue, evalue) using Karlin-Altschul formula
                              ||
                              ||
                              \/

	i) Collects form parameters from web page and processes on backend built on flutter framework. Displays result on a new page.

CODE DESCRIPTION - 

CLASS Node - Used to store nodes of Binary Tree used in efficient searching of queries K-mers.

_init__(...) 				        -> Initialises and defines all the necessary variables.
__lt__(...)  			  	      -> Used to sort the class objects on basis of hash score
__repr__(...)				        -> Used to print details of object

CLASS K_Mer - Used to store K_mer sequences with starting positions in both database and query sequence.

_init__(...) 				        -> Initialises and defines all the necessary variables.
__lt__(...)  				        -> Used to sort the class objects on basis of hash score
__repr__(...)				        -> Used to print details of object

CLASS BLAST - 

__init__(...) 				      -> Initialises and defines all the necessary variables.
kmers_positions(...)		    -> Used to form a list of kmers with their start positions in the database
kmers_hash_table(...) 		  -> Used to form list of hashcodes of kmers obtained through kmers_positions(...)
hash_code(...)				      -> Calculates and returns hash code of given sequence
make_binary_tree(...) 		  -> Used sorted array instead of Binary Tree
match_kmer_binary_tree(...)	-> Uses binary search to find kmers exists or not that are obtained through 								   recur(...) 
recur(...)                  -> Find all the possible combinations of kmer within HSSP threshold
binarysearch(...)           -> Simple binary search module
smith_waterman(...)			    -> Seeds and Extends alignment and find score and stats of the alignment.
statistics(...)				      -> Defining the statistics required(Karlin-Altschul) 

USAGE - 1) Clone the repository. 
        2) source env/bin/activate
        3) python app.py
        4) Upload a database .fasta file. (Upload sequence.fasta present in extras folder of this directory)
        5) Enter the query sequence (eg. GCCTATACAGTTGAACTCGGTACAGAAGTAAATGAGTTCGCCTGTGTTGTGGCAGATGCTGTCATAAAAACTTTGCAACCAGTATCTGAATTACTTACACCACTGGGCATTGATTTAGATGAGTGGAGTATGGCTACATACTACTTATTTGATGAGTCTGGTGAGTTTAAATTGGCTTC).
        6) Enter Other Hyperparameters.
        7) Wait for the results on new screen.

EXAMPLES - 

