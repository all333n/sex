The program simulates N diploid individuals, each individual carries 2 chromosomes along which deleterious mutations occur (at rate U per chromosome per generation), and a number nbA of loci under sexually antagonistic selection. The sex modifier is located at the mid-point of the chromosome. When an individual reproduces, the number of cross-overs at meiosis is sampled from Poisson distribution of parameter L, and cross-overs occur evenly along the chromosome.

There are first a given number of preliminary generations during which the rate of sex is fixed, then the rate of sex can evolve. The program produces a file "resultats.txt" which stores parameter values and the time length of the simulation, and another file (result_N..txt) where each line gives the mean rate of sex in females and males, mean fitness in females and males, mean number of deleterious mutations per chromosome (in females and males), genetic diversity at sexually antagonistic loci (in females and males), number of fixed mutations and number of males in the population, measured at different generations.
The program also generates two files "femalesÉtxt" and "maleÉtxt" giving the distribution of sigma values in females and males at different timepoints.

The file "sexdiff.h" provides function prototypes and global variables, in particular the structure "chr" that represents a chromosome. The main function (in "main.cpp") calls functions (defined in "files.cpp") that read parameters in an input file ("parametres.txt"), write parameters in the file "resultats.txt" and then calls the function "recursion" that does the simulation. The file "parametres.txt" provides parameter values in the order indicated at the beginning of the file; each line corresponding to a parameter set must begin with * (one can put several lines corresponding to several parameter sets, the program will do one after the other). The file "SelRec.cpp" contains functions that calculate the fitness of individuals and construct recombinant chromosomes. 

Note that the program uses the dynamic-bitset class, defined in the boost library (www.boost.org)

Using gcc, the program can by compiled using the command
g++ -O3 *.cpp

