#ifndef SEXDIFF_H
#define SEXDIFF_H

#include <vector>
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include "MersenneTwister.h"
using namespace std;


// Gobal variables:

#define fichierLecture "parametres.txt"     // name of input  
#define fichierEcriture "resultats.txt"		// and output files

// definition of structure "chr" representing a chromosome:
// "mod" is the allele at the modifier locus (here the sex rate coded by the modifier allele),
// "sel" is a vector containing the positions of deleterious alleles along the chromosome
// "ant" is a dynamic_bitset (table of zeros and ones) representing sex antagonistic loci

struct chr
{
	double mod; // modifier locus
	vector<double> sel; // loci under deleterious mutation
	boost::dynamic_bitset<> ant; // loci under antagonistic selection
};



// Prototypes of functions

void ouvrirFichierE();
void ouvrirFichierS();
void ecrireParametres(int Nv, double sigiv, double sv, double hv, double smv, double hmv, double safv, double hafv, double samv, double hamv, double Uv, 
					  double Uav, int nbAv, double Lv, double cv, double mutMv, double mutMstepv, int NbPrelimv, int NbGenv, int pasprv, int pasv);
bool lireFichier(int &Nr, double &sigir, double &sr, double &hr, double &smr, double &hmr, double &safr, double &hafr, double &samr, double &hamr,
					double &Ur, double &Uar, int &nbAr, double &Lr,  double &cr, double &mutMr, double &mutMstepr, int &NbPrelimr, int &NbGenr, 
					int &pasprr, int &pasr);
void recursion(int Nv, double sigiv, double sv, double hv, double smv, double hmv, double safv, double hafv, double samv, double hamv, double Uv, 
			   double Uav, int nbAv, double Lv, double cv, double mutMv, double mutMstepv, int NbPrelimv, int NbGenv, int pasprv, int pasv, int nov);
double gammln(const double xx);
double poisdev(const double xm);
double binldev(const double pp, const int n);
double fitness(chr &c1, chr &c2, double H, double S);
void rec(chr &res, chr &c1, chr &c2, double sz, int B);
void cntl_c_handler(int bidon);

#endif
