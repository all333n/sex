#include "sexdiff.h"
#include "MersenneTwister.h"
#include <iostream>
using namespace std;

// Random number generator:

MTRand rnd;

// Pointers on input and output files:

FILE * fichierE;
FILE * fichierS;


int main()
{	
	// Parameters:
	
	int Nt, nbA, NbPrelim, NbGen, paspr, pas;
	double sigi, s, h, sm, hm, saf, haf, sam, ham, U, Ua, L, c, mutM, mutMstep;

	// Opens input and output files:
	bool fin;
	ouvrirFichierE();
	ouvrirFichierS();
	fin = false;
	
	int no = 1;
	do
	{
		//reads parameter values from input file:
		
		fin = lireFichier(Nt, sigi, s, h, sm, hm, saf, haf, sam, ham, U, Ua, nbA, L, c, mutM, mutMstep, NbPrelim, NbGen, paspr, pas); 
					   	  // end = true if end of input file
		if (!fin)		
		{
			// Writes parameter values in output file:
			
			ecrireParametres(Nt, sigi, s, h, sm, hm, saf, haf, sam, ham, U, Ua, nbA, L, c, mutM, mutMstep, NbPrelim, NbGen, paspr, pas); 
			
			// Simulation:
			
			recursion(Nt, sigi, s, h, sm, hm, saf, haf, sam, ham, U, Ua, nbA, L, c, mutM, mutMstep, NbPrelim, NbGen, paspr, pas, no);  
				
			no++;
		}
	} while (!fin);	
	
	// closes files:
	
	fclose(fichierE);
	fclose(fichierS);		
	
	return 0 ;
}
