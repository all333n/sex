// Functions to open input and output files, 
// read parameter values and write them in output file

#include "sexdiff.h"
#include <iostream>
#include <fstream>
using namespace std;

extern FILE * fichierE;
extern FILE * fichierS;


//Opens input file:

void ouvrirFichierE()    
{						 
	fichierE = fopen(fichierLecture,"r");
}


//Opens output file:

void ouvrirFichierS()   
{
	fichierS = fopen(fichierEcriture,"a");
}


// Reads parameter values from input file,
// returns 1 if it reaches the end of input file, else returns 0:

bool lireFichier(int &Nr, double &sigir, double &sr, double &hr, double &smr, double &hmr, double &safr, double &hafr, double &samr, double &hamr,
				 double &Ur, double &Uar, int &nbAr, double &Lr,  double &cr, double &mutMr, double &mutMstepr, int &NbPrelimr, int &NbGenr, 
				 int &pasprr, int &pasr)  
{
	int x;
	bool term;
	do {x = fgetc(fichierE);} while (!((x == '*') || (x == EOF)));
		// Lines with parameter sets must begin with *
	if (x == EOF)
	{
		term = true;
	}
	else
	{
		fscanf(fichierE,"%d ",&Nr);
		fscanf(fichierE,"%lf ",&sigir);
		fscanf(fichierE,"%lf ",&sr);
		fscanf(fichierE,"%lf ",&hr);
		fscanf(fichierE,"%lf ",&smr);
		fscanf(fichierE,"%lf ",&hmr);
		fscanf(fichierE,"%lf ",&safr);
		fscanf(fichierE,"%lf ",&hafr);
		fscanf(fichierE,"%lf ",&samr);
		fscanf(fichierE,"%lf ",&hamr);
		fscanf(fichierE,"%lf ",&Ur);
		fscanf(fichierE,"%lf ",&Uar);
		fscanf(fichierE,"%d ",&nbAr);
		fscanf(fichierE,"%lf ",&Lr);
		fscanf(fichierE,"%lf ",&cr);
		fscanf(fichierE,"%lf ",&mutMr);
		fscanf(fichierE,"%lf ",&mutMstepr);
		fscanf(fichierE,"%d ",&NbPrelimr);
		fscanf(fichierE,"%d ",&NbGenr);
		fscanf(fichierE,"%d ",&pasprr);
		fscanf(fichierE,"%d ",&pasr);
		
		term = false;
	} 
	return term;
}


//Writes parameter values in output file:

void ecrireParametres(int Nv, double sigiv, double sv, double hv, double smv, double hmv, double safv, double hafv, double samv, double hamv, double Uv, 
					  double Uav, int nbAv, double Lv, double cv, double mutMv, double mutMstepv, int NbPrelimv, int NbGenv, int pasprv, int pasv)   
							
{
	fprintf(fichierS,"\n_________________________________________\n");
	fprintf(fichierS,"\nN = %d", Nv);
	fprintf(fichierS,", sig = %g", sigiv);
	fprintf(fichierS,", s = %g", sv);
	fprintf(fichierS,", h = %g", hv);
	fprintf(fichierS,", sm = %g", smv);
	fprintf(fichierS,", hm = %g", hmv);
	fprintf(fichierS,", saf = %g", safv);
	fprintf(fichierS,", haf = %g", hafv);
	fprintf(fichierS,", sam = %g", samv);
	fprintf(fichierS,", ham = %g", hamv);
	fprintf(fichierS,", U = %g", Uv);
	fprintf(fichierS,", Ua = %g", Uav);
	fprintf(fichierS,", nbA = %d", nbAv);
	fprintf(fichierS,", L = %g", Lv);
	fprintf(fichierS,", c = %g", cv);
	fprintf(fichierS,"\nmutM = %g", mutMv);
	fprintf(fichierS,", mutMstep = %g", mutMstepv);
	fprintf(fichierS,", prelim = %d", NbPrelimv);
	fprintf(fichierS,", generations = %d", NbGenv);
	fprintf(fichierS,", paspr = %d", pasprv);
	fprintf(fichierS,", pas = %d", pasv);
}
