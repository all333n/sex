#include "sexdiff.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <cmath>
#include <csignal>
using namespace std;

extern MTRand rnd;
extern FILE * fichierE;
extern FILE * fichierS;

// for stopping the program with Ctrl-C:

bool cntl_c_bool = false;

void cntl_c_handler (int bidon) 
{
	cntl_c_bool = true;
}


/* function recursion: iterates the life cycle. The rate of sex is kept constant 
 during the first NbPrelimv generations, and is free to evolve during the next
 NbGenv generations.

 Other parameters are:

 Nv: population size
 sv, hv: selection and dominance coefficients of deleterious alleles in females
 smv, hmv: selection and dominance coefficients of deleterious alleles in males
 safv, hafv: selection and dominance coefficients of sex-antagonistic alleles in females
 samv, hamv: selection and dominance coefficients of sex-antagonistic alleles in males
 Uv: deleterious mutation rate per haploid genome
 Uav: mutation rate at sexually antagonistic loci
 nbAv: number of sexually antagonistic loci
 Lv: genome map length (average number of cross-overs at meiosis)
 cv: cost of sex
 mutMv: mutation rate at sex modifier locus
 mutMstepv: step size for modifier mutation
 pasv: number of generations between measures of sex rate, mean fitness and diversities
 sampleSv: sample size (per deme) for estimations of heterosis
 */


void recursion(int Nv, double sigiv, double sv, double hv, double smv, double hmv, double safv, double hafv, double samv, 			double hamv, double Uv, double Uav, int nbAv, double Lv, double cv, double mutMv, double mutMstepv, 
		int NbPrelimv, int NbGenv, int pasprv, int pasv, int nov)
{
	// variables:

	int i, j, k, gen, mut, chr1, chr2, ind, nb, nb4, prem1, prem2, par1, par2, nbFix, nHo, nHe, sex, Nft, Nmt;
	double cmpt1, w, wm, wbarm, wbarf, lef, lem, rd, mutPos, md, divf, divm;
	double sxf, sxm;
	double wmaxf, wmaxm;
	bool stp = 0;

	// time interval to check for fixed mutations:

	int pas_c = 20;

	// fitness effects of deleterious and sexually antagonistic mutations in females and males:
	
	double Whet = 1 - hv * sv;
	double Whom = 1 - sv;
	double Whetm = 1 - hmv * smv;
	double Whomm = 1 - smv;
	
	double Whet_a = 1 + hafv * safv;
	double Whom_a = 1 + safv;
	double Whetm_a = 1 + hamv * samv;
	double Whomm_a = 1 + samv;

	// some useful quantities:
	
	double twoL = 2 * Lv;
	double invc = 1 - (1/cv);
	int N_1 = Nv - 1;
	int twoN = 2*Nv;
	int nA_1 = nbAv - 1;

	// effective rate of sex during preliminary generations:

	double sx = (2.0 * sigiv / cv) / (1 - sigiv + (2.0 * sigiv / cv));
	
	// female (pop_f) and male (pop_m) populations: tables of 2N chromosomes:

	chr * pop_f = new chr [twoN];
	chr * temp_f = new chr [twoN];
	chr * pop_m = new chr [twoN];
	chr * temp_m = new chr [twoN];

	// tables of fitness values for females and males:

	double * Wf = new double [Nv];
	double * Wm = new double [Nv];

	// table holding allele frequencies at sexually antagonistic loci:

	double * freqs = new double [nbAv];
	
	// generates the different output files:

	// file "result...txt" keeps track of the rate of sex in females and males, mean fitnes in females and males,
	// mean number of deleterious mutations in females and males, diversities at sex-antagonistic loci in
	// females and males, the number of fixed deleterious alleles and the number of males in the population
	// (one line every "paspr" generations during the preliminary generations, and then one line every 
	// "pasv" generations):

	char nomFichier[256];
	stringstream nomF;
	nomF << "result_N" << Nv << "_sig" << sigiv  << "_sf" << sv << "_hf" << hv << "_sm" << smv << "_hm" << hmv 
		 << "_saf" << safv << "_haf" << hafv << "_sam" << samv << "_ham" << hamv
		 << "_U" << Uv  << "_Ua" << Uav << "_nbA" << nbAv << "_L" << Lv << "_c" << cv << "_" << nov << ".txt";
	nomF >> nomFichier; 
	ofstream fout(nomFichier);

	// file "females...txt" holds the distribution of sigma values among females (the 0 - 1 interval of possible
	// values of sigma is subdivided into 100 classes (between 0 and 0.01, between 0.01 and 0.02...), 
	// the frequency of each class is measured every "pasv" generations):
	
	char nomFichierF[256];
	stringstream nomFm;
	nomFm << "females_N" << Nv << "_sig" << sigiv  << "_sf" << sv << "_hf" << hv << "_sm" << smv << "_hm" << hmv 
		  << "_saf" << safv << "_haf" << hafv << "_sam" << samv << "_ham" << hamv
		  << "_U" << Uv  << "_Ua" << Uav << "_nbA" << nbAv << "_L" << Lv << "_c" << cv << "_" << nov << ".txt";
	nomFm >> nomFichierF; 
	ofstream foutF(nomFichierF);

	// file "males...txt": distribution of sigma values among males:
	
	char nomFichierM[256];
	stringstream nomM;
	nomM << "males_N" << Nv << "_sig" << sigiv  << "_sf" << sv << "_hf" << hv << "_sm" << smv << "_hm" << hmv 
		 << "_saf" << safv << "_haf" << hafv << "_sam" << samv << "_ham" << hamv
		 << "_U" << Uv  << "_Ua" << Uav << "_nbA" << nbAv << "_L" << Lv << "_c" << cv << "_" << nov << ".txt";
	nomM >> nomFichierM; 
	ofstream foutM(nomFichierM);

	// tables for classes of sigma values:
	
	const int nbpts = 100;
	double * sexF = new double[nbpts];
	double * sexM = new double[nbpts];

	// keeping track of starting time:
	
	time_t debut, fin;
	struct tm *ptr;
	debut = time(0);

	// for stopping with Ctrl-C:
	
	cntl_c_bool = false;
	signal(SIGINT, cntl_c_handler);
	
	//initialization: 

	//initial numbers of males and females (N/2):
	
	int Nm = Nv / 2;
	int Nf = Nv / 2;
	int twoNm = 2*Nm;
	int twoNf = 2*Nf;
	int twoNm_1 = twoNm - 1;
	int twoNf_1 = twoNf - 1;

	// the value of each modifier allele is set to "sigiv", no deleterious mutation,
	// all sexually antagonistic loci carry allele 0:
	
	for (i = 0; i < twoNf; i++)
	{
		pop_f[i].mod = sigiv;
		pop_f[i].ant.resize(nbAv);
	}
	
	for (i = 0; i < twoNm; i++)
	{
		pop_m[i].mod = sigiv;
		pop_m[i].ant.resize(nbAv);
	}
	
	nbFix = 0;
	
	//preliminary generations: without mutation at the modifier locus
	
	for (gen = 0; gen < NbPrelimv; gen++)
	{ 
		// the following removes fixed deleterious mutations every "pas_c" generations 
		// (nbFix keeps track of number of fixed mutations):
		
		if (gen % pas_c == 0)
			for (i = pop_f[0].sel.size()-1; i >= 0 ; i--)
			{
				mutPos = pop_f[0].sel[i];
				j = 1;
				while ((j < twoNf) && 
					   (find(pop_f[j].sel.begin(), pop_f[j].sel.end(), mutPos) != pop_f[j].sel.end())) 
					j++;
				if (j == twoNf)
				{
					j = 0;
					while ((j < twoNm) && 
						   (find(pop_m[j].sel.begin(), pop_m[j].sel.end(), mutPos) != pop_m[j].sel.end())) 
						j++;
					if (j == twoNm)
					{
						for (k = 0; k < twoNf; k++)
							pop_f[k].sel.erase(find(pop_f[k].sel.begin(), pop_f[k].sel.end(), mutPos));
						for (k = 0; k < twoNm; k++)
							pop_m[k].sel.erase(find(pop_m[k].sel.begin(), pop_m[k].sel.end(), mutPos));
						nbFix++;
					}
				}
			}
		
		//mutation in females:
		
		for (i = 0; i < twoNf; i++) // for each chromosome in females
		{
			mut = poisdev(Uv); // number of new deleterious mutations on chromosome
			if (mut > 0)
			{
				// each mutation has a random position between 0 and 2L:

				for (j = 0; j < mut; j++)
					pop_f[i].sel.push_back(twoL * rnd.rand());
				sort(pop_f[i].sel.begin(), pop_f[i].sel.end());
			}
			
			mut = poisdev(Uav); // number of mutations at sexually antagonistic loci
			for (j = 0; j < mut; j++)
				pop_f[i].ant.flip(rnd.randInt(nA_1)); // the locus switches 
		}

		//mutation in males:
		
		for (i = 0; i < twoNm; i++)
		{
			mut = poisdev(Uv);
			if (mut > 0)
			{
				for (j = 0; j < mut; j++)
					pop_m[i].sel.push_back(twoL * rnd.rand());
				sort(pop_m[i].sel.begin(), pop_m[i].sel.end());
			}
			
			mut = poisdev(Uav);
			for (j = 0; j < mut; j++)
				pop_m[i].ant.flip(rnd.randInt(nA_1));
		}
		
		// fills table Wf, Wm with fitnesses of all females and males,
		// computes maximal fitness (wmax) and average fitness (wbar) in females and males:
		
		wbarf = 0;
		wbarm = 0;
		wmaxf = 0;
		wmaxm = 0;

		// female fitnesses:

		for (j = 0; j < Nf; j++) // for each female
		{
			nb = 2*j;

			// effect of deleterious mutations:
			w = fitness(pop_f[nb], pop_f[nb+1], Whet, Whom);

			// effect of sexually-antagonistic loci:
			nHo = ((pop_f[nb].ant) & (pop_f[nb+1].ant)).count();
			nHe = ((pop_f[nb].ant) ^ (pop_f[nb+1].ant)).count();
			w *= pow(Whom_a, nHo) * pow(Whet_a, nHe);

			wbarf += w;
			Wf[j] = w;
			if (wmaxf < w)
				wmaxf = w;
		}
		wbarf /= Nf; // mean female fitness

		// male fitnesses:

		for (j = 0; j < Nm; j++) // for each male
		{
			nb = 2*j;
			w = fitness(pop_m[nb], pop_m[nb+1], Whetm, Whomm);
			nHo = ((pop_m[nb].ant) & (pop_m[nb+1].ant)).count();
			nHe = ((pop_m[nb].ant) ^ (pop_m[nb+1].ant)).count();
			w *= pow(Whomm_a, nHo) * pow(Whetm_a, nHe);
			Wm[j] = w;
			wbarm += w;
			if (wmaxm < w)
				wmaxm = w;
		}
		wbarm /= Nm; // mean male fitness
		
		Nft = 0; Nmt = 0;
		
		// sampling the next generation:
		
		for (ind = 0; ind < Nv; ind++) // for each individual of the next generation:
		{
			// sampling the mother:
			
			do{
				chr1 = rnd.randInt(twoNf_1);
				prem1 = chr1/2;
				
			} while (rnd.rand() > Wf[prem1] / wmaxf);
			
			rd = rnd.rand();

			// if the mother reproduces sexually (possible only if there are males):
			
			if ((rd < sx) && (Nm > 0))
			{
				if (rnd.rand() < 0.5)
					sex = 0;	// female offspring
				else
					sex = 1;	// male offspring
				
				// produces recombinant chromosome:
				
				if (sex == 0) // female offspring: stores rec. chromosome in temp_f
				{
					if (chr1 % 2 == 0)
						rec(temp_f[2*Nft], pop_f[chr1], pop_f[chr1 + 1], Lv, nbAv);
					else
						rec(temp_f[2*Nft], pop_f[chr1], pop_f[chr1 - 1], Lv, nbAv);
				}
				else // male offspring: stores rec. chromosome in temp_m
				{
					if (chr1 % 2 == 0)
						rec(temp_m[2*Nmt], pop_f[chr1], pop_f[chr1 + 1], Lv, nbAv);
					else
						rec(temp_m[2*Nmt], pop_f[chr1], pop_f[chr1 - 1], Lv, nbAv);
				}
				
				// sampling the father:
				
				do{
					chr1 = rnd.randInt(twoNm_1);
					prem1 = chr1/2;
					
				} while (rnd.rand() > Wm[prem1] / wmaxm);
				
				// produces recombinant chromosome from the father:
				
				if (sex == 0) // female offspring
				{
					if (chr1 % 2 == 0)
						rec(temp_f[2*Nft + 1], pop_m[chr1], pop_m[chr1 + 1], Lv, nbAv);
					else
						rec(temp_f[2*Nft + 1], pop_m[chr1], pop_m[chr1 - 1], Lv, nbAv);
					Nft++;
				}
				else  // male offspring
				{
					if (chr1 % 2 == 0)
						rec(temp_m[2*Nmt + 1], pop_m[chr1], pop_m[chr1 + 1], Lv, nbAv);
					else
						rec(temp_m[2*Nmt + 1], pop_m[chr1], pop_m[chr1 - 1], Lv, nbAv);
					Nmt++;
				}
			}

			// if the mother reproduces asexually:
	
			else
			{
				par1 = 2*prem1;
				temp_f[2*Nft] = pop_f[par1];
				temp_f[2*Nft + 1] = pop_f[par1 + 1];
				Nft++;
			}	
		} 

		// numbers of females and males in the new generation:
		
		Nf = Nft; twoNf = 2*Nf; twoNf_1 = twoNf - 1;
		Nm = Nmt; twoNm = 2*Nm; twoNm_1 = twoNm - 1;

		// new population of females and males:
		
		for (i = 0; i < twoNf; i++)
			pop_f[i] = temp_f[i];
		for (i = 0; i < twoNm; i++)
			pop_m[i] = temp_m[i];
			
		// stats every "pasprv" generations:
		
		if (gen % pasprv == 0)
		{
			for (i = 0; i < nbAv; i++)
				freqs[i] = 0;
			
			// number of mutations per chromosome in females, and allele frequencies 
			// at sexually antagonistic loci:
			
			lef = 0;
			for (i = 0; i < twoNf; i++)
			{
				lef += pop_f[i].sel.size();
				
				for (j = 0; j < nbAv; j++)
					freqs[j] += pop_f[i].ant[j];
			}
			lef /= twoNf;
			
			// genetic diversity at sexually antagonistic loci among females:
			
			divf = 0;
			for (j = 0; j < nbAv; j++)
				divf += (freqs[j]/twoNf) * (1 - (freqs[j]/twoNf));
			
			for (i = 0; i < nbAv; i++)
				freqs[i] = 0;
			
			// number of mutations per chromosome in males, and allele frequencies 
			// at sexually antagonistic loci:
			
			lem = 0;
			for (i = 0; i < twoNm; i++)
			{
				lem += pop_m[i].sel.size();
				
				for (j = 0; j < nbAv; j++)
					freqs[j] += pop_m[i].ant[j];
			}
			
			// genetic diversity at sexually antagonistic loci among females:
			
			divm = 0;
			if (Nm > 0)
			{
				lem /= twoNm;
				for (j = 0; j < nbAv; j++)
					divm += (freqs[j]/twoNm) * (1 - (freqs[j]/twoNm));
			}
			
			// writes in output file:
				 
			fout << sigiv << " " << sigiv << " " << wbarf << " " << wbarm << " " << lef << " " << lem 
				 << " " << divf << " " << divm << " " << nbFix << " " << Nm << endl;
		}
		
		if (cntl_c_bool)
		{
			stp = 1;
			break;
		}
	}
	
	// letting the rate of sex evolve:
	
	if (!cntl_c_bool)
	{
		for (gen = 0; gen < NbGenv; gen++)
		{ 
			// the following removes fixed deleterious mutations every "pas_c" generations 
			// (nbFix keeps track of number of fixed mutations):
			
			if (gen % pas_c == 0)
				for (i = pop_f[0].sel.size()-1; i >= 0 ; i--)
				{
					mutPos = pop_f[0].sel[i];
					j = 1;
					while ((j < twoNf) && 
						   (find(pop_f[j].sel.begin(), pop_f[j].sel.end(), mutPos) != pop_f[j].sel.end())) 
						j++;
					if (j == twoNf)
					{
						j = 0;
						while ((j < twoNm) && 
							   (find(pop_m[j].sel.begin(), pop_m[j].sel.end(), mutPos) != pop_m[j].sel.end())) 
							j++;
						if (j == twoNm)
						{
							for (k = 0; k < twoNf; k++)
								pop_f[k].sel.erase(find(pop_f[k].sel.begin(), pop_f[k].sel.end(), mutPos));
							for (k = 0; k < twoNm; k++)
								pop_m[k].sel.erase(find(pop_m[k].sel.begin(), pop_m[k].sel.end(), mutPos));
							nbFix++;
						}
					}
				}
			
			//mutation in females:
			
			for (i = 0; i < twoNf; i++)
			{
				//modifier locus:
				if (rnd.rand() < mutMv)
				{
					rd = rnd.rand();
					if (rd < 0.5)
						pop_f[i].mod = rnd.rand(); // with prob. 0.5 the mutant allele codes for any rate of sex
					else
					{
						//else it is close to the value of the parent allele:
						pop_f[i].mod = pop_f[i].mod + rnd.rand() * 2 *mutMstepv - mutMstepv;
						if (pop_f[i].mod < 0)
							pop_f[i].mod = 0;
						else if (pop_f[i].mod > 1)
							pop_f[i].mod = 1;	
					}
				}
				
				//deleterious mutations:
				mut = poisdev(Uv);
				if (mut > 0)
				{
					for (j = 0; j < mut; j++)
						pop_f[i].sel.push_back(twoL * rnd.rand());
					sort(pop_f[i].sel.begin(), pop_f[i].sel.end());
				}
				
				//sexually antagonistic mutations:
				mut = poisdev(Uav);
				for (j = 0; j < mut; j++)
					pop_f[i].ant.flip(rnd.randInt(nA_1));
			}
			
			//mutation in males:
			for (i = 0; i < twoNm; i++)
			{
				//modifier locus:
				if (rnd.rand() < mutMv)
				{
					rd = rnd.rand();
					if (rd < 0.5)
						pop_m[i].mod = rnd.rand();
					else
					{
						pop_m[i].mod = pop_m[i].mod + rnd.rand() * 2 *mutMstepv - mutMstepv;
						if (pop_m[i].mod < 0)
							pop_m[i].mod = 0;
						else if (pop_m[i].mod > 1)
							pop_m[i].mod = 1;	
					}
				}
				
				//deleterious mutations:
				mut = poisdev(Uv);
				if (mut > 0)
				{
					for (j = 0; j < mut; j++)
						pop_m[i].sel.push_back(twoL * rnd.rand());
					sort(pop_m[i].sel.begin(), pop_m[i].sel.end());
				}
				
				//sexually antagonistic mutations:
				mut = poisdev(Uav);
				for (j = 0; j < mut; j++)
					pop_m[i].ant.flip(rnd.randInt(nA_1));
			}
			
			// fills table Wf, Wm with fitnesses of all females and males,
			// computes maximal fitness (wmax) and average fitness (wbar) in females and males:
			
			wbarf = 0;
			wbarm = 0;
			wmaxf = 0;
			wmaxm = 0;
			
			// female fitnesses:
			for (j = 0; j < Nf; j++)
			{
				nb = 2*j;
				
				// effect of deleterious mutations:
				w = fitness(pop_f[nb], pop_f[nb+1], Whet, Whom);
				
				// effect of sexually antagonistic mutations:
				nHo = ((pop_f[nb].ant) & (pop_f[nb+1].ant)).count();
				nHe = ((pop_f[nb].ant) ^ (pop_f[nb+1].ant)).count();
				w *= pow(Whom_a, nHo) * pow(Whet_a, nHe);
				wbarf += w;
				
				// effect of the cost of sex (paid only if there are males, otherwise all females reproduce asexually):
				md = (pop_f[nb].mod + pop_f[nb+1].mod) / 2.0;
				if (Nm > 0)
					w *= (1 - md + (2.0 * md / cv));
				else
					w *= (1 - md);
				Wf[j] = w;
				if (wmaxf < w)
					wmaxf = w;
			}
			wbarf /= Nf;
			
			// male fitnesses:
			for (j = 0; j < Nm; j++)
			{
				nb = 2*j;
				
				// effect of deleterious mutations:
				w = fitness(pop_m[nb], pop_m[nb+1], Whetm, Whomm);
				
				// effect of sexually antagonistic mutations:
				nHo = ((pop_m[nb].ant) & (pop_m[nb+1].ant)).count();
				nHe = ((pop_m[nb].ant) ^ (pop_m[nb+1].ant)).count();
				w *= pow(Whomm_a, nHo) * pow(Whetm_a, nHe);
				
				Wm[j] = w;
				wbarm += w;
				if (wmaxm < w)
					wmaxm = w;
			}
			wbarm /= Nm;
			
			Nft = 0; Nmt = 0;
			
			// sampling the next generation:
			
			for (ind = 0; ind < Nv; ind++)
			{
				// sampling the mother:
				
				do{
					chr1 = rnd.randInt(twoNf_1);
					prem1 = chr1/2;
					
				} while (rnd.rand() > Wf[prem1] / wmaxf);
				
				// rate of sex of the mother:
				
				par1 = 2*prem1;
				md = (pop_f[par1].mod + pop_f[par1+1].mod) / 2.0;
				if (Nm > 0)
					sx = (2.0 * md / cv) / (1 - md + (2.0 * md / cv));
				else
					sx = 0;
				rd = rnd.rand();
				
				// if sexually produced offspring:
				
				if ((rd < sx) && (Nm > 0))
				{
					if (rnd.rand() < 0.5)
						sex = 0;	// female offspring
					else
						sex = 1;	// male offspring
					
					// recombination:
					
					if (sex == 0)
					{
						if (chr1 % 2 == 0)
							rec(temp_f[2*Nft], pop_f[chr1], pop_f[chr1 + 1], Lv, nbAv);
						else
							rec(temp_f[2*Nft], pop_f[chr1], pop_f[chr1 - 1], Lv, nbAv);
					}
					else
					{
						if (chr1 % 2 == 0)
							rec(temp_m[2*Nmt], pop_f[chr1], pop_f[chr1 + 1], Lv, nbAv);
						else
							rec(temp_m[2*Nmt], pop_f[chr1], pop_f[chr1 - 1], Lv, nbAv);
					}
					
					// sampling a father:
					
					do{
						chr1 = rnd.randInt(twoNm_1);
						prem1 = chr1/2;
						
					} while (rnd.rand() > Wm[prem1] / wmaxm);
					
					// recombination:
					
					if (sex == 0)
					{
						if (chr1 % 2 == 0)
							rec(temp_f[2*Nft + 1], pop_m[chr1], pop_m[chr1 + 1], Lv, nbAv);
						else
							rec(temp_f[2*Nft + 1], pop_m[chr1], pop_m[chr1 - 1], Lv, nbAv);
						Nft++;
					}
					else
					{
						if (chr1 % 2 == 0)
							rec(temp_m[2*Nmt + 1], pop_m[chr1], pop_m[chr1 + 1], Lv, nbAv);
						else
							rec(temp_m[2*Nmt + 1], pop_m[chr1], pop_m[chr1 - 1], Lv, nbAv);
						Nmt++;
					}
				}
				
				// if asexually produced offspring:
				
				else
				{
					temp_f[2*Nft] = pop_f[par1];
					temp_f[2*Nft + 1] = pop_f[par1 + 1];
					Nft++;
				}	
			} 
			
			// numbers of females and males in the new generation:
			
			Nf = Nft; twoNf = 2*Nf; twoNf_1 = twoNf - 1;
			Nm = Nmt; twoNm = 2*Nm; twoNm_1 = twoNm - 1;
			
			// new population of females and males:
			
			for (i = 0; i < twoNf; i++)
				pop_f[i] = temp_f[i];
			for (i = 0; i < twoNm; i++)
				pop_m[i] = temp_m[i];
			
			// stats every "pasv" generations:
			
			if (gen % pasv == 0)
			{
				for (i = 0; i < nbpts; i++)
				{	
					sexF[i] = 0;
					sexM[i] = 0;
				}
				
				sxf = 0;
				sxm = 0;
				for (i = 0; i < nbAv; i++)
					freqs[i] = 0;
				
				// number of mutations per chromosome in females, sex rates distribution
				// and allele frequencies at sexually antagonistic loci:
				
				lef = 0; 
				for (i = 0; i < twoNf; i++)
				{
					lef += pop_f[i].sel.size();
					sxf += pop_f[i].mod;
					
					for (j = 0; j < nbAv; j++)
						freqs[j] += pop_f[i].ant[j];
					
					if (pop_f[i].mod < 1)
						sexF[int(floor(nbpts * pop_f[i].mod))] += 1;
					else
						sexF[nbpts - 1] += 1;
				}
				sxf /= twoNf;
				lef /= twoNf;
				
				// genetic diversity at sexually antagonistic loci in females:
				
				divf = 0;
				for (j = 0; j < nbAv; j++)
					divf += (freqs[j]/twoNf) * (1 - (freqs[j]/twoNf));
				
				for (i = 0; i < nbAv; i++)
					freqs[i] = 0;
				
				// number of mutations per chromosome in females, sex rates distribution
				// and allele frequencies at sexually antagonistic loci:
				
				lem = 0;
				for (i = 0; i < twoNm; i++)
				{
					lem += pop_m[i].sel.size();
					sxm += pop_m[i].mod;
					
					for (j = 0; j < nbAv; j++)
						freqs[j] += pop_m[i].ant[j];
					
					if (pop_m[i].mod < 1)
						sexM[int(floor(nbpts * pop_m[i].mod))] += 1;
					else
						sexM[nbpts - 1] += 1;
				}
				
				// genetic diversity at sexually antagonistic loci in females:
				
				divm = 0;
				if (twoNm > 0)
				{
					sxm /= twoNm;
					lem /= twoNm;
					for (j = 0; j < nbAv; j++)
						divm += (freqs[j]/twoNm) * (1 - (freqs[j]/twoNm));
				}
				
				// writes in output file:
				
				fout << sxf << " " << sxm << " " << wbarf << " " << wbarm << " " << lef << " " << lem 
					 << " " << divf << " " << divm << " " << nbFix << " " << Nm << endl;
				
				// writes distributions of sigma values in females and males:
				
				for (i = 0; i < nbpts; i++)
				{	
					foutF << sexF[i] / twoNf << " ";
					foutM << (twoNm > 0 ? (sexM[i] / twoNm) : 0) << " ";
				}
				foutF << endl; foutM << endl;
			}
			
			if (cntl_c_bool)
				break;
		}
	}
	
	fin = time(0);  
	
	if (cntl_c_bool)
	{
		if (stp == 1)
			fprintf(fichierS, "\n\nInterrupted by user at preliminary generation %d ", gen);
		else
			fprintf(fichierS, "\n\nInterrupted by user at generation %d ", gen);
	}
	
	// Writes results:
	fprintf(fichierS, "\n\nResults in file ");
	fprintf(fichierS, nomFichier);
	fprintf(fichierS, "\n");
	         
	// Simulation time length:
	int temps = difftime(fin, debut);
	fprintf(fichierS,
		 "\n%d generations took %d hour(s) %d minute(s) %d seconds\n",
		 NbGenv, temps / 3600, (temps % 3600) / 60, temps % 60);
		
	// Date and time:
	ptr=localtime(&fin);
	fprintf(fichierS, asctime(ptr));
	
	delete [] pop_f;
	delete [] temp_f;
	delete [] pop_m;
	delete [] temp_m;
	delete [] Wf;
	delete [] Wm;
	delete [] sexF;
	delete [] sexM;
}
