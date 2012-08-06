// Functions for fitness computation and for recombination:

#include "sexdiff.h"
#include <vector>
#include <cmath>
using namespace std;

extern MTRand rnd;

//Computes the fitness of a diploid with chromosomes c1 and c2.
//Multiplicative model, wHe = 1-hs, wHo = 1-s:

double fitness(chr &c1, chr &c2, double wHe, double wHo)
{
	int s1 = c1.sel.size() - 1;
	int s2 = c2.sel.size() - 1;
	
	double w = 1.0;
	
	while ((s1 > -1) || (s2 > -1))
	{
		if (s1 == -1)
		{
			w *= pow(wHe, s2 + 1);
			break;
		}
		
		if (s2 == -1)
		{
			w *= pow(wHe, s1 + 1);
			break;
		}
		
		// heterozygous mutation on c1:
		
		if (c1.sel[s1] > c2.sel[s2])
		{
			w *= wHe;
			s1--;
			continue;
		}
		
		// heterozygous mutation on c2:
		
		if (c1.sel[s1] < c2.sel[s2])
		{
			w *= wHe;
			s2--;
			continue;
		}
		
		// homozygous mutation:
		
		w *= wHo;
		s1--;
		s2--;
	}
	
	return w;
}




// Constructs chromosome "res" by recombining chromosomes c1 and c2
// "Sz" is the map length L (average number of cross-overs per meiosis).
// "B" is the number of sexually antagonistic loci.
// The modifier gene is inherited from chromosome c1.

void rec(chr &res, chr &c1, chr &c2, double Sz, int B)
{
	int s1 = c1.sel.size() - 1;
	int s2 = c2.sel.size() - 1;
	
	boost::dynamic_bitset<> rec;
	boost::dynamic_bitset<> off1;
	boost::dynamic_bitset<> off2;
	
	res.sel.clear();
	res.ant.clear();
	res.mod = c1.mod;
	
	int j, cmpt;
	vector<double> Co;
	
	double twoS = 2 * Sz;
	double inter = twoS / (B + 1);
	
	// number of cross-overs:
	
	int nbCo = int(poisdev(Sz));
	
	// positions of cross-overs (between 0 and 2L) are put in the vector "Co":
	
	for (j = 0; j < nbCo; j++)
		Co.push_back(twoS * rnd.rand());
	sort(Co.begin(), Co.end());
	
	// loci under deleterious mutation:
	
	// counting the number of cross-overs on the right of the modifier locus:
	
	cmpt = 0;
	j = Co.size() - 1;
	while ((j > -1) && (Co[j] > Sz))
	{
		cmpt++;
		j--;
	}
	
	// if the number of cross-overs on the right is even:
	
	if (cmpt % 2 == 0)
	{	
		// for each cross-over (starting on extreme right):
		
		for (j = 1; j <= cmpt; j++)
		{
			if (j % 2 == 1)
			{
				// all mutations on the right of cross-over on chromosome 1 are incorporated
				
				while ((s1 > -1) && (c1.sel[s1] > Co[nbCo - j]))
				{
					res.sel.push_back(c1.sel[s1]);
					s1--;
				}
				while ((s2 > -1) && (c2.sel[s2] > Co[nbCo - j]))
					s2--;
			}
			else
			{
				// all mutations on the right of cross-over on chromosome 2 are incorporated
				
				while ((s2 > -1) && (c2.sel[s2] > Co[nbCo - j]))
				{
					res.sel.push_back(c2.sel[s2]);
					s2--;
				}
				while ((s1 > -1) && (c1.sel[s1] > Co[nbCo - j]))
					s1--;
			}
		}
	}
	
	// if the number of cross-overs on the right is odd:
	
	else
	{	
		// for each cross-over (starting on extreme right):
		
		for (j = 1; j <= cmpt; j++)
		{
			if (j % 2 == 1)
			{
				// all mutations on the right of cross-over on chromosome 2 are incorporated
				
				while ((s2 > -1) && (c2.sel[s2] > Co[nbCo - j]))
				{
					res.sel.push_back(c2.sel[s2]);
					s2--;
				}
				while ((s1 > -1) && (c1.sel[s1] > Co[nbCo - j]))
					s1--;
			}
			else
			{
				// all mutations on the right of cross-over on chromosome 1 are incorporated
				
				while ((s1 > -1) && (c1.sel[s1] > Co[nbCo - j]))
				{
					res.sel.push_back(c1.sel[s1]);
					s1--;
				}
				while ((s2 > -1) && (c2.sel[s2] > Co[nbCo - j]))
					s2--;
			}
		}
	}
	
	// mutations between the modifier and the nearest cross-over on the right (on chromosome 1):
	
	while ((s1 > -1) && (c1.sel[s1] > Sz))
	{
		res.sel.push_back(c1.sel[s1]);
		s1--;
	}
	while ((s2 > -1) && (c2.sel[s2] > Sz))
		s2--;
	
	// number of cross-overs on the left of the modifier locus:
	
	int frst = nbCo - cmpt;
	
	// for each cross-over (starting with the nearest to the modifier on the left):
	
	for (j = 1; j <= frst; j++)
	{
		if (j % 2 == 1)
		{
			// all mutations on the right of cross-over on chromosome 1 are incorporated
			
			while ((s1 > -1) && (c1.sel[s1] > Co[frst - j]))
			{
				res.sel.push_back(c1.sel[s1]);
				s1--;
			}
			while ((s2 > -1) && (c2.sel[s2] > Co[frst - j]))
				s2--;
		}
		else
		{
			// all mutations on the right of cross-over on chromosome 2 are incorporated
			
			while ((s2 > -1) && (c2.sel[s2] > Co[frst - j]))
			{
				res.sel.push_back(c2.sel[s2]);
				s2--;
			}
			while ((s1 > -1) && (c1.sel[s1] > Co[frst - j]))
				s1--;
		}
	}
	
	// mutations on the left of the left-most cross-over:
	
	if (frst % 2 == 0)
		while (s1 > -1)
		{
			res.sel.push_back(c1.sel[s1]);
			s1--;
		}
	else
		while (s2 > -1)
		{
			res.sel.push_back(c2.sel[s2]);
			s2--;
		}
	
	// sorts mutations on offspring chromosome:
	
	sort(res.sel.begin(), res.sel.end());
	
	// recombination mask (for loci under antagonistic selection):
	
	for (j = 0; j < nbCo; j++)
		rec.resize(int(floor(Co[j] / inter)), (j % 2) == 0 ? 0 : 1);
	rec.resize(B, (nbCo % 2) == 0 ? 0 : 1);
	
	// loci under antagonistic selection:
	
	if ((frst % 2) == 0)
	{
		off2 = c2.ant & rec;
		rec.flip();
		off1 = c1.ant & rec;
	}
	else
	{
		off1 = (c1.ant & rec);
		rec.flip();
		off2 = (c2.ant & rec);
	}
	res.ant = (off1 | off2);
}
