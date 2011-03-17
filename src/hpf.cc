/**
   \mainpage
   The goal of this project is the computation of RNA-RNA-Interaction
   dynamics using a model of RNA that allows several interaction sites
   and assumes simple hybridization at each site.  In the model, we
   assume that single structures outside of the hybridization sites
   and the hybridization sites itself are each equilibrated.
   
   We start by defining classes for the computation of hybrid
   partition functions and joint probabilities for two unpaired sites.
   
 */



#include  <stdlib.h>
#include  <string.h>
#include  <stdio.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}


// #ifdef _OPENMP
// #include <omp.h>
// #endif


#include <LocARNA/matrices.hh>

#include "unpaired_pf.hh"

#include "hybrid_pf.hh"


int
main()
{
    // set some global variables for Vienna libRNA
    dangles=2;
    temperature = 37.;


    //                      0        1         2         3         4         5         6         7         8         9         0         1         2
    //                      12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    const std::string seqA="GUUGGGAACUAGACCGAUCGCCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUAUUGGGCACCCUGGUGGGGGUUAGUUUAGUUUUUCCCCAGG";

    //                      123456789012345678901234567
    const std::string seqB="ACAUAGCUGAUAAAUCCCGGCUUCUAU";


    UnpairedPF unpaired_pf_A(seqA);
    UnpairedPF unpaired_pf_B(seqB);
    HybridPF hybrid_pf(seqA,seqB);
    
    // print all subsequence hybrid partition functions (for test purposes)
    double kT = GASCONST/1000.0 * (K0+temperature);

    
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t j1=i1; j1<=seqA.length(); j1++) {
	    for (size_t i2=1; i2<=seqB.length(); i2++) {
		for (size_t j2=i2; j2<=seqB.length(); j2++) {
		    HybridPF::pf_t pfhyb = hybrid_pf.get_pf(i1,j1,i2,j2);
		    double ensemble_energy_hyb = - kT * log(pfhyb);
		    
		    UnpairedPF::pf_t upfA=unpaired_pf_A.get_pf(i1,j1);
		    UnpairedPF::pf_t upfB=unpaired_pf_B.get_pf(i2,j2);
		    double ensemble_energy_A = - kT * log(upfA);
		    double ensemble_energy_B = - kT * log(upfB);
		    
		    if (pfhyb >0) {
			std::cout << i1 << " " << j1 << " "
				  << i2 << " " << j2 <<" "
				  << pfhyb << " "
				  << upfA << " " 
			    	  << upfB << " " 
				  << ensemble_energy_hyb << " "
				  << ensemble_energy_A << " "
				  << ensemble_energy_B << " "
				  << std::endl;
		    }
		}
	    }
	}
    }
    
}
