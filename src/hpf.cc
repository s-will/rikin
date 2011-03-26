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

#include "hybrid_ensemble.hh"

int
main()
{
    // set some global variables for Vienna libRNA
    dangles=2;
    temperature = 37.;


    //                      0        1         2         3         4         5         6         7         8         9         0         1         2
    //                      12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    // const std::string seqA="GUUGGGAACUAGACCGAUCGCCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUAUUGGGCACCCUGGUGGGGGUUAGUUUAGUUUUUCCCCAGG";
    const std::string seqA="CCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUA";
    //                                        UAUCUU...           ...UACA
    //const std::string seqA="CCCCGGGG";
    
    //                      123456789012345678901234567
    const std::string seqB="ACAUAGCUGAUAAAUCCCGGCUUCUAU";
    
    //const std::string seqB="GGGGGCCCCC";
    
    

    UnpairedPF unpaired_pf_A(seqA);
    UnpairedPF unpaired_pf_B(seqB);
    HybridPF hybrid_pf(seqA,seqB);
    
    // print all subsequence hybrid partition functions (for test purposes)
    double kT = GASCONST/1000.0 * (K0+temperature);


    size_t min_hyb_size = 2; // minimal size of hybridisation
    
    // threshold on single site hybridization ensemble energy
    double th_hyb_energy=1;

    // Indexing for single hybridization 
    // ----\        /-----
    //     i1------j1     
    //     i2------j2     
    //  ---/        \--------


    std::cout << "Single Hybridisations:" << std::endl;
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t j1=i1+min_hyb_size-1; j1<=seqA.length(); j1++) {
	    for (size_t i2=1; i2<=seqB.length(); i2++) {
		for (size_t j2=i2+min_hyb_size-1; j2<=seqB.length(); j2++) {
		    
		    HybridPF::pf_t pfhyb = hybrid_pf.get_pf(i1,j1,i2,j2);
		    double ensemble_energy_hyb = - kT * log(pfhyb);

		    if ( ensemble_energy_hyb > th_hyb_energy ) continue;
		    
		    UnpairedPF::pf_t upfA=unpaired_pf_A.get_unpaired_prob_single(i1,j1);
		    UnpairedPF::pf_t upfB=unpaired_pf_B.get_unpaired_prob_single(i2,j2);
		    double ensemble_energy_A = - kT * log(upfA);
		    double ensemble_energy_B = - kT * log(upfB);
		    
		    double total_ensemble_energy = 
			ensemble_energy_hyb
			+ ensemble_energy_A
			+ ensemble_energy_B;
		    
		    std::cout << i1 << " " << j1 << " "
			      << i2 << " " << j2 <<" "
			      << total_ensemble_energy << " ("
			      << ensemble_energy_hyb << "+"
			      << ensemble_energy_A << "+"
			      << ensemble_energy_B << ")"
			      << std::endl;
		}
	    }
	}
    }

    
    // Indexing for double hybridization 
    // ----\        /--------\       /---------
    //     i1------j1        k1-----l1
    //     i2------j2        k2-----l2
    //  ---/        \-------/         \------------

    
    std::cout << "Double Hybridisations:" << std::endl;
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t j1=i1+min_hyb_size-1; j1<=seqA.length(); j1++) {
	    
	    for (size_t i2=1; i2<=seqB.length(); i2++) {
		for (size_t j2=i2+min_hyb_size-1; j2<=seqB.length(); j2++) {
		    
		    HybridPF::pf_t pfhyb_ij = hybrid_pf.get_pf(i1,j1,i2,j2);
		    double ensemble_energy_hyb_ij = - kT * log(pfhyb_ij);
		    
		    if ( ensemble_energy_hyb_ij > th_hyb_energy ) continue;
			    
		    for (size_t k1=j1+2; k1<=seqA.length(); k1++) {
			for (size_t l1=k1+min_hyb_size-1; l1<=seqA.length(); l1++) {
			    
			    UnpairedPF::pf_t upfA=
				unpaired_pf_A.get_unpaired_prob_single(k1,l1)
				* unpaired_pf_A.get_unpaired_prob_conditional(i1,j1,k1,l1);
			    
			    double ensemble_energy_A = - kT * log(upfA);
			    
			    for (size_t k2=j2+2; k2<=seqB.length(); k2++) {
				for (size_t l2=k2+min_hyb_size-1; l2<=seqB.length(); l2++) {
				    
				    UnpairedPF::pf_t upfB=
					unpaired_pf_B.get_unpaired_prob_single(k2,l2)
					* unpaired_pf_B.get_unpaired_prob_conditional(i2,j2,k2,l2);
				    
				    HybridPF::pf_t pfhyb_kl = hybrid_pf.get_pf(k1,l1,k2,l2);
				    double ensemble_energy_hyb_kl = - kT * log(pfhyb_kl);
				    
				    if (ensemble_energy_hyb_kl > th_hyb_energy) continue;
				    
				    double ensemble_energy_B = - kT * log(upfB);
				    
				    double total_ensemble_energy = 
					ensemble_energy_hyb_ij
					+ ensemble_energy_hyb_kl
					+ ensemble_energy_A
					+ ensemble_energy_B;
				    
				    
				    std::cout << i1 << " " << j1 << " "
					      << i2 << " " << j2 <<" "
					      << k1 << " " << l1 << " "
					      << k2 << " " << l2 <<" "
						  << total_ensemble_energy << " ("
					      << ensemble_energy_hyb_ij << "+"
					      << ensemble_energy_hyb_kl << "+"
					      << ensemble_energy_A << "+"
					      << ensemble_energy_B << ")"
					      << std::endl;
				    
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
