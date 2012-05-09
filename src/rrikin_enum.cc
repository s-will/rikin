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



#include  <iostream>

#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}


using std::cout;
using std::endl;

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_enum_cmdline.h"

#include <LocARNA/matrices.hh>
#include "hybrid_ensemble_model.hh"


//! @brief Check validity of states (for debugging)
void
check_state_validity(const HybEnsModel::StateDescription &state, 
		     const HybEnsModel &model) {
#ifndef NDEBUG
    if (not state.is_valid(model)) {
	std::cerr << "ERROR: generated state "<<state<<" is not valid in model."<<std::endl;
	abort();
    }
#endif
}

int
main(int argc, char **argv)
{
    gengetopt_args_info args_info;
    
    // get options (call gengetopt command line parser)
    if (cmdline_parser (argc, argv, &args_info) != 0)
	exit(1) ;
    
    if ( args_info.inputs_num != 2 ) {
	std::cerr << "Expect two sequences as input on command line"<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }
    
    std::string seqA(args_info.inputs[0]);
    std::string seqB(args_info.inputs[1]);
        
    // set some global variables for Vienna libRNA
    dangles=2;
    
    const double th_hyb_energy = args_info.hyb_energy_threshold_arg;

    // ------------------------------------------------------------
    // enumerate states
    
    HybEnsModel model(seqA,seqB);
    
    const size_t minsitesize=model.minsitesize();
    const size_t minsitedist=model.minsitedist();

    // 0 interaction sites
    
    HybEnsModel::StateDescription empty_state;
    
    check_state_validity(empty_state,model);
    printf("%f\n",model.energy(empty_state));
    
    // Indexing for single hybridization 
    // ----\        /-----
    //     i1------j1     
    //     i2------j2     
    //  ---/        \--------
    
    //cout << "Single Hybridisations:" << endl;
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t i2=1; i2<=seqB.length(); i2++) {
	    if ( (model.pair_type(i1,i2)==0) ) {
		continue;
	    }
	    for (size_t j1=i1+minsitesize-1; j1<=seqA.length(); j1++) {
		for (size_t j2=i2+minsitesize-1; j2<=seqB.length(); j2++) {
		    if ( (model.pair_type(j1,j2))==0 ) {
			continue;
		    }
		    
		    HybEnsModel::StateDescription state(i1,i2,j1,j2);
		    		    
		    check_state_validity(state,model);

		    double energy_hyb = 
			model.energy_hybrid(state[0]);
		    
		    if ( energy_hyb > th_hyb_energy ) continue;
		    
		    // double energy_unpair =
		    // 	model.energy_unpair(state[0]);
		    
		    // double total_energy = 
		    // 	energy_hyb
		    // 	+ energy_unpair;
		    
		    double total_energy=model.energy(state);

		    
		    // using cout<< instead of printf causes has extrem overhead
		    printf("%f %lu %lu %lu %lu\n",total_energy,i1,i2,j1,j2);
		}
	    }
	}
    }

    // Indexing for double hybridization 
    // ----\        /--------\       /---------
    //     i1------j1        k1-----l1
    //     i2------j2        k2-----l2
    //  ---/        \-------/         \------------

    
    //cout << "Double Hybridisations:" << endl;
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t i2=1; i2<=seqB.length(); i2++) {
	    if ( (model.pair_type(i1,i2)==0) ) {
		continue;
	    }
	    
	    for (size_t j1=i1+minsitesize-1; j1<=seqA.length(); j1++) {
		for (size_t j2=i2+minsitesize-1; j2<=seqB.length(); j2++) {
		    
		    if ( (model.pair_type(j1,j2)==0) ) {
			continue;
		    }
		    
		    HybEnsModel::StateDescription::ISite is1(i1,i2,j1,j2);
		    
		    double energy_hyb1 = model.energy_hybrid(is1);
		    
		    if ( energy_hyb1 > th_hyb_energy ) continue;
		    
		    for (size_t k1=j1+minsitedist+1; k1<=seqA.length(); k1++) {
			for (size_t k2=j2+minsitedist+1; k2<=seqB.length(); k2++) {
			    if ( (model.pair_type(k1,k2)==0) ) {
				continue;
			    }
			    
			    for (size_t l1=k1+minsitesize-1; l1<=seqA.length(); l1++) {
				for (size_t l2=k2+minsitesize-1; l2<=seqB.length(); l2++) {
				    if ( (model.pair_type(l1,l2)==0) ) {
					continue;
				    }
				    
				    HybEnsModel::StateDescription::ISite is2(k1,k2,l1,l2);

				    HybEnsModel::StateDescription state(i1,i2,j1,j2,k1,k2,l1,l2);

				    check_state_validity(state,model);

				    double energy_hyb2 = model.energy_hybrid(is2);
				    
				    if (energy_hyb2 > th_hyb_energy) continue;
				    
				    //double energy_unpair=model.energy_unpair(is1,is2);
				    
				    //double total_energy = energy_hyb1 + energy_hyb2 + energy_unpair;
				    
				    double total_energy=model.energy(state);
				    
				    // using cout<< instead of printf causes has extrem overhead
				    printf("%f %lu %lu %lu %lu %lu %lu %lu %lu\n",total_energy,i1,i2,j1,j2,k1,k2,l1,l2);
				    
				}
			    }
			}
		    }
		}
	    }
	}
    }

    cmdline_parser_free(&args_info);
    exit(0);
}
