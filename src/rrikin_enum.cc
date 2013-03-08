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
#include  <stdio.h>

#include  <string.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}

#include <LocARNA/stopwatch.hh>


/* control output */
bool verbose;

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


void
write_state(double energy, const HybEnsModel::StateDescription &state, bool binary) {
    if (binary) {
	printf("%.2f ",energy);
	
	HybEnsModel::StateDescription::code_t code;
	state.encode(code);
	fwrite(reinterpret_cast<char *>(&code),sizeof(char),8,stdout);
	fputc(0,stdout);
    } else {
	if (state.size()==0) {
	    printf("%6.2f\n",energy);	
	} else if (state.size()==1) {
	    printf("%6.2f %3lu %3lu %3lu %3lu\n",
		   energy,
		   state[0].i1,state[0].i2,state[0].j1,state[0].j2
		   );
	} else {
	    printf("%6.2f %3lu %3lu %3lu %3lu %3lu %3lu %3lu %3lu\n",
		   energy,
		   state[0].i1,state[0].i2,state[0].j1,state[0].j2,
		   state[1].i1,state[1].i2,state[1].j1,state[1].j2
		   );
	}
    }
}

size_t
enumerate_double_sites(const HybEnsModel &model,
			    double th_hyb_energy,
		       double th_total_energy,
		       bool binary) {
    
    // Indexing for double hybridization 
    // ----\        /--------\       /---------
    //     i1------j1        k1-----l1
    //     i2------j2        k2-----l2
    //  ---/        \-------/         \------------
    
    size_t count_double_states=0;

    const std::string &seqA = model.seqA();
    const std::string &seqB = model.seqB();
    const size_t minsitesize=model.minsitesize();
    const size_t minsitedist=model.minsitedist();

	
    //cout << "Double Hybridisations:" << endl;
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	
	double progress = (i1/(double)seqA.size());
	if (verbose) std::cerr << "\r" << (int(progress*10000)/100.0) << " %   ";

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
				    if (total_energy <= th_total_energy) {
					write_state(total_energy,state,binary);
					count_double_states++;
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

    return count_double_states;
}

int
main(int argc, char **argv)
{

    LocARNA::StopWatch stopwatch;
    stopwatch.start("total");

    gengetopt_args_info args_info;
    
    // get options (call gengetopt command line parser)
    if (cmdline_parser (argc, argv, &args_info) != 0)
	exit(1) ;

    bool homodimer=args_info.homodimer_given;
    bool antisense=args_info.antisense_given;
    
    size_t expected_sequences=(homodimer||antisense)?1:2;

    if (homodimer && antisense) {
	std::cerr << "Options homodimer and antisense are mutually exclusive."<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }
    
    if ( args_info.inputs_num != expected_sequences ) {
	std::cerr << "Expect "<<expected_sequences<<" sequence(s) on command line."<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }
    
    std::string seqA = args_info.inputs[0];
    HybEnsModel::norm_RNA_seq(seqA);
    std::string seqB = "";
    
    if (homodimer) {
	seqB = seqA;
	HybEnsModel::reverse(seqB);
    } else if (antisense) {
	seqB = seqA;
	HybEnsModel::complement(seqB);
    } else {
	seqB = args_info.inputs[1];
	HybEnsModel::norm_RNA_seq(seqB);
    }
        
    // set some global variables for Vienna libRNA
    dangles=2;
    
    const double th_hyb_energy = args_info.max_hyb_energy_arg;
    const double th_total_energy = args_info.max_total_energy_arg;
    verbose        = args_info.verbose_given;
    bool enum_double_sites        = ! args_info.no_double_sites_given;
    bool binary        = args_info.binary_given;
    bool add_open_state        = args_info.add_open_state_given;
    
    // ------------------------------------------------------------
    // enumerate states
    stopwatch.start("generate_model");

    if (verbose) std::cerr << "Generate Model (precomputing energies for sequences of length "
			   <<seqA.size()<<" and "<<seqB.size()<<")" << std::endl;
    HybEnsModel model(seqA,seqB);

    stopwatch.stop("generate_model");
    
    if (verbose) stopwatch.print_info(std::cerr);

    stopwatch.start("enumerate");

    if (verbose) std::cerr << "Enumerate states (maximal single site hybridization energy " <<th_hyb_energy <<  ", max total energy " << th_total_energy << " )" << std::endl;
    
    const size_t minsitesize=model.minsitesize();

    // 0 interaction sites
    
    if (add_open_state) {	
	HybEnsModel::StateDescription empty_state;
	if (verbose) {
	    std::cerr <<"Write open state"<<std::endl;
	}
	if (model.energy(empty_state) <= th_total_energy) {
	    check_state_validity(empty_state,model);
	    write_state(model.energy(empty_state),empty_state,binary);
	}
    }
    
    // Indexing for single hybridization 
    // ----\        /-----
    //     i1------j1     
    //     i2------j2     
    //  ---/        \--------

    
    if (verbose) {
	std::cerr << "Enumerate single hybridization site states"
		  << std::endl;
    }
    stopwatch.start("enum_single");
    size_t count_single_states=0;
    
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
		    
		    // using cout<< instead of printf causes extrem overhead
		    if (total_energy <= th_total_energy) {
			write_state(total_energy,state,binary);
			count_single_states++;
		    }
		}
	    }
	}
    }
    stopwatch.stop("enum_single");

    if (verbose) {
	std::cerr << "Enumerated "<<count_single_states<<" single site states"<<std::endl;
	stopwatch.print_info(std::cerr);
    }


    if (enum_double_sites) {
	stopwatch.start("enum_double");
		
	size_t count_double_states=0;
	if (verbose) {
	    std::cerr << "Enumerate double hybridization site states"
		      << std::endl;    
	}
	
	count_double_states=enumerate_double_sites(model,th_hyb_energy,th_total_energy,binary);
    

	if (verbose) {
	    if (verbose) std::cerr << "\r";
	    std::cerr <<"Enumerated "<<count_double_states<<" double site states"<<std::endl;
	}
    
	stopwatch.stop("enum_double");
    }
    
    stopwatch.stop("enumerate");

    cmdline_parser_free(&args_info);

    stopwatch.stop("total");

    if (verbose) {
	stopwatch.print_info(std::cerr);
    }    

    exit(0);
}
