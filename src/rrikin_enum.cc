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

#include <iostream>
#include <LocARNA/stopwatch.hh>

/* control output */
bool verbose;

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_enum_cmdline.h"
#include "rri_enumeration.hh"
#include "hybrid_ensemble_model.hh"

int
main(int argc, char **argv)
{
    // speed improvement when writing to standard streams: desync c and c++ i/o
    std::ios_base::sync_with_stdio(false);

    LocARNA::StopWatch stopwatch(false);
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
    
    const double max_hyb_energy   = args_info.max_hyb_energy_arg;
    const double max_total_energy = args_info.max_total_energy_arg;
    const size_t maxsitesize = 
	(args_info.max_hyb_length_arg>=0)
	? args_info.max_hyb_length_arg
	: std::max(seqA.length(),seqB.length());
    
    const size_t maxsitesize_diff = 
	(args_info.max_hyb_length_diff_arg>=0)
	? args_info.max_hyb_length_diff_arg
	: std::max(seqA.length(),seqB.length());

    size_t region_startA=1;
    size_t region_endA=seqA.length();
    size_t region_startB=1;
    size_t region_endB=seqB.length();
    
    parse_region(args_info.regionA_arg,region_startA,region_endA);
    parse_region(args_info.regionB_arg,region_startB,region_endB);
    
    const size_t span = 
	args_info.span_arg>=0
	? args_info.span_arg
	: std::numeric_limits<size_t>::max();
    
    const size_t window = 
	args_info.span_arg>=0
	? args_info.span_arg*2
	: std::numeric_limits<size_t>::max();

    verbose        = args_info.verbose_given;
    bool enum_double_sites        = ! args_info.no_double_sites_given;
    bool binary        = args_info.binary_given;
    bool add_open_state        = args_info.add_open_state_given;
    
    // ------------------------------------------------------------
    // enumerate states
    stopwatch.start("init_model");

    if (verbose) std::cerr << "Initialize model (precomputing energies for sequences of length "
			   <<seqA.size()<<" and "<<seqB.size()
			   << "; region A " << region_startA << "-" << region_endA
			   << "; region B " << region_startB << "-" << region_endB
			   <<")" << std::endl;
    
    HybEnsModel model(seqA,seqB,
		      maxsitesize,
		      maxsitesize_diff,
		      region_startA,
		      region_endA,
		      region_startB,
		      region_endB,
		      span,
		      window,
		      enum_double_sites);

    stopwatch.stop("init_model");

    RRIEnumeration enumeration(model,
			       max_hyb_energy,
			       max_total_energy,
			       binary);
    
    //if (verbose) stopwatch.print_info(std::cerr);

    stopwatch.start("enumerate");

    if (verbose) std::cerr << "Enumerate states ( max hyb length "<<maxsitesize
			   << ", max hyb length diff "<<maxsitesize_diff
			   << ", max ss hyb energy " <<max_hyb_energy
			   << ", max total energy " << max_total_energy 
			   << " )" << std::endl;
    
    // 0 interaction sites
    
    if (add_open_state) {	
	HybEnsModel::StateDescription empty_state;
	if (verbose) {
	    std::cerr <<"Write open state"<<std::endl;
	}
	if (model.energy(empty_state) <= max_total_energy) {
	    enumeration.write_state(std::cout,model.energy(empty_state),
				    empty_state,
				    binary);
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

    //stopwatch.start("enum_single");
    size_t count_single_states =
	enumeration.enumerate_single_sites(false);
    //stopwatch.stop("enum_single");

    if (verbose) {
	std::cerr << "Enumerated "<<count_single_states<<" single site states"<<std::endl;
	//stopwatch.print_info(std::cerr);
    }

    if (enum_double_sites) {
	stopwatch.start("enum_double");
		
	size_t count_double_states=0;
	if (verbose) {
	    std::cerr << "Enumerate double hybridization site states"
		      << std::endl;
	}
	
	count_double_states = 
	    enumeration.enumerate_double_sites(false);
	
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
