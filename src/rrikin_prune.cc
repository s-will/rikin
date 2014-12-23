/**
 * @file rrikin_prune.cc
 *
 * Defines main() of the program rrikin_prune
 */

#include "basin_transition.hh"
#include "barrier_graph.hh"

#include <LocARNA/stopwatch.hh>


#include  <sstream>

#include  <stdlib.h>
#include  <string>
#include  <math.h>
#include <assert.h>
#include  <fstream>
#include  <iomanip>
#include  <stack>

#include <zlib.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}



bool consider_double_sites;

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_prune_cmdline.h"

int
main(int argc, char **argv)
{
    LocARNA::StopWatch stopwatch;
    stopwatch.start("total");

    gengetopt_args_info args_info;
    
    // get options (call gengetopt command line parser)
    if (cmdline_parser (argc, argv, &args_info) != 0)
	exit(1) ;

    if ( args_info.inputs_num != 1 ) {
	std::cerr << "Expect name of input file on command line."<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }

    std::string inputfile = args_info.inputs[0];
    
    // names of molecules for rxns and spcs output
    std::string nameA = "A";
    std::string nameB = "B";
    
    if (args_info.nameA_given) {
	nameA=args_info.nameA_arg;
    }
    
    if (args_info.nameB_given) {
	nameB=args_info.nameB_arg;
    }
    
    bool nameAB_given=args_info.nameAB_given;
    std::string nameAB;
    if (args_info.nameAB_given) {
	nameAB=args_info.nameAB_arg;
    }

    double max_outflow  = args_info.max_outflow_arg;
    double min_pequ    = args_info.min_pequ_arg;
    
    double qu_outflow = 
	args_info.qu_outflow_given
	? args_info.qu_outflow_arg
	: std::numeric_limits<double>::max();

    double qu_pequ = 
	args_info.qu_pequ_given
	? args_info.qu_pequ_arg
	: std::numeric_limits<double>::max();

    size_t num_outflow = 
	args_info.num_outflow_given
	? args_info.num_outflow_arg
	: std::numeric_limits<size_t>::max();

    size_t num_pequ = 
	args_info.num_pequ_given
	? args_info.num_pequ_arg
	: std::numeric_limits<size_t>::max();  

    
    double min_rate = args_info.min_rate_arg;
    
    bool special_first_state = ! args_info.no_special_first_state_given;

    /* control behavior */
    bool simplify_graph      = ! args_info.dont_simplify_graph_given;
    
    /* control output */
    bool verbose             = args_info.verbose_given;
    bool debug_out           = args_info.debug_given;

    /* tracking */
    bool track = args_info.track_given;
    std::string track_file;
    if (track) { track_file = args_info.track_arg; }


    std::set<size_t> to_keep_set;
    for(size_t i=0; i<args_info.to_keep_given;++i) {
	to_keep_set.insert(args_info.to_keep_arg[i]-1);
    }

    std::string ratesfile;
    if (args_info.ratesfile_given) {
	ratesfile = args_info.ratesfile_arg;
    }

    std::string barfile;
    if (args_info.barfile_given) {
    	barfile = args_info.barfile_arg;
    }

    std::string pffile;
    if (args_info.pffile_given) {
	pffile = args_info.pffile_arg;
    }

    std::string rxnsfile;
    if (args_info.rxns_given) {
        rxnsfile = args_info.rxns_arg;
    }

    std::string spcsfile;
    if (args_info.spcs_given) {
        spcsfile = args_info.spcs_arg;
    }
    
    cmdline_parser_free(&args_info);
    
    //stopwatch.start("read");

    // construct barrier graph
    std::ifstream in(inputfile.c_str(), std::ios::in | std::ios::binary);
    BarrierGraph bg(in,
		    min_rate,
		    special_first_state,
		    verbose,
		    debug_out,
		    track);
    in.close();
    
    //stopwatch.stop("read");
    
    
    size_t num_total_basins = bg.num_basins();
    
    if (verbose) {
	std::cerr << "Read in "<<num_total_basins<<" states." << std::endl;
	bg.print_stats(std::cerr);
	//stopwatch.print_info(std::cerr);
    }

    if (simplify_graph) {
	
	// if requested, determine max-outflow and min-pequ by quantile or basin number
	if ( 0<=qu_outflow && qu_outflow<=100 ) {
	    max_outflow = std::min( max_outflow,
				    bg.max_outflow_by_quantile(qu_outflow) );
	}
	if ( num_outflow <= num_total_basins ) {
	    max_outflow = std::min( max_outflow,
				    bg.max_outflow_by_number(num_outflow) );
	}
	if ( 0<=qu_pequ && qu_pequ<=100 ) {
	    min_pequ = std::max( min_pequ, 
				 bg.min_pequ_by_quantile(qu_pequ) );
	}
	if ( num_pequ <= num_total_basins ) {
	    min_pequ = std::max( min_pequ, 
				 bg.min_pequ_by_number(num_pequ) );
	}
	
	if (verbose) {
	    std::cerr << "Prune basins with outflow larger " << max_outflow 
		      << " or equilibrium probability smaller "<< min_pequ << std::endl;
	    std::cerr << "Remove rates smaller than " << min_rate << std::endl;	    
	}
	
	stopwatch.start("prune");
	
	bg.prune(max_outflow,min_pequ,min_rate);
	
	stopwatch.stop("prune");
	
	if (verbose) {
	    std::cerr << "Dissolved "<<num_total_basins-bg.num_basins()
		      <<" basins resulting in "<<bg.num_basins()<<" states." << std::endl;
	    bg.print_stats(std::cerr);
	}
	
    }

    if (verbose) {
	std::cerr << "Reindex" << std::endl;
    }
    bg.reindex();

    if (special_first_state) {
	// Exchange basin indices 1 and 2 (1-based) such that the global minimum is state 1
	// and the first state is state 2
    
	// if (verbose) {
	//     std::cerr << "Move global minimum to smallest index; first state to second index." << std::endl;
	// }
	bg.swap_indices(0,1);
    }
    
    // handle user-defined basin merging
    
    if (!to_keep_set.empty()) {
	if (verbose) {
	    std::cerr << "Keep only the "<<to_keep_set.size()<<" specified basins." << std::endl;
	}
	
	bg.reduce_basin_set(to_keep_set,min_rate);
	
	if (verbose) {
	    std::cerr << "Reindex" << std::endl;
	}
	bg.reindex();

	if (verbose) {
	    bg.print_stats(std::cerr);
	}
    }
    
    std::vector<size_t> components;
    std::vector<size_t> component_sizes=bg.connected_components(components);
    if (component_sizes.size()>1) {
    	if (verbose) {
	    std::cerr << "Components: #="<<component_sizes.size()<<" sizes: ";
	    for (size_t i=0; i<component_sizes.size(); i++)
		std::cerr << component_sizes[i]<<" ";
	    std::cerr <<std::endl;
	
	    std::cerr << "Keep only the first component." 
		      << std::endl;
	}
	
	bg.keep_single_component(1,components); // note: this reindexes internally
	
	if (verbose) {
	    bg.print_stats(std::cerr);
	}
    }
    


    // print basins of barrier graph
    bg.print_basins(std::cout);    

    std::cout << std::endl
	      << std::endl;
    bg.print_barrier_graph(std::cout);


	    
    double max_diff = bg.check_rates();
    if (max_diff > 1e-12) {
	std::cerr << "WARNING: maximal deviation from detailed balance:"<<max_diff<<std::endl;
    }
    
    if (barfile != "") {
    	if (verbose) {
    	    std::cerr << "Write pseudo-bar file for treekin '"<<barfile<<"'."<<std::endl;	
    	}
	
    	std::ofstream fout(barfile.c_str());
    	if (fout.good()) {
    	    bg.print_treekin_barriers(fout, "seqname", 1.0 /* RT */);
    	    fout.close();
    	} else {
    	    std::cerr << "Cannot write barriers file."<<std::endl;
    	}
    }
    
        
    if (ratesfile != "") {
	if (verbose) {
	    std::cerr << "Write rates matrix to file '"<<ratesfile<<"'."<<std::endl;
	}
	if (track_file.substr(ratesfile.length()-3,3)==".gz") {
	    // write gzip'd
	    gzFile fh=gzopen(ratesfile.c_str(),"wb");
	    if (fh==NULL) {
		std::cerr << "Cannot write rates matrix to file "<<ratesfile<<"."<<std::endl;
	    } else {
		bg.gzwrite_treekin_ratesmatrix(fh);
	    }
	    gzclose(fh);
	} else {
	    // write uncompressed
	    std::ofstream fout(ratesfile.c_str());
	    if (fout.good()) {
		bg.write_treekin_ratesmatrix(fout);
	    } else {
		std::cerr << "Cannot write rates file."<<std::endl;
	    }
	    fout.close();
	}
    }

    if (pffile != "") {
    
	if (verbose) {
	    std::cerr << "Write partition functions of basins and transition states to file '"<<pffile<<"'."<<std::endl;
	}
	std::ofstream fout(pffile.c_str(),std::ios::out | std::ios::binary);
	if (fout.good()) {
	    bg.print_pfs(fout,true);
	} else {
	    std::cerr << "Cannot write partition functions to file."<<std::endl;
	}
	fout.close();
	
    }

    //handle names of molecules in rxns/spcs files
    if (!nameAB_given) {
	nameAB=nameA+nameB;
    }

    if (rxnsfile != "") {
    
	if (verbose) {
	    std::cerr << "Write reactions to file '"<<rxnsfile<<"'."<<std::endl;
	}
	std::ofstream fout(rxnsfile.c_str(),std::ios::out | std::ios::binary);
	if (fout.good()) {
	    bg.print_rxns(fout,nameA,nameB,nameAB);
	} else {
	    std::cerr << "Cannot write reactions to file."<<std::endl;
	}
	fout.close();
	
    }

    if (spcsfile != "") {
    
	if (verbose) {
	    std::cerr << "Write reactions to file '"<<spcsfile<<"'."<<std::endl;
	}
	std::ofstream fout(spcsfile.c_str(),std::ios::out | std::ios::binary);
	if (fout.good()) {
	    bg.print_spcs(fout,nameA,nameB,nameAB);
	} else {
	    std::cerr << "Cannot write reactions to file."<<std::endl;
	}
	fout.close();
	
    }


    if (track) {
	bool sparse_prune_track=true;
	stopwatch.start("write_pruning_info");
	if (track_file.substr(track_file.length()-3,3)==".gz") {
	    // write gzip'd
	    gzFile fh=gzopen(track_file.c_str(),"wb");
	    if (fh==NULL) {
		std::cerr << "Cannot write prune track to file "<<track_file<<"."<<std::endl;
	    } else {
		bg.gzwrite_pruning_track(fh,sparse_prune_track);
	    }
	    gzclose(fh);
	} else {
	    // write uncompressed
	    try {
		std::ofstream out(track_file.c_str());
		bg.write_pruning_track(out,sparse_prune_track);
		out.close();
	    } catch(const std::ofstream::failure &e) {
		std::cerr << "Cannot write pruning track to file "<<track_file<<". "<<e.what()<<std::endl;
	    }
	}
	
	stopwatch.stop("write_pruning_info");
    }
    
    
    if (verbose) {
	stopwatch.print_info(std::cerr);
    }

    exit(0);
}
