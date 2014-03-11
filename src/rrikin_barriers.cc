/**
 * @todo Try the use of a lru cache
 * (e.g. http://code.google.com/p/lru-cache-cpp/). Check first if the
 * hash causes our space problems with large instances.
 *
 * @todo what happens if the input list is incomplete? I.e. there are
 * lower energy neighbors which are not in the hash? Check how exactly
 * the tests are performed.
 *
 * @todo change the transitions data structure to array of hashs
 */

/**
 * @file rrikin_barriers.cc
 *
 * Defines main() of the program rrikin_barriers
 *
 * rrikin_barriers constructs the barrier tree/graph and rate matrix for the
 * macro state process that moves between basins in the energy
 * landscape.
 *
 * @note Due to our non-standard definition of transition energies, it
 * becomes necessary to rethink the definition of local minima and and
 * basins in the energy landscape!  Usually these are defined in terms
 * of (non-transition) state energies, comparing the energies of
 * source and target state in each move/transition.  In general, these
 * notions should be defined in terms of transition state energies,
 * since these energies govern the speed of transitions.  We define a
 * neighborship relation x ->_t y with the semantics that there is a
 * move from x to y with transition state t.  For a state (set of
 * states/ensemble) s, E_s denotes it's energy (ensemble energy) and
 * Z_s denotes it's Boltzmann weight (partition function).  A state x
 * is a local minimum, iff forall x ->_t y: E_x <= E_t. Note that the
 * common definition puts E_y in place of E_t. Analogously, the
 * criterion for adaptive walks and steepest descent walk (i.e. basin
 * assignment) is modified to compare to the transition state energy
 * in place of the target state energy.
 *
 * IS THIS TRUE?: Note that both (i.e., old and new) local minimum
 * definitions are equivalent for Kawasaki and Metropolis rates. The
 * steepest descent definition is equivalent for the case of Kawasaki
 * rates but not *equivalent* for Metropolis. This could indicate a
 * problem of the new (or the old) definition!?: For x->_t y, where
 * E_t<=E_x, which implies E_y<=E_x, the comparison to E_t does not
 * differentiate in the case of Metropolis rates, since there E_t=E_x!
 *
 * @note In our setting, for all x->_t y the property holds that E_t
 * >= min(E_x,E_y). [Note that this is a consequence of either the
 * source or the target being a relaxation of the transition state and
 * Prop.: For energies of ensembles x and y, there holds that x
 * subset y implies E_x >= E_y.] This important property allows us to
 * use a barrier-like algorithm, where states are processed in the
 * order of increasing energy, for identifying local minima and
 * assigning basin membership to all states.
 *
 * @note For the computation of macro-state rates for the transitions
 * between basins it suffices to accumulate the partition functions of
 * basins and transition states between basins. Note that the
 * definition e.g. given by Flamm, Hofacker, 2007 is equivalent to the
 * Arrhenius rate derived from these partition functions, i.e.
 * k(a->b) = sum_x in a,y in b Pr[x|a]k(x->y) = Z^trans_ab / Z_a,
 * where Z^trans_ab = sum_x in a,y in b: Z^trans_xy.
 *
 * @note Although the model is degenerate, steepest descent walks are 
 * uniquely defined, since the order on moves is fixed
 *
 * @note Barriers/transition states between basins: the energy of the
 * barrier has to be the maximum of source state energy, transition
 * state energy, and target state energy
 * (cf. transition_t::barrier_energy()). Otherwise, the barrier can be
 * lower than one of the explicit states. In this case the barrier
 * would not reflect that the energy has to be raised to this explicit
 * state for the transition between the basins.
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

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}

/* control output */
bool debug_out;
bool verbose;

/* control behavior */
bool simplify_graph;

bool consider_double_sites;

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_barriers_cmdline.h"

/* Methods of basin */

void Basin::print_header(std::ostream &out) const {
    printf("%5s %-32s %10s %6s %6s",
	   "idx",
	   "description",
	   "n_s",
	   "ensE",
	   "minE"
	   );
}

void
Basin::print(std::ostream &out, const HybEnsModel &model) const {
    printf("%5lu %-32s %10.2f %6.2f %6.2f",
	   index_,
	   local_minimum.toString().c_str(),
	   states,
	   - model.RT() * log(Z),
	   minimum_energy
	   );
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

    
    // names of molecules for rxns and spcs output
    std::string nameA = "A";
    std::string nameB = "B";
    
    if (args_info.nameA_given) {
	nameA=args_info.nameA_arg;
    }
    
    bool nameB_given=args_info.nameB_given;
    if (args_info.nameB_given) {
	nameB=args_info.nameB_arg;
    }
    
    bool nameAB_given=args_info.nameAB_given;
    std::string nameAB;
    if (args_info.nameAB_given) {
	nameAB=args_info.nameAB_arg;
    }
    
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

    double max_outflow  = args_info.max_outflow_arg;
    double min_rate     = args_info.min_rate_arg;
    double min_p_equ    = args_info.min_p_equ_arg;
    bool binary         = args_info.binary_given;
    bool special_open_state = ! args_info.no_special_open_state_given;
    bool gradient       = ! args_info.no_gradient_given;
    
    
    consider_double_sites = ! args_info.no_double_sites_given;

    simplify_graph      = ! args_info.dont_simplify_graph_given;
    verbose             = args_info.verbose_given;
    debug_out           = args_info.debug_given;

    size_t to_keep_given       = args_info.to_keep_given;
    std::set<size_t> to_keep_set;
    for(size_t i=0; i<to_keep_given;++i) {
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


    std::string seqA = args_info.inputs[0];
    HybEnsModel::norm_RNA_seq(seqA);
    std::string seqB = "";
    
    if (homodimer) {
	seqB = seqA;
	HybEnsModel::reverse(seqB);
    }  else if (antisense) {
	seqB = seqA;
	HybEnsModel::complement(seqB);
    } else {
	seqB = args_info.inputs[1];
	HybEnsModel::norm_RNA_seq(seqB);
    }


    if (verbose) {
	std::cerr << "seqA="<<seqA<<", seqB="<<seqB << std::endl;
    }

    cmdline_parser_free(&args_info);
    
    // global settings for Vienna libRNA
    dangles=2;
    
    if (verbose) {
	std::cerr << "Construct barrier graph." << std::endl;
    }

    stopwatch.start("construct");

    // construct barrier graph
    BarrierGraph bg(seqA,seqB,binary,
		    special_open_state,
		    consider_double_sites,
		    gradient,
		    verbose,
		    debug_out);
    
    stopwatch.stop("construct");
    
    size_t num_total_basins = bg.num_basins();
    

    if (verbose) {
	if (bg.model().is_homodimer()) {
	    std::cerr << "Model homodimer."<<std::endl;
	}
	std::cerr << "Generated "<<num_total_basins<<" basins." << std::endl;
	bg.print_stats(std::cerr);
	stopwatch.print_info(std::cerr);
    }

    
    if (simplify_graph) {

	if (verbose) {
	    std::cerr << "Merge basins with outflow larger " << max_outflow << " or equilibrium probability smaller "<< min_p_equ << std::endl;
	    std::cerr << "Remove rates smaller than " << min_rate << std::endl;	    
	}
	
	stopwatch.start("merge");
	
	bg.merge_basins_by_outflow(max_outflow,min_p_equ,min_rate);
	
	stopwatch.stop("merge");
	
	if (verbose) {
	    std::cerr << "Merged "<<num_total_basins-bg.num_basins()
		      <<" basins resulting in "<<bg.num_basins()<<" states." << std::endl;
	    bg.print_stats(std::cerr);
	}
	
    }

    if (verbose) {
	std::cerr << "Reindex" << std::endl;
    }
    bg.reindex();

    if (special_open_state) {
	// Exchange basin indices 1 and 2 such that the global minimum is state 1
	// and the open state is state 2
    
	if (verbose) {
	    std::cerr << "Move global minimum to smallest index; open state to second index." << std::endl;
	}
	bg.swap_indices(0,1);
    }
    
    if (verbose) {
	bg.print_stats(std::cerr);
    }
    
    /* handle use-defined basin merging */
    
    if (to_keep_given) {

	if (verbose) {
	    std::cerr << "Keep only the "<<to_keep_given<<" specified basins." << std::endl;
	}
	
	bg.reduce_basin_set(to_keep_set,min_rate);

	if (verbose) {
	    std::cerr << "Reindex" << std::endl;
	}
	bg.reindex();
	
    }


    // print basins of barrier graph
    bg.print_basins(std::cout);    

    std::cout << std::endl
	      << std::endl;
    bg.print_barrier_graph(std::cout);
    
    std::vector<size_t> components;
    std::vector<size_t> component_sizes=bg.connected_components(components);
    if (component_sizes.size()>1) {
    	if (verbose) {
	    std::cerr << "Components: #="<<component_sizes.size()<<" sizes: ";
	    for (size_t i=0; i<component_sizes.size(); i++)
		std::cerr << component_sizes[i]<<" ";
	    std::cerr <<std::endl;
	
	    std::cerr << "Keep only the first component (which contains the open state.)" 
		      << std::endl;
	}
	
	bg.keep_single_component(1,components);
	
	if (verbose) {
	    bg.print_stats(std::cerr);
	}
    }
    
	    
    if (verbose) {
	double max_diff = bg.check_rates();
	if (max_diff > 1e-12) {
	    std::cerr << "WARNING: maximal deviation from detailed balance:"<<max_diff<<std::endl;
	}
    }

    if (barfile != "") {
	if (verbose) {
	    std::cerr << "Write bar file for treekin '"<<barfile<<"'."<<std::endl;	
	}
	
	std::ofstream fout(barfile.c_str());
	if (fout.good()) {
	    bg.print_barriers(fout);
	    fout.close();
	} else {
	    std::cerr << "Cannot write barriers file."<<std::endl;
	}
    }
    
        
    if (ratesfile != "") {
	if (verbose) {
	    std::cerr << "Write rates matrix to file '"<<ratesfile<<"'."<<std::endl;
	}
	
	std::ofstream fout(ratesfile.c_str());
	if (fout.good()) {
	    bg.print_treekin_ratesmatrix(fout);
	} else {
	    std::cerr << "Cannot write rates file."<<std::endl;
	}
	fout.close();
	
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
    if (!nameB_given && bg.model().is_homodimer()) {
	nameB=nameA;
    }    
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


        
    if (verbose) {
	stopwatch.print_info(std::cerr);
    }

    exit(0);
}
