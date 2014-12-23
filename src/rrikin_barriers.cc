/**
 * @todo what happens if the input list is incomplete? I.e. there are
 * lower energy neighbors which are not in the hash? Check how exactly
 * the tests are performed.
 *
 * @todo change the transitions data structure to array of hashs
 *
 * @todo (low priority) Try the use of a lru cache
 * (e.g. http://code.google.com/p/lru-cache-cpp/). Check first if the
 * hash causes space problems with large instances.
 */

/**
 * @file rrikin_barriers.cc
 *
 * Defines main() of the program rrikin_barriers
 *
 * rrikin_barriers constructs the barrier tree/graph and rate matrix for the
 * macro state process that moves between gradient basins in the energy
 * landscape.
 *
 * @todo check notes
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
#include "rri_barrier_graph.hh"

#include <LocARNA/stopwatch.hh>


#include <cstdlib>
#include <cmath>
#include <cassert>

#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>

#include  <zlib.h>

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


int
main(int argc, char **argv)
{
    LocARNA::StopWatch stopwatch;
    stopwatch.start("total");

    gengetopt_args_info args_info;
    
    // get options (call gengetopt command line parser)
    if (cmdline_parser (argc, argv, &args_info) != 0)
	exit(1);
    
    bool homodimer=args_info.homodimer_given;
    bool antisense=args_info.antisense_given;

    
    // names of molecules for rxns and spcs output
    std::string nameA = "A";
    std::string nameB = "B";
        
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

    bool binary         = args_info.binary_given;
    bool special_open_state = ! args_info.no_special_open_state_given;
    bool gradient       = ! args_info.no_gradient_given;

    consider_double_sites = ! args_info.no_double_sites_given;

    verbose             = args_info.verbose_given;
    debug_out           = args_info.debug_given;

    const std::string inputfile = args_info.input_arg;
    const std::string outputfile = args_info.output_arg;


    bool track = args_info.track_given;
    std::string track_file;
    if (track) { track_file = args_info.track_arg; }
    
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

    const size_t maxsitesize = 
    	(args_info.max_hyb_length_arg>=0)
    	? args_info.max_hyb_length_arg
    	: std::max(seqA.length(),seqB.length());    

    const size_t maxsitesize_diff = 
	(args_info.max_hyb_length_diff_arg>=0)
	? args_info.max_hyb_length_diff_arg
	: std::max(seqA.length(),seqB.length());

    const double max_recover_energy=args_info.max_recover_energy_arg;
    
    size_t region_startA=1;
    size_t region_endA=seqA.length();
    size_t region_startB=1;
    size_t region_endB=seqB.length();
    
    if (args_info.regionA_arg>0) {
	region_endA=std::min(seqA.length(),(size_t)args_info.regionA_arg);
    } else if (args_info.regionA_arg<0) {
	if ((size_t)(-args_info.regionA_arg)<seqA.length()) {
	    region_startA=seqA.length()+args_info.regionA_arg+1;
	}
    }

    if (args_info.regionB_arg>0) {
	region_endB=std::min(seqB.length(),(size_t)args_info.regionB_arg);
    } else if (args_info.regionB_arg<0) {
	if ((size_t)(-args_info.regionB_arg)<seqB.length()) {
	    region_startB=seqB.length()+args_info.regionB_arg+1;
	}
    }

    const size_t span = 
	args_info.span_arg>=0
	? args_info.span_arg
	: std::numeric_limits<size_t>::max();
    
    const size_t window = 
	args_info.span_arg>=0
	? args_info.span_arg*2
	: std::numeric_limits<size_t>::max();


    double min_rate     = args_info.min_rate_arg;

    // if (verbose) {
    // 	std::cerr << "seqA="<<seqA<<", seqB="<<seqB << std::endl;
    // }

    cmdline_parser_free(&args_info);
    
    // global settings for Vienna libRNA
    dangles=2;
    
    if (verbose) {
	std::cerr << "Initialize model." << std::endl;
    }

    stopwatch.start("initialize");

    // construct barrier graph
    RRIBarrierGraph bg(seqA,seqB,
		       special_open_state,
		       maxsitesize,
		       maxsitesize_diff,
		       max_recover_energy,
		       region_startA,
		       region_endA,
		       region_startB,
		       region_endB,
		       span,
		       window,
		       consider_double_sites,
		       gradient,
		       verbose,
		       debug_out);
    
    stopwatch.stop("initialize");

    if (track) {
	bg.track_basins();
    }

    if (verbose) {
	std::cerr << "Construct barrier graph." << std::endl;
    }
    
    stopwatch.start("construct");
    
    try {
	std::ifstream in(inputfile.c_str(),std::ios::in); // | std::ios::binary
	if (!in) {throw std::ifstream::failure("Cannot open file.");}
	bg.read_states(in,binary);
	in.close();
    } catch(std::ifstream::failure &e) {
	std::cerr << "ERROR: Cannot read input from file "<<inputfile<<". "<<e.what()<<std::endl;
	exit(-1);
    }

    stopwatch.stop("construct");
    
    size_t num_total_basins = bg.num_basins();
    

    if (verbose) {
	if (bg.model().is_homodimer()) {
	    std::cerr << "Model homodimer."<<std::endl;
	}
	std::cerr << "Generated "<<num_total_basins<<" basins." << std::endl;
	bg.print_stats(std::cerr);
    }
    
    //stopwatch.start("write");
    if (verbose) {
	std::cerr<<"Write output to "<<outputfile<<std::endl;
    }
    try {
	std::ofstream out(outputfile.c_str(),std::ios::out | std::ios::binary);
	bg.write_binary(out, min_rate);
	out.close();
    } catch(const std::ofstream::failure &e) {
	std::cerr << "Cannot write output to file "<<outputfile<<". "<<e.what()<<std::endl;
    }
    //stopwatch.stop("write");

    if (track) {
	if (track_file.substr(track_file.length()-3,3)==".gz") {
	    // write gzip'd
	    gzFile fh=gzopen(track_file.c_str(),"wb");
	    if (fh==NULL) {
		std::cerr << "Cannot write basin track to file "<<track_file<<"."<<std::endl;
	    } else {
		bg.gzwrite_basin_track(fh);
	    }
	    gzclose(fh);
	} else {
	    // write uncompressed
	    try {
		std::ofstream out(track_file.c_str());
		bg.write_basin_track(out);
		out.close();
	    } catch(const std::ofstream::failure &e) {
		std::cerr << "Cannot write basin track to file "<<track_file<<". "<<e.what()<<std::endl;
	    }
	}
    }

    if (verbose) {
	stopwatch.print_info(std::cerr);
    }

    exit(0);
}
