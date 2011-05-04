#include  <iostream>
#include  <sstream>

#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}


// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_barriers_cmdline.h"

#include <LocARNA/matrices.hh>
#include "hybrid_ensemble_model.hh"

#include <tr1/unordered_map>
typedef std::tr1::unordered_map< std::string, size_t > state_hash_t;

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
    
    std::string seqA=args_info.inputs[0];
    std::string seqB=args_info.inputs[1];
    
    // set some global variables for Vienna libRNA
    dangles=2;
    
    HybEnsModel model(seqA,seqB);  
    
    // read states from standard input
    
    std::string line;
    
    // keys are compressed states!
    // values are basins
    
    state_hash_t state_hash;
    
    size_t state_index=0; // count states
    size_t basin_index=0; // count basins
    
    HybEnsModel::StateDescription state;

    while (getline(std::cin,line)) {
	std::istringstream linein(line);
	double energy;
	linein >> energy;
	
	std::vector<size_t> state_vec;
	for (size_t i; linein >> i;) {
	    state_vec.push_back(i);
	}
	
	switch ( state_vec.size() ) {
	case 0:
	    state = HybEnsModel::StateDescription();
	    break;
	case 4:
	    state = HybEnsModel::StateDescription(state_vec[0],
						  state_vec[1],
						  state_vec[2],
						  state_vec[3]);
	    break;
	case 8:
	    state = HybEnsModel::StateDescription(state_vec[0],
						  state_vec[1],
						  state_vec[2],
						  state_vec[3],
						  state_vec[4],
						  state_vec[5],
						  state_vec[6],
						  state_vec[7]);
	    break;
	default:
	    std::cerr << "Warning: ignored invalid input line" << std::endl;
	    continue;
	}
	
	//std::cout << "read " << state_index << " " << energy << " " << " "  << state << std::endl;
	
	
	// enumerate moves/neighbors of the state
	
	std::vector<size_t> neighbor_basins; // already seen neighbor basins of state 
	
	HybEnsModel::MoveIterator mi(state,model);
	for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    
	    
	    //std::cout << " move "; move->print(std::cout); std::cout<<std::endl;
		
	    HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	    if (tE < 1e6) { // if there is a valid tranisition
		
		HybEnsModel::StateDescription state2=state;
		move->apply(state2);
		
		std::string code; // string for holding codes
		state2.encode(code);
		
		//std::cout << " neighbor "<<state2<<std::endl;

		// // test encoding again
		// HybEnsModel::StateDescription state3;
		// state3.decode(code);
		// //std::cout << "Test " << state2 << " " << state3<< std::endl;
		// assert(state2==state3);
		
		state_hash_t::const_iterator it = state_hash.find(code);
		if (state_hash.end() != it) {
		    //found => belongs to basin of already seen local minimum
		    neighbor_basins.push_back(it->second);
		}
	    }
	}
	
	if (neighbor_basins.size()==0) {
	    //not found => new local minimum
	    basin_index++;
	    std::cout << state << " is local minimum " << basin_index << std::endl;
	    
	    // put state into hash
	    state_hash[state.encode()] = basin_index;
	} else {
	    // state connects into at least one basin
	    
	    std::cout << state << " connects into basin(s) ";
	    
	    for (std::vector<size_t>::iterator it=neighbor_basins.begin(); it!=neighbor_basins.end(); it++) {
		std::cout << *it << " ";
	    }
	    std::cout << std::endl;
	
	    // need to determine steepest descent neighbor or some other criterion to assign
	    // the state to one of the neighbor basins
	    // ATTENTION: arbitrarily assign to the first in the list!
	    
	    state_hash[state.encode()] = neighbor_basins[0];
	    
	}
	
	state_index++;
    }

    cmdline_parser_free(&args_info);    
    exit(0);
}
