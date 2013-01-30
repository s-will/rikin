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
 *
 * @todo CHECK: is the process irreducibel?
 */

#include "rrikin_barriers.hh"

#include  <sstream>

#include  <stdlib.h>
#include  <string>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_barriers_cmdline.h"

/* Methods of basin */

void Basin::print_header(std::ostream &out) const {
    printf("%5s %-30s %4s %6s %5s %6s %6s\n",
	   "idx",
	   "description",
	   "n_s",
	   "ensE",
	   "ct",
	   "minE",
	   "h"
	   );
}

void
Basin::print(std::ostream &out, const HybEnsModel &model) const {
    printf("%5lu %-30s %6.2f %6.2f %5lu %6.2f %6.2f\n",
	   basin_index,
	   local_minimum.toString().c_str(),
	   states,
	   - model.RT() * log(Z),
	   connect_to,
	   minimum_energy,
	   barrier-minimum_energy
	   );
}

/* Output transitions */

std::ostream &
operator <<(std::ostream &out,const BarrierGraph::transition_t &t) {
    out << "Transition to "<<t.target_basin_index
	<< " (sourceE="<<t.source_state_energy
	<< " tgtE="<<t.target_state_energy
	<< " transE="<<t.transition_energy
	<< " transE-sourceE="<<(t.transition_energy-t.source_state_energy)<<")";
    return out;
}

/* Methods of BarrierGraph */

bool 
BarrierGraph::read_state(std::istream &in, 
			 HybEnsModel::StateDescription &state,
			 double &energy, 
			 size_t lineno) const {
    std::string line;
	
    if (!getline(in,line)) return false;
	
    std::istringstream linein(line);
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
	std::cerr << "ERROR: invalid input line "<< lineno<<"." << std::endl;
	exit(-1);
    }
        
    if (not state.is_valid(model)) {
	std::cerr << "ERROR: read state "<<state<<" at line "<<lineno<<" is not valid in model."<<std::endl;
	exit(-1);
    }
        
    return true;
}


void
BarrierGraph::add_transition( const transition_t &tr, const HybEnsModel &model) {
    Basin &src_basin = basins[tr.source_basin_index];
	
    // keep track of minimal energy transition from this basin to some smaller basin
    if ( (tr.source_basin_index > tr.target_basin_index) && tr.barrier_energy() < src_basin.barrier ) {
	if (debug_out) std::cerr << "  New barrier for basin "<<tr.source_basin_index<<" to basin "<< tr.target_basin_index  <<": "<< tr.barrier_energy()<<std::endl;
	src_basin.barrier = tr.barrier_energy();
	src_basin.connect_to = tr.target_basin_index;
    }
	
    // add transition to partition function for the transition
    // between the source and target basin
    if (transitions[tr.source_basin_index].find(tr.target_basin_index)
	== transitions[tr.source_basin_index].end()) {
	transitions[tr.source_basin_index][tr.target_basin_index] = BasinTransition(tr.barrier_energy(),model);
    } else {
	transitions[tr.target_basin_index][tr.target_basin_index].update(tr.barrier_energy(),model);
    }
}
    

void
BarrierGraph::process_state(const HybEnsModel::StateDescription &source_state, double source_energy) {
	
    // minimum transition state energy from source state source_state to a neighbor
    double min_transition_energy = std::numeric_limits<double>::infinity();
	
    // index of the basin of the neighbor state with minimum transition energy
    size_t min_tE_target_basin_index = std::numeric_limits<size_t>::max();

    std::vector<transition_t> trans; //! transitions from source to neighbor states

    size_t moves_counter=0;

    // ----------------------------------------
    // Enumerate neighbors of source_state
    //
    // Register all transitions from source_state to previously read neighbor states.
    // Keep track of neighbor state with smallest transition energy 
    //
    HybEnsModel::MoveIterator mi(source_state,model);
    for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    
	//std::cout << " move "; move->print(std::cout); std::cout<<std::endl;
	    
	moves_counter++;
	    
	HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	HybEnsModel::StateDescription neigh_state=source_state;
	move->apply(neigh_state);
	    
	// #ifndef NDEBUG
	// 	    // check whether generated neighbor state is a valid state in the model
	// 	    if( not neigh_state.is_valid(model) ) {
		
	// 		std::cerr
	// 		    << "ERROR: Move produces invalid state" << std::endl 
	// 		    << source_state << " == ";
	// 		move->print(std::cerr); 
	// 		std::cerr << " ==> " << neigh_state << std::endl;
	// 		abort();
	// 	    }
	// #endif
	    
	// encode neighbor and search neighbor code in hash
	    
	std::string neigh_code; // string for holding code
	neigh_state.encode(neigh_code);
	    
	state_hash_t::const_iterator it = state_hash.find(neigh_code);
	if (state_hash.end() != it) {
	    //found => belongs to basin of already seen local minimum
		
	    // NOTE: in the transition from source_state to
	    // neigh_state, neigh_state has lower or equal energy
	    // than the source_state! The transition state has higher
	    // or equal energy than neigh_state (and possibly source_state)
		
	    // NOTE: the following differs from the standard
	    // barriers algorithm, since there transition energies
	    // are sorted in the same order as target state
	    // energies
		
	    size_t target_basin_index=it->second;
		
	    // if transition energy to neigh_state is smaller than the former minimum
	    if ( tE < min_transition_energy ) {
		// record new minimal energy transition
		min_transition_energy = tE;
		min_tE_target_basin_index = target_basin_index;
	    }
		
	    // compute energy of neighbor state
	    double neigh_energy = model.energy(neigh_state);
		
	    // record new transition.  In this way, we collect all
	    // transitions from the source state to energetically
	    // lower target states Explicitly, we don't record the
	    // transitions to higher energy states. Since
	    // transitions are symmetric, we can later add those
	    // transitions.
	    trans.push_back(transition_t(source_energy,      // energy of source state
					 target_basin_index, // index of target basin
					 neigh_energy,       // energy of target state
					 tE)                 // energy of the transition state
			    );
		
	    if (debug_out) {
		std::cerr <<"\t"<< trans[trans.size()-1]
			  << " by ";
		move->print(std::cerr);
		std::cerr << std::endl;
	    }
	}
    }

    if (debug_out) std::cerr << "  " << trans.size() << " transitions, "
			     << moves_counter << " moves" <<std::endl;
	


    if (moves_counter==0) {
	if (debug_out) std::cerr << "Ignore frozen state"<<std::endl;
	return;
    }

    // index of the basin of the processed source state
    size_t source_basin_index;
	
    // ------------------------------------------------------------
    // Perform basin assignment of source_state.
    //
    // Either construct a new basin with local minimum source_state
    // or assign source_state to an existing basin
    //
    if ( min_transition_energy > source_energy ) {
	// no transition state to an energetically lower target
	// state is energetically lower than the source state
	//
	// Consequently, source_state is a new local minimum
	    
	source_basin_index = basins.size();
	    
	if (debug_out) std::cerr << "  New basin "<<source_basin_index<<std::endl;

	// put state into hash
	state_hash[source_state.encode()] = source_basin_index;
	    
	// generate new basin and put into object's basin list
	Basin new_basin(source_basin_index,source_state,source_energy,model);
	basins.push_back(new_basin);
	    
    } else {
	// source_state is not a local minimum but belongs to basin
	// min_tE_target_basin_index, which is the basin that is reached
	// with the lowest transition energy.
	    
	// handle case where basin with index source_basin_index
	// was merged before
	source_basin_index = min_tE_target_basin_index;

	if (debug_out) std::cerr << "  Assign to basin "<<source_basin_index<<std::endl;
	    
	// assign basin index source_basin_index to source_state and register
	// source_state as new member of the basin
	state_hash[source_state.encode()] = source_basin_index;
	basins[source_basin_index].add_state(source_energy,model);
    }
	
    // ----------------------------------------
    // Register all transitions from source_basin_index to other basins.
    //
    for(std::vector<transition_t>::iterator it=trans.begin(); 
	trans.end()!=it; ++it) {
	if (it->target_basin_index != source_basin_index) {
		
	    it->source_basin_index = source_basin_index;
		
	    if (debug_out)
		std::cerr << "add transition "<<source_basin_index<<" -> "
			  << it->target_basin_index << " "
			  << it->transition_energy << " "
			  << it->barrier_energy() << " "
			  << std::endl;
		
	    add_transition(*it,model);
	    add_transition(it->reverse(),model);
		
	}
    }	
}

size_t
BarrierGraph::num_basins() const {
    size_t n=0;
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    n++;
	}
    }
    return n;
}


size_t
BarrierGraph::num_transitions() const {
    size_t n=0;
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    const transitions_map_t::const_iterator &it = transitions.find(i);
	    n += it->second.size();
	}
    }
    return n;
}

void
BarrierGraph::merge_basins(double max_outflow, double min_rate) {
	
    // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
    
    
    // run through sorted basins and merge
    for (size_t i=0; i<basins.size(); ++i) {
	
	Basin &x0 = basins[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions[x0.basin_index];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.basin_index)) {
	    trs_x0.erase(x0.basin_index);
	}

	// remove transitions wiht rates lower than min_rate
	for (transitions_map_row_t::const_iterator it=trs_x0.begin();
	     trs_x0.end()!=it; ) {
	    
	    //if (basins[it->first].merged()) continue;
	    assert(!basins[it->first].merged());
	    
	    if (it->second.get_Z() / std::min(x0.get_Z(), basins[it->first].get_Z())
		< min_rate) {
		
		//remove transition
		//std::cerr << "Remove transition between "<<x0.basin_index<<" and "<<it->first<<std::endl; 
		transitions[it->first].erase(x0.basin_index);
		it=trs_x0.erase(it);
		
	    } else {
		++it;
	    }
        }
	//std::cerr << "DONE"<<std::endl; 
		
	
        // compute total outflow
        double total_out = 0.0;
        for (transitions_map_row_t::const_iterator it=trs_x0.begin();
	     trs_x0.end()!=it; ++it) {
	    
	    if (basins[it->first].merged()) continue;
	    if (it->first==x0.basin_index) continue;
	    
	    total_out += it->second.get_Z();
        }

	// first filter outgoing rates



	    
	if (debug_out) {
	    std::cerr << i << " " << sorted_basin_idxs[i] << " " << x0.get_Z() << " " << total_out << " r=" << (total_out/x0.get_Z());
	    
	    std::cerr <<"  \t";
	}
        for (transitions_map_row_t::const_iterator it=trs_x0.begin();
	     trs_x0.end()!=it; ++it) {

	    if (basins[it->first].merged()) continue;
	    
	    if (debug_out) {
		std::cerr << " r_"<<it->first<<"="<< (it->second.get_Z()/x0.get_Z());
	    }
        }
	    
	if (debug_out) {
	    std::cerr << std::endl;
	}

        if (total_out/x0.get_Z() > max_outflow) {

	    if (debug_out) {
		std::cerr << "Dissolve basin "<< x0.basin_index << std::endl;
	    }

	    for (transitions_map_row_t::iterator it=trs_x0.begin();
		 trs_x0.end()!=it; ++it) {
		
		Basin &y = basins[it->first];
		if (y.merged()) continue;
		
		if (x0.basin_index==y.basin_index) continue;
		
		// 1) distribute the basin's partition function to its neighbors
		double Z_yx0 =  it->second.get_Z();
		double fraction = Z_yx0/total_out;
		
		y.merge_in(x0,fraction);
		
		if (debug_out) {
		    std::cerr << "Transfer a fraction of " << fraction << " of " << x0.basin_index
			      << "'s partfunc to " << y.basin_index<<std::endl;
		}
		
		// 2) distribute the partition function of the transitions to this basin to its neighbors

		for (transitions_map_row_t::iterator it2=trs_x0.begin();
		     trs_x0.end()!=it2; ++it2) {
		    
		    Basin &x = basins[it2->first];
		    if (x.merged()) continue;
				    
		    if (x.basin_index==y.basin_index) continue;
		    if (x0.basin_index==x.basin_index) continue;
		    
		    // update the transition from x to y (via x0)
		    
		    double Z_xx0 = transitions[x.basin_index][x0.basin_index].get_Z();
		    
		    if (debug_out) {
			std::cerr << "Transfer a fraction of " << fraction << " of transition "
				  << x.basin_index <<"->"<< x0.basin_index << " to " << x.basin_index 
				  << "->"<< y.basin_index<<std::endl;
		    }
		    
		    transitions[x.basin_index][y.basin_index].update(Z_xx0 * fraction);		    
		} // end for it2 (over neighbors of x0)
		
	    } // end for it (over neighbors of x0)
	    
	    // finally, mark as merged
	    if (debug_out) {
		std::cerr << "Mark "<<x0.basin_index<<" as merged."<<std::endl;
	    }
	    x0.merged_ = true;
	    
	    // and release all transitions from and to the merged basin
	    for (transitions_map_row_t::iterator it=trs_x0.begin();
		 trs_x0.end()!=it; ++it) {
		transitions[it->first].erase(x0.basin_index);
	    }
	    transitions.erase(x0.basin_index);
	    
        }
	
    }
}

void
BarrierGraph::filter_rates(double min_rate) {
    // run through basins
    for (size_t i=0; i<basins.size(); ++i) {
	
	Basin &x0 = basins[i];
	
	if (x0.merged()) continue;
		
	transitions_map_row_t &trs_x0=transitions[i];
	
	// assert((trs_x0.end()==trs_x0.find(i)));
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(i)) {
	    //std::cerr << "Basin "<<i<<" has transition to itself: " <<transitions[i][i].get_Z() << "."<<std::endl;
	    trs_x0.erase(i);
	}
	
	// remove transitions wiht rates lower than min_rate
	for (transitions_map_row_t::const_iterator it=trs_x0.begin();
	     trs_x0.end()!=it; ) {
	    
	    //if (basins[it->first].merged()) continue;
	    assert(!basins[it->first].merged());
	    
	    if (it->second.get_Z() / std::min(x0.get_Z(), basins[it->first].get_Z())
		< min_rate) {
		
		//remove transition
		//std::cerr << "Remove transition between "<<i<<" and "<<it->first<<std::endl; 
		transitions[it->first].erase(i);
		it=trs_x0.erase(it);
		
	    } else {
		++it;
	    }
        }
	//std::cerr << "DONE"<<std::endl; 
    }    
}


BarrierGraph::BarrierGraph(const std::string &seqA, const std::string &seqB)
    : model(seqA,seqB)
{
        
    HybEnsModel::StateDescription orig_state;
    double energy;
        
    // counter for states that are read for construction of graph
    size_t state_counter;

    state_counter=0;
    
    // use to check input
    double last_energy=-std::numeric_limits<double>::max();
    int line=1;

    while (read_state(std::cin,orig_state,energy,line)) {
	if (energy < last_energy) {
	    std::cerr << "ERROR: input states have to be sorted by increasing energy (at line "<<line<<": "<<energy<<"<"<<last_energy << " ).\n";
	    exit(-1);
	}

	if (debug_out) std::cerr << "read " << state_counter << " " << energy << " " << " "  << orig_state << std::endl;
            
	process_state(orig_state,energy);
            
	state_counter++;

	last_energy=energy;
	line++;
    }
}

void
BarrierGraph::print_edges(std::ostream &out,const Basin &b) const {
    out << b.basin_index
	//<< " (" << b.get_Z() << ") "
	<< "  -> ";
	
    transitions_map_t::const_iterator
	ts = transitions.find(b.basin_index);
	
    if (ts == transitions.end()) return;
	
    for(transitions_map_row_t::const_iterator it=ts->second.begin();
	it != ts->second.end();
	++it) {
	    
	if ((it->first!=b.basin_index) && (!basins[it->first].merged())) {
	    out << " " << it->first
		//<< " (" << it->second.get_Z() << ")"
		<< " " << (it->second.get_Z() / b.get_Z());
	}
    }
    out << std::endl;
}

void 
BarrierGraph::print_barrier_graph(std::ostream &out) const {
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    print_edges(out,basins[i]);
	}
    }	
}


void
BarrierGraph::print_basins(std::ostream &out) const {
    basins[0].print_header(out);
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    basins[i].print(out,model);
	} else {
	    // printf("merged to %4lu",basins[i].merged_to); basins[i].print(std::cout,model);
	}
    }
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
    
    if ( args_info.inputs_num != 2 ) {
	std::cerr << "Expect two sequences as input on command line"<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }
    
    std::string seqA=args_info.inputs[0];
    std::string seqB=args_info.inputs[1];

    cmdline_parser_free(&args_info);
    
    // global settings for Vienna libRNA
    dangles=2;
    
    if (verbose) {
	std::cerr << "Construct barrier graph." << std::endl;
    }

    stopwatch.start("construct");

    // construct barrier graph
    BarrierGraph barriers(seqA,seqB);

    stopwatch.stop("construct");
    
    size_t num_total_basins = barriers.num_basins();
    

    double max_outflow=0.7;
    double min_rate=1e-5;    
    
    if (verbose) {
	std::cerr << "Generated "<<num_total_basins<<" basins." << std::endl;
	stopwatch.print_info(std::cerr);
	std::cerr << "Merge basins with outflow larger " << max_outflow << std::endl;
	std::cerr << "Remove rates smaller than " << min_rate << std::endl;
	
	barriers.print_stats(std::cerr);
    }

    stopwatch.start("merge");
    
    barriers.merge_basins(max_outflow,min_rate);
    
    stopwatch.stop("merge");

    if (verbose) {
	std::cerr << "Merged "<<num_total_basins-barriers.num_basins()
		  <<" basins resulting in "<<barriers.num_basins()<<" states." << std::endl;
	barriers.print_stats(std::cerr);
    }

    if (verbose) {
	std::cerr << "Remove rates smaller than " << min_rate << " (again)."<< std::endl;
    }
    
    barriers.filter_rates(min_rate);

    if (verbose) {
	barriers.print_stats(std::cerr);
    }
    
    // print basins of barrier graph
    barriers.print_basins(std::cout);    

    std::cout << std::endl
	      << std::endl;
    barriers.print_barrier_graph(std::cout);

    if (verbose) {
	stopwatch.print_info(std::cerr);
    }    

    exit(0);
}
