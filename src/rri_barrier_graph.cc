#include "rri_barrier_graph.hh"


#include <iostream>
#include <iomanip>

const bool RECOVER_MISSING_STATES=true;

bool 
RRIBarrierGraph::read_state(std::istream &in, 
			 HybEnsModel::StateDescription &state,
			 double &energy, 
			 size_t lineno,
			 bool binary) const {
    
    if(binary) {
	in >> energy;
	if (in.eof()) return false;
	
	if(in.get()!=' ') {
	    std::cerr << "expected blank after energy at line "<< lineno <<std::endl;
	    return false; 
	}
	
	if ( !state.read_binary(in) || in.get()!=0) { 
	    std::cerr << "Error while parsing state at "<< lineno <<"."<<std::endl;
	    return false; 
	}
	
    }else{
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
	
    }

    if (not state.is_valid(model_)) {
	std::cerr << "ERROR: read state "<<state<<" at line "<<lineno<<" is not valid in model."<<std::endl;
	exit(-1);
    }
        
    return true;
}

void
RRIBarrierGraph::process_move(const HybEnsModel::Move *move,
			      const HybEnsModel::StateDescription &source_state,
			      const HybEnsModel::energy_t &source_energy,
			      HybEnsModel::energy_t &min_transition_energy,
			      size_t &min_transE_target_basin_index,
			      std::vector<transition_t> &transitions
			      ) {
    //std::cout << " move "; move->print(std::cout); std::cout<<std::endl;
        
    HybEnsModel::energy_t transE=move->transitionEnergy(); // -- energy of transition state
    
    HybEnsModel::StateDescription neigh_state=source_state;
    move->apply(neigh_state);
    
    if (!consider_double_sites_ && neigh_state.size()==2) {
	return;
    }
    
    // encode neighbor and search neighbor code in hash
    
    HybEnsModel::StateDescription::code_t neigh_code; // string for holding code
    neigh_state.encode(neigh_code);
    
    state_hash_t::const_iterator it = state_hash_.find(neigh_code);
    if (state_hash_.end() != it) {
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
	if ( transE < min_transition_energy ) {
	    // record new minimal energy transition
	    min_transition_energy = transE;
	    min_transE_target_basin_index = target_basin_index;
	}
	
	// compute energy of neighbor state
	double neigh_energy = model_.energy(neigh_state);
	
	// record new transition.  In this way, we collect all
	// transitions from the source state to energetically lower
	// target states (recall: input states are sorted by energy).
	// Notably, we don't record the transitions to higher energy
	// states (recall: these states are not in the hash). Since
	// transitions are symmetric, we add those transitions later.
	transitions.push_back(transition_t(source_energy,
				     target_basin_index,
				     neigh_energy,
				     transE)
			);
	
	if (debug_out_) {
	    std::cerr <<"\t"<< transitions[transitions.size()-1]
		      << " by ";
	    move->print(std::cerr);
	    std::cerr << std::endl;
	}
    }
}

void
RRIBarrierGraph::assign_to_basin(const HybEnsModel::StateDescription &state,
				 const HybEnsModel::energy_t &energy,
				 size_t basin_index
				 ) {
    if (debug_out_) std::cerr << "  Assign to basin "<<basin_index<<std::endl;
    
    // assign basin index source_basin_index to source_state and register
    // source_state as new member of the basin
    state_hash_[state.encode()] = basin_index;
    basins_[basin_index].add_state(model_.boltzmann_weight(energy));
}

size_t
RRIBarrierGraph::create_new_basin(const HybEnsModel::StateDescription &state,
				  const HybEnsModel::energy_t &energy
				  ) {
    size_t basin_index = basins_.size();
    
    if (debug_out_) std::cerr << "  New basin "<<basin_index<<std::endl;
    
    // put state into hash
    state_hash_[state.encode()] = basin_index;
    
    // generate new basin and put into object's basin list
    basins_.push_back(Basin(basin_index,model_.boltzmann_weight(energy)));
    
    return basin_index;
}

void 
RRIBarrierGraph::register_transitions(size_t source_basin_index,
				      std::vector<transition_t> &trans
				      ) {
    // ----------------------------------------
    // Register all transitions from source_basin_index to other basins.
    //
    for(std::vector<transition_t>::iterator it=trans.begin(); 
	trans.end()!=it; ++it) {
	if (it->target_basin_index != source_basin_index) {
		
	    it->source_basin_index = source_basin_index;
		
	    if (debug_out_)
		std::cerr << "add transition "<<source_basin_index<<" <-> "
			  << it->target_basin_index << " "
		    //<< it->transition_energy << " "
			  << it->barrier_energy() << " "
			  << std::endl;
	

	    // note: each transition between source and target is
	    // registered only once, namely when the source index is
	    // larger, therefore we add transitions in both directions
	    // (finally, this produces a symmetric matrix <transitions>
	    add_transition(*it);
	    add_transition(it->reverse());
	    
	}
    } // end iterate trans
}

void
RRIBarrierGraph::process_state(const HybEnsModel::StateDescription &source_state, 
			       double source_energy) {
    
    // minimum transition state energy from source state source_state to a neighbor
    double min_transition_energy = std::numeric_limits<double>::infinity();
	
    // index of the basin of the neighbor state with minimum transition energy
    size_t min_transE_neighbor_basin_index = std::numeric_limits<size_t>::max();

    std::vector<transition_t> transitions; //! transitions from source to neighbor states

    size_t moves_counter=0;

    // ----------------------------------------
    // Enumerate neighbors of source_state
    //
    // Register all transitions from source_state to previously read neighbor states.
    // THIS IMPLIES that each transition between x and y is registered only once.
    // 
    // Keep track of neighbor state with smallest transition energy. 
    //
    //
    HybEnsModel::MoveIterator mi(source_state,model_,consider_double_sites_);
    for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	moves_counter++;

	process_move(move,
		     source_state,
		     source_energy,
		     min_transition_energy,       // <-- all following params are output parameters!
		     min_transE_neighbor_basin_index,
		     transitions
		     );

    }
   
    if (debug_out_) std::cerr << "  " << transitions.size() << " transitions, "
			     << moves_counter << " moves" <<std::endl;
    

    if (moves_counter==0) {
	if (debug_out_) std::cerr << "Ignore frozen state " << source_state << "." << std::endl;
	return;
    }

    // index of the basin of the source state
    size_t basin_index=std::numeric_limits<size_t>::max();
 
    // in case of homodimer check for symmetry, and assign to basin of previously found symmetric state 
    if (model_.is_homodimer()) {
	HybEnsModel::StateDescription symmetric_state = source_state.symmetric_state(model_.seqA().length());
	
	if (symmetric_state == source_state) {
	    // if state is symmetric to itself, then do nothing
	    if (debug_out_) {
		std::cerr << "  State " << model_.to_dotbracket(source_state)
			  <<" is self-symmetric."
			  <<std::endl;
	    }
	} else {
	    // Test whether the symmetric state of source_state has been seen before.
	    // In this case, the two symmetric states are identified.
	    state_hash_t::const_iterator symmetric_state_it=state_hash_.find(symmetric_state.encode());
	    if (symmetric_state_it != state_hash_.end()) {
		if (debug_out_) {
		    std::cerr << "  State " << model_.to_dotbracket(source_state)
			      <<" is equivalent to previously found "
			      << model_.to_dotbracket(symmetric_state)
			      <<" of basin "<<symmetric_state_it->second<<"."<<std::endl;
		}
		
		basin_index = symmetric_state_it->second;
		if (debug_out_) std::cerr << "  Assign to basin "<<basin_index<<std::endl;
		state_hash_[source_state.encode()] = basin_index;
		basins_[basin_index].add_state(model_.boltzmann_weight(source_energy));
	    }
	}
    } // end special homodimer symmetry handling
    
    
    // assign basin
    if ( basin_index == std::numeric_limits<size_t>::max()) { //unless assigned by symm handling
	
	// assign source state to either new basin or the basin of its deepest neighor; the latter, if barrier<=0
	
	if ( !gradient_ || min_transition_energy > source_energy ) {
	    // no transition state to an energetically lower target
	    // state is energetically lower than the source state
	    //
	    // Consequently, source_state is a new local minimum (among
	    // the previously seen states)

	    
	    // with state recovery: check whether there is a gradient walk;
	    // if yes, create missing states on walk
	    // otherwise, create new basin
	    if (RECOVER_MISSING_STATES 
		&& !(special_first_state_ && source_state.size()==0)
		&& source_energy <= max_recover_energy_
		) {
		basin_index = create_gradient_walk(source_state,source_energy);
	    } else {
		basin_index = create_new_basin(source_state,source_energy);
	    }
	} else {
	    // source_state is not a local minimum but belongs to basin
	    // min_transE_target_basin_index, which is the basin that is reached
	    // with the lowest transition energy.
	    
	    assign_to_basin(source_state,source_energy,min_transE_neighbor_basin_index);
	    basin_index = min_transE_neighbor_basin_index;
	}
    }
    
    // register the transitions from the source basin to other states in the hash
    register_transitions(basin_index,transitions);

}

size_t
RRIBarrierGraph::create_gradient_walk(const HybEnsModel::StateDescription &state,
				      double energy) {
    
    if (debug_out_) 
	std::cerr << "Walk along from "<<state<<" "<<energy<<std::endl;
    
    HybEnsModel::StateDescription::code_t state_code = state.encode();
    state_hash_t::const_iterator it = state_hash_.find(state_code);

    if (state_hash_.end() != it) {
	// hey, surprise! the state is already hashed

	if (debug_out_) 
	    std::cerr << "... hits the hash at "<<it->second<<" " << state << "."<<std::endl;
	
	// do nothing
	return it->second;
    }
    
    // ------------------------------------------------------------
    // Determine the minimum transition energy neighbor
    //
    
    // minimum transition state energy from source state state to a neighbor
    double min_transE = std::numeric_limits<double>::infinity();
    
    // index of the basin of the neighbor state with minimum transition energy
    HybEnsModel::StateDescription gradient_state;
    // ... and the energy of this state
    double gradient_state_energy = std::numeric_limits<double>::infinity();
    
    size_t moves_counter=0;

    // Enumerate neighbors of state to get minimum transition
    // energy neighbor
    HybEnsModel::MoveIterator mi(state,model_,consider_double_sites_);
    for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	moves_counter++;
	
        HybEnsModel::MoveIterator mi(state,model_,consider_double_sites_);
	for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	    HybEnsModel::StateDescription neigh_state = state;
	    move->apply(neigh_state);
	    
	    double neigh_energy=model_.energy(neigh_state);
	    
	    if (tE<min_transE && !isinf(neigh_energy)) {
		min_transE = tE;
		
		gradient_state = neigh_state;
		gradient_state_energy = neigh_energy;
	    }
	}
    }

    if (debug_out_)  std::cerr << "... min tE: "<<min_transE-energy<<std::endl;
    //
    // ------------------------------------------------------------
    
    size_t basin_index;
    if (min_transE >= energy) {
	// the source state is a local minimum, create basin (and terminate)
	
	if (debug_out_)
	    std::cerr << "... create new basin at "<<state<<": "<< energy <<std::endl;
	
	basin_index = create_new_basin(state,energy);
	
    } else {
	// there is a gradient neighbor, walk along to find the final basin index
	
	// recursive call
	basin_index = create_gradient_walk(gradient_state,gradient_state_energy);
	
	basins_[basin_index].add_state(model_.boltzmann_weight(energy));
	state_hash_[state_code] = basin_index;
    }
    
    return basin_index;
}

void
RRIBarrierGraph::read_states(std::istream &in,
			     bool binary) {

    HybEnsModel::StateDescription source_state;
    double energy;
        
    // counter for states that are read for construction of graph
    size_t state_counter;

    state_counter=0;
    
    // use to check input
    double last_energy=-std::numeric_limits<double>::infinity();
    int line=1;
    
    if (special_first_state_) {
	// add the empty/open state first, to guarantee to create a new basin 

	HybEnsModel::StateDescription open_state;

	if (debug_out_) { 
	    std::cerr << "add open state " << open_state << std::endl;
	}
	
	process_state(open_state, -4.1);
	
	state_counter++;
    }
    

    while (read_state(in,source_state,energy,line,binary)) {
    	
	if (verbose_ && (state_counter%5000==0)) {
	    std::cerr << "\r" << state_counter<< " \tb:"<<basins_.size()<<" \th:"<<state_hash_.size()<<"\t";
	}
	if (energy < last_energy) {
	    std::cerr << "ERROR: input states have to be sorted by increasing energy (at line "
		      <<line<<": "<<energy<<"<"<<last_energy << " ).\n";
	    exit(-1);
	}

	if (debug_out_) {
	    std::cerr << "read " << state_counter << " " << energy << " "
		      << " "  << source_state << std::endl;
	}
	
	if (special_first_state_ && source_state.size()==0) {
	    // if (verbose_) {
	    // 	std::cerr << std::endl << "Ignore open state in input." <<std::endl;
	    // }
	    continue;
	}
	
	if (!consider_double_sites_ && source_state.size()==2) {
	    if (debug_out_) {
		std::cerr << "Ignore double site state " << source_state <<"."<<std::endl;
	    }
	    continue;
	}
	
	process_state(source_state,energy);
            
	state_counter++;
	
	last_energy=energy;
	line++;
    }
    std::cerr << "\r";

}


RRIBarrierGraph::RRIBarrierGraph(const std::string &seqA, 
				 const std::string &seqB,
				 bool special_open_state,
				 size_t maxsitesize,
				 size_t maxsitesizediff,
				 double max_recover_energy,
				 bool consider_double_sites,
				 bool gradient,
				 bool verbose,
				 bool debug_out)
    :
    BarrierGraph(special_open_state,verbose,debug_out),
    model_(seqA,seqB,std::max(seqA.length(),seqB.length()),maxsitesizediff,consider_double_sites),
    max_recover_energy_(max_recover_energy),
    consider_double_sites_(consider_double_sites),
    gradient_(gradient)
{
        
}


void
RRIBarrierGraph::add_transition( const transition_t &tr ) {
    // add transition to partition function for the transition
    // between the source and target basin
    if (transitions_[tr.source_basin_index].find(tr.target_basin_index)
	== transitions_[tr.source_basin_index].end()) {
	transitions_[tr.source_basin_index][tr.target_basin_index] = 
	    BasinTransition( model_.boltzmann_weight(tr.barrier_energy()));
    } else {
	transitions_[tr.source_basin_index][tr.target_basin_index].update(model_.boltzmann_weight(tr.barrier_energy()));
    }
}


/* Output microstate transitions */
std::ostream &
operator <<(std::ostream &out,const RRIBarrierGraph::transition_t &t) {
    out << "Transition to "<<t.target_basin_index
	<< " (sourceE="<<t.source_state_energy
	<< " tgtE="<<t.target_state_energy
	<< " transE="<<t.transition_energy
	<< " transE-sourceE="<<(t.transition_energy-t.source_state_energy)<<")";
    return out;
}
