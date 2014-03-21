#include "basin_transition.hh"
#include "barrier_graph.hh"

#include <stack>
#include <iostream>
#include <stdio.h>
#include <iomanip>
#include <limits>
#include <cmath>
#include <set>
#include <sstream>


const bool GRAD_WALK_MISSING=true; //whether to complete gradient walks for missing states

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
			 size_t lineno,
			 bool binary) const {
    
    if(binary) {
    in >> energy;
    if (in.eof()) return false;

    if(in.get()!=' ') {
	std::cerr << "expected blank after energy at line "<< lineno <<std::endl;
	return false; 
    }
    
    HybEnsModel::StateDescription::code_t code;
    
    char *codebuf = reinterpret_cast<char *>(&code);

    in.read(codebuf,8);
    
    // std::cerr << lineno << " " << energy << std::endl;
    // for (size_t i=0; i<8; i++) {
    // 	std::cerr << (int)codebuf[i] << " " ;
    // }
    // std::cerr << std::endl;
    
    if(in.get()!=0) {
	std::cerr << "expected zero delimiter after state code at line "<< lineno <<std::endl;
	return false; 
    }
    
    state.decode(code);
    
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
BarrierGraph::add_transition( const transition_t &tr, const HybEnsModel &model) {
    // add transition to partition function for the transition
    // between the source and target basin
    if (transitions_[tr.source_basin_index].find(tr.target_basin_index)
	== transitions_[tr.source_basin_index].end()) {
	transitions_[tr.source_basin_index][tr.target_basin_index] = BasinTransition(tr.barrier_energy(),model);
    } else {
	transitions_[tr.source_basin_index][tr.target_basin_index].update(tr.barrier_energy(),model);
    }
}
    

void
BarrierGraph::process_move(const HybEnsModel::Move *move,
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
    } else if (GRAD_WALK_MISSING
	       &&
	       neigh_state.max_site_size() > maxsitesize_
	       && !isinf(model_.energy(neigh_state))
	       ) {
	// Attention: there are moves (shift moves) where the target state is invalid, but
	// the transition energy is finite!

	
	// The neighbor state is not in the hash, because neighbor site length is beyond maxsitesize

	// compute energy of neighbor state
	
	// skip open state if open state is special
	if (!(special_open_state_ && source_state.size()==0)) {

	    // if transition energy to neigh_state is smaller than the former minimum
	    if ( transE < min_transition_energy) {
		
		// std::cerr
		//     << ((model_.energy(neigh_state)>=source_energy)?"* ":"  ")
		//     << source_state << " " << source_energy 
		//     << " " << neigh_state
		//     << " " << model_.energy(neigh_state) << " " << transE << std::endl;
		
		// record new minimal energy transition
		min_transition_energy = transE;
		min_transE_target_basin_index = std::numeric_limits<size_t>::max(); // signal mimimum not in hash
	    
		// transitions to states beyond maxsitesize are not recorded
	    }
	}
    }
    // else {
    // 	if (debug_out_ && !(special_open_state_ && source_state.size()==0) && model_.energy(neigh_state)<source_energy) {
    // 	    std::cerr << "Neighbor "<< neigh_state<<" of "<<source_state<<" is missing."<<std::endl; 
    // 	}
    // }
}

size_t
BarrierGraph::assign_to_basin(const HybEnsModel::StateDescription &source_state,
			      const HybEnsModel::energy_t &source_energy,
			      const HybEnsModel::energy_t &min_transition_energy,
			      size_t min_transE_target_basin_index
			      ) {
    size_t source_basin_index;

    // ------------------------------------------------------------
    // Perform basin assignment of source_state.
    //
    // Either construct a new basin with local minimum source_state
    // or assign source_state to an existing basin
    //
    if ( !gradient_ || min_transition_energy > source_energy ) {
	// no transition state to an energetically lower target
	// state is energetically lower than the source state
	//
	// Consequently, source_state is a new local minimum (among
	// the previously seen states)
	
	source_basin_index = basins_.size();
	
	if (debug_out_) std::cerr << "  New basin "<<source_basin_index<<std::endl;

	// put state into hash
	state_hash_[source_state.encode()] = source_basin_index;
	
	// generate new basin and put into object's basin list
	Basin new_basin(source_basin_index,source_state,source_energy,model_);
	basins_.push_back(new_basin);
	    
    } else {
	// source_state is not a local minimum but belongs to basin
	// min_transE_target_basin_index, which is the basin that is reached
	// with the lowest transition energy.
	    
	// handle case where basin with index source_basin_index
	// was merged before
	source_basin_index = min_transE_target_basin_index;
	
	if (debug_out_) std::cerr << "  Assign to basin "<<source_basin_index<<std::endl;
	    
	// assign basin index source_basin_index to source_state and register
	// source_state as new member of the basin
	state_hash_[source_state.encode()] = source_basin_index;
	basins_[source_basin_index].add_state(source_energy,model_);
    }
    // end assign to basin

    return source_basin_index;
}

void 
BarrierGraph::register_transitions(size_t source_basin_index,
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
	    add_transition(*it,model_);
	    add_transition(it->reverse(),model_);
	    
	}
    } // end iterate trans
}

void
BarrierGraph::process_state(const HybEnsModel::StateDescription &source_state, 
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
		     min_transition_energy,       // <-- the last three params are output parameters!
		     min_transE_neighbor_basin_index,
		     transitions
		     );
    }
   
    if (debug_out_) std::cerr << "  " << transitions.size() << " transitions, "
			     << moves_counter << " moves" <<std::endl;
    
    // handle case "minimum neighbor state is not in hash"
    if (GRAD_WALK_MISSING
	&&
	(min_transition_energy < source_energy) 
	&&
	(min_transE_neighbor_basin_index == std::numeric_limits<size_t>::max())) {
	
	// std::cerr 
	//     << "min outside "
	//     << source_state << " " << source_energy 
	//     << " -> "
	//     << min_transition_energy << std::endl; 
	
	// TODO this is not working properly without gradient basins 
	assert(gradient_/* case !gradient not handled properly */);
	
	if (debug_out_) 
	    std::cerr << source_state << " walks outside maxsitesize boundaries." << std::endl;	

	// do gradient walk from source_state
	gradient_walk(source_state, source_energy);
	

	//FIX LATER
	min_transition_energy = std::numeric_limits<double>::infinity();
    }

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
		std::cerr << "  State "<<to_dotbracket(source_state)
			  <<" is self-symmetric."
			  <<std::endl;
	    }
	} else {
	    // Test whether the symmetric state of source_state has been seen before.
	    // In this case, the two symmetric states are identified.
	    state_hash_t::const_iterator symmetric_state_it=state_hash_.find(symmetric_state.encode());
	    if (symmetric_state_it != state_hash_.end()) {
		if (debug_out_) {
		    std::cerr << "  State "<<to_dotbracket(source_state)
			      <<" is equivalent to previously found "
			      <<to_dotbracket(symmetric_state)
			      <<" of basin "<<symmetric_state_it->second<<"."<<std::endl;
		}
		
		basin_index = symmetric_state_it->second;
		if (debug_out_) std::cerr << "  Assign to basin "<<basin_index<<std::endl;
		state_hash_[source_state.encode()] = basin_index;
		basins_[basin_index].add_state(source_energy,model_);
	    }
	}
    } // end special homodimer symmetry handling
    
    
    // assign basin
    if ( basin_index == std::numeric_limits<size_t>::max()) { //unless assigned by symm handling
	// assign source state to either new basin or the basin of its deepest neighor; the latter, if barrier<=0
	basin_index 
	    = assign_to_basin(source_state,
			      source_energy,
			      min_transition_energy,
			      min_transE_neighbor_basin_index);
    }
    
    // register the transitions from the source basin to other states in the hash
    register_transitions(basin_index,transitions);

}

size_t BarrierGraph::gradient_walk(const HybEnsModel::StateDescription &source_state,
				   double source_energy) {
    
    if (debug_out_) std::cerr << "Walk along from "<<source_state<<" "<<source_energy<<std::endl;
    
    HybEnsModel::StateDescription::code_t source_state_code = source_state.encode();
    state_hash_t::const_iterator it = state_hash_.find(source_state_code);

    if (state_hash_.end() != it) {
	// hey, surprise! the state is already hashed

	if (debug_out_) std::cerr << "... hits the hash at "<<it->second<<"."<<std::endl;
	
	// do nothing but return the basin index of the min tE neighbor
	return it->second;
    }
    
    // minimum transition state energy from source state source_state to a neighbor
    double min_transE = std::numeric_limits<double>::infinity();
    
    // index of the basin of the neighbor state with minimum transition energy
    HybEnsModel::StateDescription min_transE_state;
    
    size_t moves_counter=0;

    // Enumerate neighbors of source_state to get minimum transition
    // energy neighbor
    HybEnsModel::MoveIterator mi(source_state,model_,consider_double_sites_);
    for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	moves_counter++;
	
        HybEnsModel::MoveIterator mi(source_state,model_,consider_double_sites_);
	for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	    HybEnsModel::StateDescription neigh_state = source_state;
	    move->apply(neigh_state);

	    if (tE<min_transE && !isinf(model_.energy(neigh_state))) {
		min_transE = tE;

		min_transE_state = neigh_state;
	    }
	}
    }

    if (debug_out_)  std::cerr << "... min tE: "<<min_transE-source_energy<<std::endl;
    
    if (moves_counter==0) return std::numeric_limits<size_t>::max();
    
    size_t basin_index;
    if (min_transE >= source_energy) {
	// the source state is a local minimum, create basin and terminate
	
	basin_index = basins_.size();
	    
	if (debug_out_) std::cerr << "... new basin "<<basin_index<<"."<<std::endl;
 	
	// generate new basin and put into object's basin list
	Basin new_basin(basin_index,source_state,source_energy,model_);
	basins_.push_back(new_basin);
    } else {
	// there is a gradient neighbor, walk along to find the final basin index
	
	basin_index = gradient_walk(min_transE_state,min_transE);
	
	basins_[basin_index].add_state(source_energy,model_);
    }

    // put state into hash
    state_hash_[source_state_code] = basin_index;
    
    return basin_index;
}

size_t
BarrierGraph::num_basins() const {
    size_t n=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    n++;
	}
    }
    return n;
}


size_t
BarrierGraph::num_transitions() const {
    size_t n=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    const transitions_map_t::const_iterator &it = transitions_.find(i);
	    if(it != transitions_.end()) {
		n += it->second.size();
	    }
	}
    }
    return n;
}

void
BarrierGraph::filter_basin_transitions(const Basin &x0, double min_rate) {
    transitions_map_row_t &trs_x0=transitions_[x0.idx()];
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ) {
	
	//if (basins_[it->first].merged()) continue;
	assert(!basins_[it->first].merged());
	    
	if (it->second.get_Z() / std::min(x0.get_Z(), basins_[it->first].get_Z())
	    < min_rate) {
	    
	    //remove transition
	    //std::cerr << "Remove transition between "<<x0.idx()<<" and "<<it->first<<std::endl; 
	    transitions_[it->first].erase(x0.idx());
	    it=trs_x0.erase(it);
	    
	} else {
	    ++it;
	}
    }
    //std::cerr << "DONE"<<std::endl;
}

void
BarrierGraph::filter_transitions(double min_rate) {
    // run through basins
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[i];
	
	filter_basin_transitions(x0,min_rate);
    }    
}


void
BarrierGraph::merge_basin(Basin &x0) {

    // compute total outflow
    double total_out = outflow_pf(x0);
    
    transitions_map_row_t &trs_x0=transitions_[x0.idx()];
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ++it) {
	
	Basin &y = basins_[it->first];
	if (y.merged()) continue;
	
	if (x0.idx()==y.idx()) continue;
		
	// 1) distribute the basin's partition function to its neighbors
	double Z_yx0 =  it->second.get_Z();
	double fraction = Z_yx0/total_out;
		
	y.merge_in(x0,fraction);
		
	if (debug_out_) {
	    std::cerr << "Transfer a fraction of " << fraction << " of " << x0.idx()
		      << "'s partfunc to " << y.idx()<<std::endl;
	}
		
	// 2) distribute the partition function of the transitions to this basin to its neighbors
		
	for (transitions_map_row_t::iterator it2=trs_x0.begin();
	     trs_x0.end()!=it2; ++it2) {
		    
	    Basin &x = basins_[it2->first];
	    if (x.merged()) continue;
				    
	    if (x.idx()==y.idx()) continue;
	    if (x0.idx()==x.idx()) continue;
		    
	    // update the transition from x to y (via x0)
		    
	    double Z_xx0 = transitions_[x.idx()][x0.idx()].get_Z();
		    
	    if (debug_out_) {
		std::cerr << "Transfer a fraction of " << fraction << " of transition "
			  << x.idx() <<"->"<< x0.idx() << " to " << x.idx() 
			  << "->"<< y.idx()<<std::endl;
	    }
		    
	    transitions_[x.idx()][y.idx()].update(Z_xx0 * fraction);		    
	} // end for it2 (over neighbors of x0)
		
    } // end for it (over neighbors of x0)
	    
    // finally, mark as merged
    if (debug_out_) {
	std::cerr << "Mark "<<x0.idx()<<" as merged."<<std::endl;
    }
    x0.mark_merged();
	    
    // and release all transitions from and to the merged basin
    for (transitions_map_row_t::iterator it=trs_x0.begin();
	 trs_x0.end()!=it; ++it) {
	transitions_[it->first].erase(x0.idx());
    }
    transitions_.erase(x0.idx());
    
}


bool
BarrierGraph::is_to_be_merged(const Basin &x0,double max_outflow,double min_p_equ) const {
    double total_Z = compute_Z();
    
    // compute total outflow
    double total_out = outflow_pf(x0);
    
    // probability in equilibrium
    double p_equ = x0.get_Z() / total_Z;
    
    bool to_be_merged = 
	(p_equ < min_p_equ) 
	|| (total_out/x0.get_Z() > max_outflow)
	;
    
    if (debug_out_) {
	if (to_be_merged) {
	    std::cerr << "Select basin "<<x0.idx()<<" for merge because ";
	    if (total_out/x0.get_Z()) {
		std::cerr << "total outflow is "<<total_out/x0.get_Z();
	    }
	    if (p_equ < min_p_equ) {
		std::cerr << "equilibrium probability is "<<p_equ;
	    }
	    std::cerr<<std::endl;
	}
    }
    
    if (special_open_state_ && x0.get_local_minimum().size()==0) {
	if (verbose_) {
	    std::cerr << "Suppress merge of open state basin."<<std::endl;
	}
	to_be_merged=false;
    }
    
    return to_be_merged;
}

void
BarrierGraph::merge_basins_by_outflow(double max_outflow,double min_p_equ, double min_rate) {
    // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins_.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
        
    // run through sorted basins and merge
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions_[x0.idx()];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.idx())) {
	    trs_x0.erase(x0.idx());
	}
	
	// remove transitions wiht rates lower than min_rate
	filter_basin_transitions(x0,min_rate);
	            
        if (is_to_be_merged(x0,max_outflow,min_p_equ)) {
	    
	    if (debug_out_) {
		std::cerr << "Dissolve basin "<< x0.idx() << " with outflow "<< outflow_pf(x0)/x0.get_Z() << std::endl;
	    }
	    
	    merge_basin(x0);
	
	}
    }
}


void
BarrierGraph::reduce_basin_set(const std::set<size_t> &to_keep, double min_rate) {
    // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins_.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
        
    // run through sorted basins and merge
    for (size_t i=0; i<basins_.size(); ++i) {
	
	Basin &x0 = basins_[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions_[x0.idx()];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.idx())) {
	    trs_x0.erase(x0.idx());
	}
	
	// remove transitions wiht rates lower than min_rate
	filter_basin_transitions(x0,min_rate);
	
        if (to_keep.count(x0.idx()) == 0) {
	    
	    if (verbose_) {
		std::cerr << "Dissolve basin "<< x0.idx() << std::endl;
	    }
	    
	    merge_basin(x0);
	
	}
    }
}


BarrierGraph::BarrierGraph(const std::string &seqA, 
			   const std::string &seqB,
			   std::istream &in,
			   bool binary,
			   bool special_open_state,
			   bool consider_double_sites,
			   bool gradient,
			   size_t maxsitesize,
			   bool verbose,
			   bool debug_out)
    :
    model_(seqA,seqB,std::max(seqA.length(),seqB.length()),true), // note: not restricted by maxsitesize
    special_open_state_(special_open_state),
    consider_double_sites_(consider_double_sites),
    gradient_(gradient),
    maxsitesize_(maxsitesize),
    verbose_(verbose),
    debug_out_(debug_out)
{
        
    HybEnsModel::StateDescription source_state;
    double energy;
        
    // counter for states that are read for construction of graph
    size_t state_counter;

    state_counter=0;
    
    // use to check input
    double last_energy=-std::numeric_limits<double>::infinity();
    int line=1;
    
    if (special_open_state_) {
	// add the empty/open state first, to guarantee to create a new basin 

	HybEnsModel::StateDescription open_state;

	if (debug_out_) { 
	    std::cerr << "add open state " << open_state << std::endl;
	}
	
	process_state(open_state, -4.1);
	
	state_counter++;
    }
    

    while (read_state(in,source_state,energy,line,binary)) {
	if (verbose_ && (state_counter%5000==0)) std::cerr << "\r" << state_counter;

	if (energy < last_energy) {
	    std::cerr << "ERROR: input states have to be sorted by increasing energy (at line "
		      <<line<<": "<<energy<<"<"<<last_energy << " ).\n";
	    exit(-1);
	}

	if (debug_out_) {
	    std::cerr << "read " << state_counter << " " << energy << " "
		      << " "  << source_state << std::endl;
	}
	
	if (special_open_state_ && source_state.size()==0) {
	    if (verbose_) {
		std::cerr << std::endl << "Ignore open state in input." <<std::endl;
	    }
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

void
BarrierGraph::print_edges(std::ostream &out,const Basin &b) const {
    out << b.idx()
	// << " (" << b.get_Z() << ") "
	<< "  -> ";
	
    transitions_map_t::const_iterator
	ts = transitions_.find(b.idx());
	
    if (ts == transitions_.end()) return;
	
    for(transitions_map_row_t::const_iterator it=ts->second.begin();
	it != ts->second.end();
	++it) {
	    
	if ((it->first!=b.idx()) && (!basins_[it->first].merged())) {
	    out << " " << it->first
		// << " (" << it->second.get_Z() << ")"
		<< " " << (it->second.get_Z() / b.get_Z());
	}
    }
    out << std::endl;
}

void 
BarrierGraph::print_barrier_graph(std::ostream &out) const {
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    print_edges(out,basins_[i]);
	}
    }	
}

double BarrierGraph::compute_Z() const {
    double total_Z=0.0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    total_Z += basins_[i].get_Z();
	}
    }
    return total_Z;
}


void
BarrierGraph::print_basins(std::ostream &out) const {
    basins_[0].print_header(out);
    out<< " \tmax_out \tensE_out \ttotal_out \tp_equ";
    out<<std::endl;

    double total_Z = compute_Z();

    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    basins_[i].print(out,model_);
	    
	    double max_out=max_outflow(basins_[i]);
	    
	    double ensE_max_out = - model_.RT() * log(max_out);
	    
	    double p_equ = basins_[i].get_Z() / total_Z;
	    
	    out << " \t" << max_out<< " \t" << ensE_max_out<<" \t"<<outflow(basins_[i])
		<<" \t" << p_equ;
	    
	    out<<std::endl;
	} else {
	    // printf("merged to %4lu",basins_[i].merged_to); basins_[i].print(std::cout,model);
	}
    }
}

std::string
BarrierGraph::format_rate_for_treekin(double x)  {
    // using sprintf here is somwehat ugly, but there seem to be no
    // easy workarounds (i.e. without boost)
    const size_t bufsiz=256;
    char buf[bufsiz];
    snprintf(buf,bufsiz,"%10.4g",x);
    buf[bufsiz-1]=0;
    return std::string(buf);
}


std::string
BarrierGraph::to_dotbracket(const HybEnsModel::StateDescription &sd) const {
    size_t n = model_.seqA().length();
    size_t m = model_.seqB().length();
    
    std::string s(n+m+1,'.');
    
    s[n] = '&';
    
    for (size_t k=0; k<sd.size(); k++) {
	s[sd[k].i1-1]='(';
	s[sd[k].j1-1]=')';
	s[n+1+sd[k].i2-1]='(';
	s[n+1+sd[k].j2-1]=')';
    }
    
    return s;
}


void
BarrierGraph::print_barriers(std::ostream &out) const {
    out << "     " << model_.seqA() << "&" << model_.seqB() << std::endl;
    
    size_t count=0;
    for(size_t i=0; i<basins_.size(); ++i) {
	if (basins_[i].merged()) continue;
	count++;
	size_t ow = out.width(4);

	std::ostringstream energy;
	energy.setf(std::ios_base::fixed, std::ios_base::floatfield);
	energy << std::setprecision(2)
	       << std::setw(6)
	       << (- model_.RT() * log(basins_[i].get_Z()));

	out << count << std::setw(ow) 
	    << " " << to_dotbracket(basins_[i].get_local_minimum())
	    << " " << energy.str() << std::endl;
    }
}

void
BarrierGraph::print_treekin_ratesmatrix(std::ostream &out) const
{
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    
	    double total=0.0; // row total rate (i!=j)
	    
	    std::vector<double> row(basins_.size());
	    
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		// "bucket sort"
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    double rate = (it->second.get_Z() / basins_[i].get_Z());
		    row[it->first] = rate;
		    total += rate;
		}
	    }
	    
	    row[i] = -total;

	    for(size_t j=0; j<row.size();++j) {
		if (!basins_[j].merged()) {
		    out << format_rate_for_treekin(row[j]) << " ";
		}
	    }
	    out << "\n";
	}
    }
}


void
BarrierGraph::print_pfs(std::ostream &out,bool binary) const {
    if (binary) {
	// write dimension
	size_t dim=0;
	for(size_t i=0; i<basins_.size(); ++i) {
	    if (!basins_[i].merged()) {
		dim++;
	    }
	}
	out.write(reinterpret_cast<const char *>(&dim),sizeof(dim));
    }
    
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    
	    std::vector<double> row(basins_.size());
	    
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		// "bucket sort"
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    row[it->first] = it->second.get_Z();
		}
	    }
	    
	    // set diagonal
	    row[i] = basins_[i].get_Z();
	    
	    for(size_t j=0; j<row.size();++j) {
		if (!basins_[j].merged()) {
		    if (binary) {
			out.write(reinterpret_cast<const char *>(&row[j]),sizeof(row[j]));
		    } else {
			out << row[j] << " ";
		    }
		}
	    }
	    if (!binary) {out << "\n";}
	}
    }
}

void
BarrierGraph::print_rxns(std::ostream &out,
			 const std::string &nameA,
			 const std::string &nameB,
			 const std::string &nameAB) const {
    const size_t openstate_index=1;
   
    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged()) {
	    const transitions_map_t::const_iterator &ts = transitions_.find(i);
	    if (ts != transitions_.end()) {
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    size_t j=it->first;
		    double rate = it->second.get_Z()/basins_[i].get_Z();
		    
		    std::ostringstream subs1;
		    std::ostringstream subs2;
		    std::ostringstream prods1;
		    std::ostringstream prods2;
		    
		    if (i!=openstate_index) {
			subs1<<nameAB<<i+1;
			subs2<<nameAB<<i+1;
		    } else {
			subs1<<nameA<<"+"<<nameB;
			subs2<<nameA<<" "<<nameB;
		    }
		    
		    if (j!=openstate_index) {
			prods1<<nameAB<<j+1;
			prods2<<nameAB<<j+1;
		    } else {
			prods1<<nameA<<"+"<<nameB;
			prods2<<nameA<<" "<<nameB;
		    }

		    out << "Rxn " << subs1.str() << "->" << prods1.str() <<std::endl
			<< "Subs " << subs2.str() << std::endl
			<< "Prods " << prods2.str() <<std::endl
			<< "Rate " << rate <<std::endl;    
		    
		    out << std::endl;
		}
	    }
	}
    }
}

void
BarrierGraph::print_spcs(std::ostream &out,
			 const std::string &nameA,
			 const std::string &nameB,
			 const std::string &nameAB) const {
    const size_t openstate_index=1;

    if (nameA!=nameB) {
	out << nameA << " " << 1 << std::endl;
	out << nameB << " " << 1 << std::endl;
    }

    for(size_t i=0; i<basins_.size(); ++i) {
	if (!basins_[i].merged() && (i!=openstate_index)) {
	    out << nameAB << (i+1) << " " << 0<< std::endl;
	}
    }
}

void
BarrierGraph::reindex() {
    std::vector<bool> keep(basins_.size());
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	keep[it->idx()] = !it->merged();
    }
    reindex(keep);
}

void
BarrierGraph::keep_single_component(size_t c,const std::vector<size_t> &components) {
    std::vector<bool> keep(basins_.size());
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	keep[it->idx()] = components[it->idx()]==c;
    }
    reindex(keep);
}


void
BarrierGraph::reindex(const std::vector<bool> &keep) {
    
    std::vector<size_t> old2new(basins_.size());
    
    
    // phase 0: remove transitions of merged basins_
    for (std::vector<Basin>::iterator it=basins_.begin();
	 basins_.end()!=it; ++it) {
	if (!keep[it->idx()]) {
	    transitions_.erase(it->idx());
	}
    }

    for(transitions_map_t::iterator it=transitions_.begin();
	it != transitions_.end(); ++it) {
	for (std::vector<Basin>::iterator it2=basins_.begin();
	     basins_.end()!=it2; ++it2) {
	    if (!keep[it2->idx()]) {
		it->second.erase(it2->idx());
	    }
	}
    }

    // phase 1: reindex basins_; copy in place
    {
	size_t i=0;
	for (std::vector<Basin>::iterator it=basins_.begin();
	     basins_.end()!=it; ++it) {
	    if (keep[it->idx()]) {
		old2new[it->idx()]=i;
		if (it->idx()!=i) {
		    it->set_idx(i);
		    basins_[i]=*it;
		}
		i++;
	    }
	}
	basins_.resize(i);
    }
    
    // phase 2; remap transitions
    transitions_map_t new_transitions;

    for(transitions_map_t::const_iterator it=transitions_.begin();
	it != transitions_.end(); ++it) {
	for(transitions_map_row_t::const_iterator it2=it->second.begin();
	    it2 != it->second.end(); ++it2) {
	    new_transitions[old2new[it->first]][old2new[it2->first]] = it2->second;
	}
    }

    
    transitions_.clear();
    transitions_ = new_transitions;

    for (size_t i=0; i<basins_.size(); i++) {
	if (transitions_.find(i) == transitions_.end()) {
	    //if (verbose_) std::cerr << "WARNING: no transitions for state "<<i<<std::endl;
	    transitions_[i] = transitions_map_row_t();
	}
    }
}

double
BarrierGraph::check_rates() const
{
    // check whether transition rate matrix is symmetric 
    // (i.e. whether process is in detailed balance)

    double max_diff=0.0;

    for(size_t i=0; i<basins_.size(); ++i) {
     	std::vector<double> row(basins_.size());
	
	const transitions_map_t::const_iterator &ts = transitions_.find(i);
	if (ts != transitions_.end()) {
	    for(transitions_map_row_t::const_iterator it=ts->second.begin();
		it != ts->second.end(); ++it) {
		size_t j=it->first;
		
		//warning: the following finds can fail if data structures are not symmetric
		double diff = abs(it->second.get_Z() -
				  transitions_.find(j)->second.find(i)->second.get_Z());
		
		diff /= it->second.get_Z();
		
		max_diff = std::max(max_diff,diff);
	    }
	}
    }
    
    return max_diff;
}

std::vector<size_t>
BarrierGraph::connected_components(std::vector<size_t> &components) const
{
    // components[i] is the number of the component of basin i
    // or 0 if the basin is still unassigned
    components.clear();
    components.assign(basins_.size(),0);
    
    std::vector<size_t> component_size;
    
    // current component number
    size_t c=0;
    
    typedef std::pair<size_t,transitions_map_row_t::const_iterator> stack_entry;
    std::stack<stack_entry> stack;
    
    for (size_t i=0; i<basins_.size(); i++) {
	
	if (components[i] > 0) continue;
	c++;

	stack.push(stack_entry(i,transitions_.find(i)->second.begin()));
	components[i]=c;
	component_size.push_back(1);

	while (!stack.empty()) {
	    stack_entry &top = stack.top();
	    
	    while (top.second != transitions_.find(top.first)->second.end()
		   &&
		   components[top.second->first]>0) 
		top.second++;
	    
	    if (top.second == transitions_.find(top.first)->second.end()) {
		stack.pop();
		continue;
	    }
	    
	    size_t j=top.second->first;
	    stack.push(stack_entry(j,transitions_.find(j)->second.begin()));
	    components[j]=c;
	    component_size[c-1]++;
	}
    }
    

    return component_size;
}

void
BarrierGraph::swap_indices(size_t x, size_t y) {
	
    // first swap basins
    std::swap(basins_[x],basins_[y]);
    // (don't forget to swap indices of basins)
    size_t tmp = basins_[x].idx();
    basins_[x].set_idx(basins_[y].idx());
    basins_[y].set_idx(tmp);
    
    // then swap transitions:
    // -- first swap the rows
    std::swap(transitions_[x],transitions_[y]);
    
    // -- then columns
    for (transitions_map_t::iterator it=transitions_.begin(); transitions_.end()!=it; ++it) {
	transitions_map_row_t &row = it->second;
	
	transitions_map_row_t::iterator itx = row.find(x);
	transitions_map_row_t::iterator ity = row.find(y);
	
	if (itx!=row.end() && ity!=row.end()) {
	    std::swap(row[x],row[y]);
	} else if (itx==row.end() && ity!=row.end()) {
	    row[x] = ity->second;
	    row.erase(y);
	} else if (itx!=row.end() && ity==row.end()) {
	    row[y] = itx->second;
	    row.erase(x);
	}
    }
}
