#include "basin_transition.hh"
#include "barrier_graph.hh"

#include <stack>
#include <iostream>
#include <stdio.h>
#include <iomanip>

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

    if (not state.is_valid(model)) {
	std::cerr << "ERROR: read state "<<state<<" at line "<<lineno<<" is not valid in model."<<std::endl;
	exit(-1);
    }
        
    return true;
}


void
BarrierGraph::add_transition( const transition_t &tr, const HybEnsModel &model) {
    // add transition to partition function for the transition
    // between the source and target basin
    if (transitions[tr.source_basin_index].find(tr.target_basin_index)
	== transitions[tr.source_basin_index].end()) {
	transitions[tr.source_basin_index][tr.target_basin_index] = BasinTransition(tr.barrier_energy(),model);
    } else {
	transitions[tr.source_basin_index][tr.target_basin_index].update(tr.barrier_energy(),model);
    }
}
    

void
BarrierGraph::process_state(const HybEnsModel::StateDescription &source_state, double source_energy) {

    
    if (!consider_double_sites && source_state.num_sites()==2) {
	if (debug_out) {
	    std::cerr << "Ignore double site state " << source_state <<std::endl;
	}
	return;
    }
    
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
    // THIS IMPLIES that each transition between x and y is registered only once.
    // 
    // Keep track of neighbor state with smallest transition energy. 
    //
    //
    HybEnsModel::MoveIterator mi(source_state,model);
    for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    
	//std::cout << " move "; move->print(std::cout); std::cout<<std::endl;
	    
	moves_counter++;
	    
	HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	HybEnsModel::StateDescription neigh_state=source_state;
	move->apply(neigh_state);

	if (!consider_double_sites && neigh_state.num_sites()==2) {
	    continue;
	}
	
	
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
	    
	HybEnsModel::StateDescription::code_t neigh_code; // string for holding code
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
    } // end iterate moves/neighbors

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
    // end assign to basin


    // ----------------------------------------
    // Register all transitions from source_basin_index to other basins.
    //
    for(std::vector<transition_t>::iterator it=trans.begin(); 
	trans.end()!=it; ++it) {
	if (it->target_basin_index != source_basin_index) {
		
	    it->source_basin_index = source_basin_index;
		
	    if (debug_out)
		std::cerr << "add transition "<<source_basin_index<<" <-> "
			  << it->target_basin_index << " "
		    //<< it->transition_energy << " "
			  << it->barrier_energy() << " "
			  << std::endl;
	

	    // note: each transition between source and target is
	    // registered only once, namely when the source index is
	    // larger, therefore we add transitions in both directions
	    // (finally, this produces a symmetric matrix <transitions>
	    add_transition(*it,model);
	    add_transition(it->reverse(),model);
	    
	}
    } // end iterate trans	
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
	    if(it != transitions.end()) {
		n += it->second.size();
	    }
	}
    }
    return n;
}

void
BarrierGraph::merge_basins(double max_outflow, double min_p_equ, double min_rate) {
	
    // // sort basins increasing by their partition function
    std::vector<size_t> sorted_basin_idxs;

    for (size_t i=0; i<basins.size(); ++i) sorted_basin_idxs.push_back(i);
    sort(sorted_basin_idxs.begin(),sorted_basin_idxs.end(),compBasinIdxs(*this));
    
    // decide on merge non-recursively
    
    double total_Z = compute_Z();
    
    // first determine for each basin, whether it should be merged
    std::vector<bool> to_be_merged(basins.size());
    for (size_t i=0; i<basins.size(); ++i) {
	Basin &x0 = basins[i];
	// compute total outflow
        double total_out = outflow_pf(x0);
	
	// probability in equilibrium
	double p_equ = basins[i].get_Z() / total_Z;
	
       	to_be_merged[i] = 
	    (total_out/x0.get_Z() > max_outflow)
	    || (p_equ < min_p_equ);
        
	if (x0.get_local_minimum().num_sites()==0) {
	    if (verbose) {
		std::cerr << "WARNING: Suppress merge of open state basin."<<std::endl;
	    }
	    to_be_merged[i]=false;
	}
	
	if (debug_out) {
	    std::cerr << i << " " << sorted_basin_idxs[i] << " " << x0.get_Z() << " " << total_out << " r=" << (total_out/x0.get_Z());
	    
	    std::cerr <<"  \t";
	}
    }
    
    // run through sorted basins and merge
    for (size_t i=0; i<basins.size(); ++i) {
	
	Basin &x0 = basins[sorted_basin_idxs[i]];
	
	//if (x0.merged()) continue;
	assert(!x0.merged());
	
	
	transitions_map_row_t &trs_x0=transitions[x0.idx()];
	
	// remove transitions from state x0 to itself
	if (trs_x0.end()!=trs_x0.find(x0.idx())) {
	    trs_x0.erase(x0.idx());
	}

	// remove transitions wiht rates lower than min_rate
	for (transitions_map_row_t::iterator it=trs_x0.begin();
	     trs_x0.end()!=it; ) {
	    
	    //if (basins[it->first].merged()) continue;
	    assert(!basins[it->first].merged());
	    
	    if (it->second.get_Z() / std::min(x0.get_Z(), basins[it->first].get_Z())
		< min_rate) {
		
		//remove transition
		//std::cerr << "Remove transition between "<<x0.idx()<<" and "<<it->first<<std::endl; 
		transitions[it->first].erase(x0.idx());
		it=trs_x0.erase(it);
		
	    } else {
		++it;
	    }
        }
	//std::cerr << "DONE"<<std::endl;
	
	
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

        if (to_be_merged[x0.idx()]) {

	    // compute total outflow
	    double total_out = outflow_pf(x0);
	    
	    if (debug_out) {
		std::cerr << "Dissolve basin "<< x0.idx() << std::endl;
	    }

	    for (transitions_map_row_t::iterator it=trs_x0.begin();
		 trs_x0.end()!=it; ++it) {
		
		Basin &y = basins[it->first];
		if (y.merged()) continue;
		
		if (x0.idx()==y.idx()) continue;
		
		// 1) distribute the basin's partition function to its neighbors
		double Z_yx0 =  it->second.get_Z();
		double fraction = Z_yx0/total_out;
		
		y.merge_in(x0,fraction);
		
		if (debug_out) {
		    std::cerr << "Transfer a fraction of " << fraction << " of " << x0.idx()
			      << "'s partfunc to " << y.idx()<<std::endl;
		}
		
		// 2) distribute the partition function of the transitions to this basin to its neighbors

		for (transitions_map_row_t::iterator it2=trs_x0.begin();
		     trs_x0.end()!=it2; ++it2) {
		    
		    Basin &x = basins[it2->first];
		    if (x.merged()) continue;
				    
		    if (x.idx()==y.idx()) continue;
		    if (x0.idx()==x.idx()) continue;
		    
		    // update the transition from x to y (via x0)
		    
		    double Z_xx0 = transitions[x.idx()][x0.idx()].get_Z();
		    
		    if (debug_out) {
			std::cerr << "Transfer a fraction of " << fraction << " of transition "
				  << x.idx() <<"->"<< x0.idx() << " to " << x.idx() 
				  << "->"<< y.idx()<<std::endl;
		    }
		    
		    transitions[x.idx()][y.idx()].update(Z_xx0 * fraction);		    
		} // end for it2 (over neighbors of x0)
		
	    } // end for it (over neighbors of x0)
	    
	    // finally, mark as merged
	    if (debug_out) {
		std::cerr << "Mark "<<x0.idx()<<" as merged."<<std::endl;
	    }
	    x0.mark_merged();
	    
	    // and release all transitions from and to the merged basin
	    for (transitions_map_row_t::iterator it=trs_x0.begin();
		 trs_x0.end()!=it; ++it) {
		transitions[it->first].erase(x0.idx());
	    }
	    transitions.erase(x0.idx());
	    
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


BarrierGraph::BarrierGraph(const std::string &seqA, const std::string &seqB,
			   bool binary,
			   bool consider_double_sites_)
    :
    debug_out(false),
    verbose(false),
    model(seqA,seqB),
    consider_double_sites(consider_double_sites_)
{
        
    HybEnsModel::StateDescription orig_state;
    double energy;
        
    // counter for states that are read for construction of graph
    size_t state_counter;

    state_counter=0;
    
    // use to check input
    double last_energy=-std::numeric_limits<double>::infinity();
    int line=1;

    while (read_state(std::cin,orig_state,energy,line,binary)) {
	if (verbose && (state_counter%5000==0)) std::cerr << "\r" << state_counter;

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
    std::cerr << "\r";
}

void
BarrierGraph::print_edges(std::ostream &out,const Basin &b) const {
    out << b.idx()
	// << " (" << b.get_Z() << ") "
	<< "  -> ";
	
    transitions_map_t::const_iterator
	ts = transitions.find(b.idx());
	
    if (ts == transitions.end()) return;
	
    for(transitions_map_row_t::const_iterator it=ts->second.begin();
	it != ts->second.end();
	++it) {
	    
	if ((it->first!=b.idx()) && (!basins[it->first].merged())) {
	    out << " " << it->first
		// << " (" << it->second.get_Z() << ")"
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

double BarrierGraph::compute_Z() const {
    double total_Z=0.0;
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    total_Z += basins[i].get_Z();
	}
    }
    return total_Z;
}


void
BarrierGraph::print_basins(std::ostream &out) const {
    basins[0].print_header(out);
    out<< " \tmax_out \tensE_out \ttotal_out \tp_equ";
    out<<std::endl;

    double total_Z = compute_Z();

    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    basins[i].print(out,model);
	    
	    double max_out=max_outflow(basins[i]);
	    
	    double ensE_max_out = - model.RT() * log(max_out);
	    
	    double p_equ = basins[i].get_Z() / total_Z;
	    
	    out << " \t" << max_out<< " \t" << ensE_max_out<<" \t"<<outflow(basins[i])
		<<" \t" << p_equ;
	    
	    out<<std::endl;
	} else {
	    // printf("merged to %4lu",basins[i].merged_to); basins[i].print(std::cout,model);
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
    size_t n = model.seqA().size();
    size_t m = model.seqB().size();
    
    std::string s(n+m+1,'.');
    
    s[n] = '&';
    
    for (size_t k=0; k<sd.num_sites(); k++) {
	s[sd[k].i1-1]='(';
	s[sd[k].j1-1]=')';
	s[n+1+sd[k].i2-1]='(';
	s[n+1+sd[k].j2-1]=')';
    }
    
    return s;
}


void
BarrierGraph::print_barriers(std::ostream &out) const {
    out << "     " << model.seqA() << "&" << model.seqB() << std::endl;
    
    size_t count=0;
    for(size_t i=0; i<basins.size(); ++i) {
	if (basins[i].merged()) continue;
	count++;
	size_t ow = out.width(4);

	std::ostringstream energy;
	energy.setf(std::ios_base::fixed, std::ios_base::floatfield);
	energy << std::setprecision(2)
	       << std::setw(6)
	       << (- model.RT() * log(basins[i].get_Z()));

	out << count << std::setw(ow) 
	    << " " << to_dotbracket(basins[i].get_local_minimum())
	    << " " << energy.str() << std::endl;
    }
}

void
BarrierGraph::print_treekin_ratesmatrix(std::ostream &out) const
{
    for(size_t i=0; i<basins.size(); ++i) {
	if (!basins[i].merged()) {
	    
	    double total=0.0; // row total rate (i!=j)
	    
	    std::vector<double> row(basins.size());
	    
	    const transitions_map_t::const_iterator &ts = transitions.find(i);
	    if (ts != transitions.end()) {
		// "bucket sort"  
		for(transitions_map_row_t::const_iterator it=ts->second.begin();
		    it != ts->second.end(); ++it) {
		    double rate = (it->second.get_Z() / basins[i].get_Z());
		    row[it->first] = rate;
		    total += rate;
		}
	    }
	    
	    row[i] = -total;

	    for(size_t j=0; j<row.size();++j) {
		if (!basins[j].merged()) {
		    out << format_rate_for_treekin(row[j]) << " ";
		}
	    }
	    out << "\n";
	}
    }
}


void
BarrierGraph::reindex() {
    std::vector<bool> keep(basins.size());
    for (std::vector<Basin>::iterator it=basins.begin();
	 basins.end()!=it; ++it) {
	keep[it->idx()] = !it->merged();
    }
    reindex(keep);
}

void
BarrierGraph::keep_single_component(size_t c,const std::vector<size_t> &components) {
    std::vector<bool> keep(basins.size());
    for (std::vector<Basin>::iterator it=basins.begin();
	 basins.end()!=it; ++it) {
	keep[it->idx()] = components[it->idx()]==c;
    }
    reindex(keep);
}


void
BarrierGraph::reindex(const std::vector<bool> &keep) {
    
    std::vector<size_t> old2new(basins.size());
    
    
    // phase 0: remove transitions of merged basins
    for (std::vector<Basin>::iterator it=basins.begin();
	 basins.end()!=it; ++it) {
	if (!keep[it->idx()]) {
	    transitions.erase(it->idx());
	}
    }

    for(transitions_map_t::iterator it=transitions.begin();
	it != transitions.end(); ++it) {
	for (std::vector<Basin>::iterator it2=basins.begin();
	     basins.end()!=it2; ++it2) {
	    if (!keep[it2->idx()]) {
		it->second.erase(it2->idx());
	    }
	}
    }

    // phase 1: reindex basins; copy in place
    {
	size_t i=0;
	for (std::vector<Basin>::iterator it=basins.begin();
	     basins.end()!=it; ++it) {
	    if (keep[it->idx()]) {
		old2new[it->idx()]=i;
		if (it->idx()!=i) {
		    it->set_index(i);
		    basins[i]=*it;
		}
		i++;
	    }
	}
	basins.resize(i);
    }
    
    // phase 2; remap transitions
    transitions_map_t new_transitions;

    for(transitions_map_t::const_iterator it=transitions.begin();
	it != transitions.end(); ++it) {
	for(transitions_map_row_t::const_iterator it2=it->second.begin();
	    it2 != it->second.end(); ++it2) {
	    new_transitions[old2new[it->first]][old2new[it2->first]] = it2->second;
	}
    }

    
    transitions.clear();
    transitions = new_transitions;

    for (size_t i=0; i<basins.size(); i++) {
	if (transitions.find(i) == transitions.end()) {
	    //if (verbose) std::cerr << "WARNING: no transitions for state "<<i<<std::endl;
	    transitions[i] = transitions_map_row_t();
	}
    }
}

void
BarrierGraph::check_rates() const
{
    std::vector<std::vector<double> > mat(basins.size());
    
    for(size_t i=0; i<basins.size(); ++i) {
	std::vector<double> row(basins.size());
	
	double total=0.0; // row total rate (i!=j)
	
	const transitions_map_t::const_iterator &ts = transitions.find(i);
	if (ts != transitions.end()) {
	    for(transitions_map_row_t::const_iterator it=ts->second.begin();
		it != ts->second.end(); ++it) {
		double rate = (it->second.get_Z() / basins[i].get_Z());
		row[it->first] = rate;
		total += rate;
	    }
	}
	
	row[i] = 1.0-total;
	
    	mat[i] = row;
    }

    double max_diff=0.0;
    for(size_t i=0; i<basins.size(); ++i) {
	for(size_t j=0; j<basins.size(); ++j) {
	    double diff = fabs(basins[i].get_Z()*mat[i][j] - basins[j].get_Z()*mat[j][i]);
	    
	    max_diff = std::max(max_diff,diff);
	}
    }
    std::cerr << "Maximal deviation from detailed balance: "<<max_diff<<std::endl;
}

std::vector<size_t>
BarrierGraph::connected_components(std::vector<size_t> &components) const
{
    // components[i] is the number of the component of basin i
    // or 0 if the basin is still unassigned
    components.clear();
    components.assign(basins.size(),0);
    
    std::vector<size_t> component_size;
    
    // current component number
    size_t c=0;
    
    typedef std::pair<size_t,transitions_map_row_t::const_iterator> stack_entry;
    std::stack<stack_entry> stack;
    
    for (size_t i=0; i<basins.size(); i++) {
	
	if (components[i] > 0) continue;
	c++;

	stack.push(stack_entry(i,transitions.find(i)->second.begin()));
	components[i]=c;
	component_size.push_back(1);

	while (!stack.empty()) {
	    stack_entry &top = stack.top();
	    
	    while (top.second != transitions.find(top.first)->second.end()
		   &&
		   components[top.second->first]>0) 
		top.second++;
	    
	    if (top.second == transitions.find(top.first)->second.end()) {
		stack.pop();
		continue;
	    }
	    
	    size_t j=top.second->first;
	    stack.push(stack_entry(j,transitions.find(j)->second.begin()));
	    components[j]=c;
	    component_size[c-1]++;
	}
    }
    

    return component_size;
}

