#include "basin_info.hh"

BasinInfo::BasinInfo(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy)
    : min_state_(state), min_energy_(energy) {
    
    states_.push_back(state);
}


void
BasinInfo::add_state(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy) {
    states_.push_back(state);

    if (energy < min_energy_) {
	min_state_  = state;
	min_energy_ = energy;
    }
}

void
BasinInfo::interaction_pair_pfs(const HybEnsModel &model, 
				pf_matrix_t &pfs) const {
    pfs.clear();
    
    // partition function of this basin
    // double basin_pf=0.0;
    
    for ( auto &code: states_ ) {
	HybEnsModel::StateDescription sd(code);
	// partition function of the current state
	double state_pf = model.partition_function(sd);
	// basin_pf += state_pf;
	
	// iterate over all possible interactions (k1,k2),
	// calculate probability Pr[(k1,k2)|Basin]; 
	// store in pfs
	
	for (auto &isite: sd) {
	    
	    // Get slice of pf matrix for common left ends i1,i2.
	    // This allows us to make use of the sparsity of this
	    // matrix.
	    const auto &pf_matrix_slice =
		model.hybrid_pf().hybrid_pfs_common_left_ends(isite.i1,isite.i2);
	    
	    for ( auto &x : pf_matrix_slice ) {
		size_t k1 = x.first.first;
		size_t k2 = x.first.second;
		pfs(k1,k2) += state_pf * model.interaction_probability(k1,k2,sd);
	    }
	}
    }
    
    // // this would transform the partition functions to probabilities:
    // // normalize all entries by basin_pf
    // for ( auto &x: pfs ) {
    // 	pfs(x.first.first,x.first.second) = x.second / basin_pf; 
    // }
}

std::ostream &operator << (std::ostream &out, const BasinInfo &bi) {
    
    HybEnsModel::StateDescription minstate(bi.min_state());
    out << minstate << "\t" << bi.min_energy() << "\t";
    
    for (BasinInfo::state_vec_t::const_iterator it=bi.states_.begin(); bi.states_.end()!=it; ++it) {
	if (*it != bi.min_state_) {
	    HybEnsModel::StateDescription state(*it);
	    out << state << " ";
	}
    }
    
    return out;
}
