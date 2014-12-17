#include "basin_info.hh"

BasinInfo::BasinInfo(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy)
    : min_state_(state), min_energy_(energy) {
    
    states_.push_back(state);
}


void
BasinInfo::update(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy) {
    states_.push_back(state);

    if (energy < min_energy_) {
	min_state_  = state;
	min_energy_ = energy;
    }
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
