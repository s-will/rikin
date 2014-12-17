#ifndef BASIN_INFO_HH
#define BASIN_INFO_HH

#include <iostream>
#include "hybrid_ensemble_model.hh"

//! Additional information about a basin
class BasinInfo {
    HybEnsModel::StateDescription::code_t min_state_; //!< minimum state in basin
    HybEnsModel::energy_t min_energy_;                //!< energy of minimum state

    typedef std::vector<HybEnsModel::StateDescription::code_t> state_vec_t; 
    
    state_vec_t states_;
    
public:

    /**
     * @brief constructor 
     * 
     * @param state  state
     * @param energy state's energy
     */    
    BasinInfo(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy);

    /** 
     * @brief update info after adding state to basin 
     * 
     * @param state  added state
     * @param energy state's energy
     */    
    void
    update(const HybEnsModel::StateDescription::code_t &state, HybEnsModel::energy_t energy);
    
    /** 
     * @brief access minimum state
     * 
     * @return minimum state 
     */
    const HybEnsModel::StateDescription::code_t &    
    min_state() const {
	return min_state_;
    }
    
    /** 
     * @brief access minimum energy
     * 
     * @return minimum energy 
     */
    HybEnsModel::energy_t min_energy() const {
	return min_energy_;
    }
    
    friend std::ostream &operator << (std::ostream &out, const BasinInfo &bi);

};


#endif // BASIN_INFO_HH
