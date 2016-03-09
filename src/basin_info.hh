#ifndef BASIN_INFO_HH
#define BASIN_INFO_HH

#include <iostream>
#include "hybrid_ensemble_model.hh"

/**
   Additional information about a basin
   
   Stores basin information going beyond the data in class Basin. The
   information here is not needed for the construction of the basin
   graph and computation of transition rates, however here we keep
   track of the states in the basin.
*/
class BasinInfo {
    HybEnsModel::StateDescription::code_t min_state_; //!< minimum state in basin
    HybEnsModel::energy_t min_energy_;                //!< energy of minimum state

    typedef std::vector<HybEnsModel::StateDescription::code_t> state_vec_t; 
    
    // vector of (codes of) all states in the basin
    state_vec_t states_;
    
    typedef state_vec_t::const_iterator const_iterator;
    
public:

    /**
     * @brief construct with first state of basin 
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
    add_state(const HybEnsModel::StateDescription::code_t &state, 
	      HybEnsModel::energy_t energy);
    
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
    HybEnsModel::energy_t
    min_energy() const {
	return min_energy_;
    }

    typedef LocARNA::SparseMatrix<double> pf_matrix_t; 

    /**
     * @brief Compute partition functions for interaction pairs in basin
     * @param model the hybrid ensemble model
     * @param[out] pfs matrix of partition functions Z^bp[ (k1,k2) | Basin ]
     */
    void
    interaction_pair_pfs(const HybEnsModel &model,
			 pf_matrix_t &pfs) const;
    
    friend std::ostream &operator << (std::ostream &out, const BasinInfo &bi);
    
    /**
     * @brief begin of states vector
     *
     * @return constant iterator pointing to begin of vector of states
     * (i.e. codes of states) in the basin
     */
    const_iterator
    begin() const {return states_.begin();}
    
    /**
     * @brief end of states vector
     *
     * @return constant iterator pointing to end of vector of states
     * in the basin
     * @see begin()
     */
    const_iterator
    end() const {return states_.end();}

};


#endif // BASIN_INFO_HH
