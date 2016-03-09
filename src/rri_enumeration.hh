#ifndef RRI_ENUMERATION_HH
#define RRI_ENUMERATION_HH

#include <iostream>

#include <string>
#include <cmath>

#include <limits>

//#include <cassert>

#include "util.hh"
#include "hybrid_ensemble_model.hh"

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}

/**
 * Enumeration of rrikin states: generate, write and read states
 *
 * @todo Implement storing of states in memory; read and write
 * methods; in memory sorting
 */
class RRIEnumeration {
    const HybEnsModel &model_;
    double max_hyb_energy_;
    double max_total_energy_;
    bool binary_;

    
public:
    
    /** 
     * @brief Construct with all parameters
     * @param out output stream
     * @param model the hybrid ensemble model
     * @param max_hyb_energy maximum energy of a hybridisation site
     * @param max_total_energy maximum total energy of a state
     */
    RRIEnumeration(const HybEnsModel &model,
		     double max_hyb_energy,
		     double max_total_energy,
		     bool binary)
	:
	model_(model),
	max_hyb_energy_(max_hyb_energy),
	max_total_energy_(max_total_energy),
	binary_(binary)
    {}
    

    /**
     * @brief Write state
     * @param ostream output stream
     * @param energy energy of state
     * @param state description of state
     */
    static
    void
    write_state(std::ostream &out,
		double energy, 
		const HybEnsModel::StateDescription &state,
		bool binary);
    
    /** 
     * @brief Read state
     *
     * Read a state description line from input stream
     * 
     * @param in input stream 
     * @param[out] state state description 
     * @param[out] energy state energy
     * @param lineno line number (for error output)
     * 
     * @return whether valid state description line could be read
     * @note Exit on invalid input!
     */
    static
    bool
    read_state(std::istream &in, HybEnsModel::StateDescription &state, 
	       double &energy, size_t lineno,bool binary);

    /**
     * @brief Enumerate single site states
     * @param in_memory whether to store in memory
     * @param ostream output stream
     * @return number of sites
     */
    size_t
    enumerate_single_sites(bool in_memory, std::ostream &out=std::cout) const;
    
    /**
     * @brief Enumerate double site states
     * @param in_memory whether to store in memory
     * @param ostream output stream
     * @return number of sites
     */
    size_t
    enumerate_double_sites(bool in_memory, std::ostream &out=std::cout) const;

    // /** @brief sort states by energy (best first)
    //  */
    // void
    // sort_by_energy();

    // /**
    //    @brief Read states
    // */
    // std::istream &
    // read_states(std::istream &in);
    
    // /**
    //    @brief Write states
    // */
    // std::ostream &
    // write_states(std::ostream &out);
    
};

#endif // RRI_ENUMERATION_HH
