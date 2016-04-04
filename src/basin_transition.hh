#ifndef BASIN_TRANSITION
#define BASIN_TRANSITION

#include <tr1/unordered_map>
#include <iostream>
//#include <limits>

/**
 * @brief A transition between two basins
 *
 * @note Use this class to store information (barrier, Z) for
 * transitions between two basins in the map of transitions stored in
 * each basin. In contrast, Basin::transition_t stores transitions
 * between (micro)states.
 */
class BasinTransition {
    
    double Z_;       //! partition function of macro transition state

public:
    BasinTransition():
	Z_(0) {}
    
    /** 
     * Construct with transition energy and model
     * 
     * @param transition_weight boltzmann weight of transition energy 
     */
    BasinTransition(double transition_weight)
	:
	Z_(transition_weight)
    {}

    /** 
     * Update transition by partition function of additional transition
     *
     * @param Z_ partition functions
     *
     */
    // obsolete note: does not correctly update the barrier
    void
    update(double Z) {
	Z_ += Z;
    }

    /** 
     * Multiply transition partition function
     *
     * @param f multiplication factor
     *
     */
    void
    multiply(double f) {
	Z_ *= f;
    }

    
    /** 
     * Get partition function
     * 
     * @return partition function of transition
     */	
    double
    Z() const {return Z_;}
};

#endif
