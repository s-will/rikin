#ifndef BASIN_TRANSITION
#define BASIN_TRANSITION

#include <tr1/unordered_map>
#include <iostream>
//#include <limits>

/**
 * @brief A transition between two basins
 *
 * Encapsulates the partition function of the transition state between two 'basins'
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
     * Adds a term to the transition by partition function
     *
     * @param Z additional term for partition function
     */
    void
    add(double Z) {
	Z_ += Z;
    }

    /** 
     * Multiply transition partition function
     *
     * @param f multiplication factor
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
