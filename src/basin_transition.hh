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
    
    
    // /** @brief construct from stream
    //  */
    // BasinTransition(std::istream &in);

    // /** @brief write to stream (binary)
    //  */
    // std::ostream &
    // write_binary(std::ostream &out);

    
    // /** 
    //  * Update transition by energy and partition function of additional transition
    //  *
    //  * @param transition_energy transition/barrier energy
    //  * @param Z partition functions
    //  */
    // void
    // update(double transition_energy,double Z) {
    // 	Z_ += Z;
    // }

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
     * Get partition function
     * 
     * @return partition function of transition
     */	
    double
    Z() const {return Z_;}
};

#endif
