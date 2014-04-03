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
    get_Z() const {return Z_;}
};



/**
 * @brief A basin in the energy landscape
 * 
 * Represents a node of the barrier graph together with outgoing edges
 * 
 */
class Basin {
    
    size_t index_; //<! index of basin 
    double Z;      //!< partition function

    double states; //!< size in number of states   
    
    bool merged_;
    
    // double minimum_energy; //!< energy of the local minimum
    // HybEnsModel::StateDescription local_minimum; //!< local minimum
    
public:
       
    /** 
     * Construct a new basin with given index and local minimum
     * 
     * @param index The index of the basin
     * @param boltzmann_weight the initial Boltzmann weight
     */
    Basin(size_t index,
	  double boltzmann_weight)
	: index_(index),
	  Z(boltzmann_weight),
	  states(1),
	  merged_(false)
    {
    }
    
    
    /** 
     * Construct undefined basin
     */
    Basin() : merged_(true)
    {
    }

    //! @brief Get index of basin
    //! @return index of basin
    size_t idx() const {
	return index_;
    }

    //! @brief Set index of basin
    //! @param idx index of basin
    void
    set_idx(size_t idx) {
	index_=idx;
    }
    
    /**
     * @brief add a state to basin
     * @param weight Boltzmann weight of the state
     */
    void
    add_state(double weight) {
	//assert(energy>=minimum_energy);
	
	states++;
	Z += weight;
    }
    
    void 
    merge_in(const Basin &x, double fraction) {
	Z+=x.Z*fraction;
	states+=x.states*fraction;
	// note: if we merge in a basin with lower minimum it still
	// does not always make sense to set its local minimum as new
	// minimum, if the fraction is too small.  To handle this
	// correctly, we could compare the ensemble energy to -RTln of the fraction
	// times Z of the merged in basin.
	// TOO SIMPLISTIC is therefore:
	// if (x.minimum_energy < minimum_energy)
	//    local_minimum=x.local_minimum;
	// minimum_energy=std::min(minimum_energy,x.minimum_energy);
    }

    double
    get_Z() const {
	return Z;
    }
    
    bool
    merged() const {
	return merged_;
    }

    /**
     * @brief Is basin mergeable
     * @note use this to forbid merge of special states, e.g. open state
     */
    bool
    mergeable() const;

    /** @brief mark as merged
     */
    void
    mark_merged() {
	merged_=true;
    }
    
    size_t
    number_of_states() const {
	return states; 
    }

    void
    print_header(std::ostream &out) const;

    void
    print(std::ostream &out) const;

};

#endif
