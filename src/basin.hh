#ifndef BASIN_HH
#define BASIN_HH

#include <tr1/unordered_map>
#include <iostream>
//#include <limits>

/**
 * @brief A basin in the energy landscape
 * 
 * Represents a node of the barrier graph together with outgoing edges
 * 
 * @note A basin is not necessarily a gradient basin, but more
 * general.  Even the number of states in a basin can be fractional,
 * i.e. basins are not necessarily discrete sets of states. (It may be
 * appropriate to rename the class from Basin to MacroState.)
 */
class Basin {
    
    size_t index_; //<! index of basin 
    double Z_;      //!< partition function

    double states_; //!< size in number of states   
    
    bool merged_;
    
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
	  Z_(boltzmann_weight),
	  states_(1),
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
	
	states_++;
	Z_ += weight;
    }
    
    void 
    merge_in(const Basin &x, double fraction) {
	Z_      += x.Z_ * fraction;
	states_ += x.states_ * fraction;
    }

    double
    Z() const {
	return Z_;
    }
    
    bool
    merged() const {
	return merged_;
    }


    /** @brief mark as merged
     */
    void
    mark_merged() {
	merged_ = true;
    }
    
    size_t
    number_of_states() const {
	return states_; 
    }

    void
    print_header(std::ostream &out) const;

    void
    print(std::ostream &out) const;

};

#endif // BASIN_HH
