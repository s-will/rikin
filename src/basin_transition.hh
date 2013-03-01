#ifndef BASIN_TRANSITION
#define BASIN_TRANSITION

#include "hybrid_ensemble_model.hh"
#include <tr1/unordered_map>
#include  <iostream>
#include <limits>


class BarrierGraph;

/**
 * @brief A transition between two basins
 *
 * @note Use this class to store information (barrier, Z) for
 * transitions between two basins in the map of transitions stored in
 * each basin. In contrast, Basin::transition_t stores transitions
 * between (micro)states.
 */
class BasinTransition {
    
    //double barrier; //! minimal energy of a transition state
    double Z;       //! partition function of macro transition state
public:
    BasinTransition(): 
	//barrier(std::numeric_limits<double>::infinity()),
	Z(0) {}
    
    /** 
     * Construct with transition energy and model
     * 
     * @param transition_energy 
     * @param model 
     */
    BasinTransition(double transition_energy,const HybEnsModel &model)
	:
	//barrier(transition_energy),
	Z(model.boltzmann_weight(transition_energy))
    {}
    
    /** 
     * Update transition with energy and model
     * 
     * @param transition_energy transition/barrier energy
     * @param model hybrid ensemble model
     */
    void
    update(double transition_energy,const HybEnsModel &model) {
	//barrier = std::min(barrier,transition_energy);
	Z+=model.boltzmann_weight(transition_energy);
    }
    
    /** 
     * Update transition by energy and partition function of additional transition
     *
     * @param transition_energy transition/barrier energy
     * @param Z_ partition functions
     */
    void
    update(double transition_energy,double Z_) {
	//barrier = std::min(barrier,transition_energy);
	Z+=Z_;
    }

    /** 
     * Update transition by partition function of additional transition
     *
     * @param Z_ partition functions
     *
     */
    // obsolete note: does not correctly update the barrier
    void
    update(double Z_) {
	Z+=Z_;
    }

    
    /* 
     * Get barrier
     *  
     * @return barrier of transition
     
    double
    get_barrier() const {return barrier;}
    */

    /** 
     * Get partition function
     * 
     * @return partition function of transition
     */	
    double
    get_Z() const {return Z;}
};



/**
 * @brief A basin in the energy landscape
 * 
 * Represents a node of the barrier graph together with outgoing edges
 * 
 */
class Basin {
    
    //friend class BarrierGraph;
    
    size_t index_; //<! index of basin 
    HybEnsModel::StateDescription local_minimum; //!< local minimum
    double minimum_energy; //!< energy of the local minimum
    double Z;      //!< partition function

    double states; //!< size in number of states   
    
    bool merged_;
    
public:
       
    /** 
     * Construct a new basin with given index and local minimum
     * 
     * @param index The index of the basin
     * @param local_minimum_ The local minimum of the basin
     * @param energy_of_minimum The energy of the local minimum
     */
    Basin(size_t index, const HybEnsModel::StateDescription &local_minimum_,double energy_of_minimum, const HybEnsModel &model)
	: index_(index),
	  local_minimum(local_minimum_),
	  minimum_energy(energy_of_minimum),
	  Z(0.0),
	  states(0),
	  merged_(false)
    {
	add_state(energy_of_minimum, model);
    }
    
    
    /** 
     * Construct undefined xbasin
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
    
    void
    add_state(double energy, const HybEnsModel &model) {
	//assert(energy>=minimum_energy);
	
	states++;
	Z += model.boltzmann_weight(energy);
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

    HybEnsModel::StateDescription
    get_local_minimum() const {
	return local_minimum;
    }

    double
    get_Z() const {
	return Z;
    }
	
    double
    get_minimum_energy() const {
	return minimum_energy;
    }
    
    bool
    merged() const {
	return merged_;
    }

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
    print(std::ostream &out, const HybEnsModel &model) const;
};

#endif
