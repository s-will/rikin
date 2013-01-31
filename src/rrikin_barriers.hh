#ifndef RRIKIN_BARRIERS_HH
#define RRIKIN_BARRIERS_HH

#include  <iostream>
#include  <string>
//#include <LocARNA/matrices.hh>
#include <LocARNA/stopwatch.hh>

#include <tr1/unordered_map>

#include <limits>

#include "hybrid_ensemble_model.hh"

/* control output */
const bool debug_out=false;
const bool verbose=true;

/* control behavior */
bool simplify_graph=true;


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
    friend class BarrierGraph;

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
    
    friend class BarrierGraph;
    
    size_t basin_index; //<! index of basin 
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
	: basin_index(index),
	  local_minimum(local_minimum_),
	  minimum_energy(energy_of_minimum),
	  Z(0.0),
	  states(0),
	  merged_(false)
    {
	add_state(energy_of_minimum, model);
    }

    //! @brief Get index of basin
    //! @return index of basin
    size_t idx() const {
	return basin_index;
    }
    
    void
    add_state(double energy, const HybEnsModel &model) {
	assert(energy>=minimum_energy);
	
	states++;
	Z += model.boltzmann_weight(energy);
    }
    
    void 
    merge_in(const Basin &x, double fraction) {
	Z+=x.Z*fraction;
	states+=x.states*fraction;
	if (x.minimum_energy < minimum_energy) local_minimum=x.local_minimum;
	minimum_energy=std::min(minimum_energy,x.minimum_energy);
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
    
    size_t
    number_of_states() const {
	return states; 
    }

    
    void
    print_header(std::ostream &out) const;

    void
    print(std::ostream &out, const HybEnsModel &model) const;
};


/**
 * Implements generation of barrier graph and maintains the graph
 * 
 * @note An adaptive walk is defined as a walk, where the energy of
 * each target state is smaller or equal to the energy of its source
 * state and the corresponding transition state energy is not larger
 * than the source state energy.
 */
class BarrierGraph
{

private:
    
    typedef std::tr1::unordered_map<size_t,BasinTransition>
    transitions_map_row_t;
    
    
    //! map from basin indices to basin transition define as map of
    //! maps, such that we can iterate over first index (equivalently,
    //! second index, in the case of a symmetric matrix; see below!)
    typedef
    std::tr1::unordered_map<size_t,transitions_map_row_t>
    transitions_map_t;
    
    //! map of all transitions from this basin to other basins,
    //! keys=basin indices of source and target @note this matrix is
    //! symmetric, we maintain the equivalence of symmetric pairs
    //! transitions[i][j] and transitions[j][i].
    transitions_map_t transitions;
    
    //! type of state hash, key=encoded state, value=index of assigned basin
    typedef std::tr1::unordered_map< std::string, size_t > state_hash_t;
    
    //! hybrid ensemble model
    HybEnsModel model;
    
    //! vector storing basins of all local minima
    std::vector<Basin> basins;
    
    //! hash of the state basin assignment indexed by encoded states
    state_hash_t state_hash;
    
    /** 
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
    bool
    read_state(std::istream &in, HybEnsModel::StateDescription &state, 
	       double &energy, size_t lineno) const;

public:
    /**
     * A transition between two microstates
     * 
     */
    struct transition_t {
	size_t source_basin_index;  //!< index of source basin
	double source_state_energy; //!< energy of source state
	size_t target_basin_index;  //!< index of target basin
	double target_state_energy; //!< energy of target state
	double transition_energy;   //!< energy of transition
	
	transition_t(double source_state_energy_,
		     size_t target_basin_index_,
		     double target_state_energy_,
		     double transition_energy_)
	    :source_basin_index(-1),                       // init with invalid source basin index
	     source_state_energy(source_state_energy_),
	     target_basin_index(target_basin_index_),
	     target_state_energy(target_state_energy_),
	     transition_energy(transition_energy_)
	{}

	transition_t(size_t source_basin_index_,
		     double source_state_energy_,
		     size_t target_basin_index_,
		     double target_state_energy_,
		     double transition_energy_)
	    :source_basin_index(source_basin_index_),
	     source_state_energy(source_state_energy_),
	     target_basin_index(target_basin_index_),
	     target_state_energy(target_state_energy_),
	     transition_energy(transition_energy_)
	{}


	transition_t
	reverse() const {
	    return transition_t(target_basin_index,
				target_state_energy,
				source_basin_index,
				source_state_energy,
				transition_energy);
	}
	
	double
	barrier_energy() const {
	    return std::max(source_state_energy,std::max(target_state_energy,transition_energy));
	}
    };
    

private:	

    /** 
     * @brief Register a new transition to a target basin in the source basin 
     *
     * 1) Determines the barrier energy and neighbor basin where the
     * basin connects to with lowest barrier.
     *
     * 2) Determines the partition functions of transitions in the barrier graph 
     *
     * @param tr transition (containing 
     */
    void
    add_transition( const transition_t &tr, const HybEnsModel &model);
    
public:


    /** 
     * @brief Generate barrier graph
     * 
     * Reads an energy-sorted list of states from std::cin and
     * constructs the barrier graph for these states. Merges basins
     * that are separated by a barrier less than min_barrier_height
     * above of one of the minima.  Note that the merging can be done
     * during the construction of the barrier graph, since barriers
     * are maxima of origin, target and transition state
     * energies. This causes barrier energies to be sorted.
     *
     * @param seqA sequence A
     * @param seqB sequence B
     *
     * @param min_barrier_height_ initializes BarrierGraph::min_barrier_height 
     *
     */
    BarrierGraph(const std::string &seqA, const std::string &seqB);

    /** @brief total outflow of basin as partition function
	@param x basin
	@return  sum over transition state partition functions
     */
    double
    outflow_pf(const Basin &x) const {
	const transitions_map_row_t &trs_x=transitions.find(x.idx())->second;
	double total_out=0;
	
	for (transitions_map_row_t::const_iterator it=trs_x.begin();
	     trs_x.end()!=it; ++it) {
	    
	    if (basins[it->first].merged()) continue;
	    if (it->first==x.idx()) continue;
	    
	    total_out += it->second.get_Z();
        }
	return total_out;
    }

    /** @brief total outflow as rate
	@param x basin
	@return sum of outflow rates
     */
    double
    outflow(const Basin &x) const {
	return outflow_pf(x)/x.get_Z();
    }

    /** @brief maximum outflow rate
	@param x basin
	@return maximum outflow rate of basin
    */
    double
    max_outflow(const Basin &x) const {
	const transitions_map_row_t &trs_x=transitions.find(x.idx())->second;
	double max_outflow = - std::numeric_limits<double>::infinity();
	
	for (transitions_map_row_t::const_iterator it=trs_x.begin();
	     trs_x.end()!=it; ++it) {
	    
	    if (basins[it->first].merged()) continue;
	    if (it->first==x.idx()) continue;
	    
	    max_outflow = std::max(max_outflow, it->second.get_Z());
        }
	
	return max_outflow/x.get_Z();
    }
    
    /**
     * @brief Process a single state in the construction of the barrier graph
     * 
     * Assumes that states are processed in ascending order of their energy.
     * Generates the neighbors of source_state and in this way determines,
     * whether the state is a local minimum or belongs to a known basin.
     *
     * @param source_state State that is processed
     * @param energy     Energy of source_state
     */
    void
    process_state(const HybEnsModel::StateDescription &source_state, double source_energy);

private:
    class compBasinIdxs {
	const BarrierGraph &bg_;
    public:
	compBasinIdxs(const BarrierGraph &bg): bg_(bg) {}
	
	bool
	operator() (size_t a,size_t b) const {
	    return bg_.basins[a].get_Z() < bg_.basins[b].get_Z();
	}
    };
public:
    
    /** @brief merge basins with large outflow to their neighbors
     *
     * @param max_outflow maximal outflow (sum of rates) where basin
     * is retained; otherwise it is distributed to its neighbors
     *
     * @param min_rate minimum rate; transition with lower rate (in
     * both directions!) are removed during the merging
     *
     * @note ATTENTION: the min_rate criterion is applied during the
     * merging phase, where rates can still be increased due to
     * merging of larger states.  However merges from larger states
     * will usually not significantly increase rates, such that the
     * error is low. Nevertheless, this is only a heuristic.
     *
     * @note merging of larger states can still produce transitions
     * with rates lower min_rate between smaller states.
     *
     */
    void
    merge_basins(double max_outflow, double min_rate);

    /**
     * @brief filter rates by minimum rate, also remove self-transitions
     * @param min_rate minimum rate; transition with lower rate (in
     * both directions!) are removed
     *
     */
    void filter_rates(double min_rate);
    
    
    /** 
     * @brief Print the list of all basins of the barriers graph 
     * 
     * Output is sorted by basin index
     *
     * @param out output stream
     */
    void
    print_basins(std::ostream &out) const;

    void
    print_edges(std::ostream &out,const Basin &b) const;

    /**
       @brief print representation of barrier graph with weights of basins and transitions
       @param ostream output stream
    */
    void 
    print_barrier_graph(std::ostream &out) const;
    
    /**
       @brief number of (non-merged) basins
    */
    size_t
    num_basins() const;

    /**
       @brief number of transitions (only from non-merged states)
    */
    size_t
    num_transitions() const;

    
    ~BarrierGraph() {
    }

    /** @brief print some statistics of the barrier graph
	@param ostream output stream
     */ 
    void
    print_stats(std::ostream &out) {
	out << "basins:            "<<num_basins()<<std::endl;
	out << "mean transitions:  "<<num_transitions()/num_basins()<<std::endl;
	
	double total_outflow=0.0;
	for(size_t i=0; i<basins.size(); ++i) {
	    if (!basins[i].merged()) {
		total_outflow += outflow(basins[i]);
	    }
	}
	out << "mean outflow:      "<<total_outflow/num_basins()<<std::endl;
    }
    
};



#endif //  RRIKIN_BARRIERS_HH
