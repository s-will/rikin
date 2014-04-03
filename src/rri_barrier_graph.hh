#ifndef RRI_BARRIER_GRAPH_HH
#define RRI_BARRIER_GRAPH_HH

#include "hybrid_ensemble_model.hh"
#include "barrier_graph.hh"

/**
 * @brief Generation of RRI barrier graph
 *
 * The RRI barrier graph is the graph over gradient basins of
 * RRI-mesostates in the hybrid ensemble model.
 *
 * Implements a barriers-like algorithm to create the graph from an
 * energy sorted list of the RRI-mesostates.
 * 
 * @note An adaptive walk is defined as a walk, where the energy of
 * each target state is smaller or equal to the energy of its source
 * state and the corresponding transition state energy is not larger
 * than the source state energy.
 */
class RRIBarrierGraph: public BarrierGraph {

private:
    //! type of state hash, key=encoded state, value=index of assigned basin
    typedef std::tr1::unordered_map< HybEnsModel::StateDescription::code_t,
				     size_t,
				     HybEnsModel::StateDescription::code_t_hash
				     > state_hash_t;
    
    //! hybrid ensemble model
    HybEnsModel model_;

    //! maximum recover energy
    double max_recover_energy_;

    //! whether double sites are supported
    bool consider_double_sites_;

    //! whether gradient walks are used to combine states into basins
    bool gradient_;

    //! hash of the state basin assignment indexed by encoded states
    state_hash_t state_hash_;


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
	       double &energy, size_t lineno,bool binary) const;

public:


    /** 
     * @brief Construct and initialize (e.g. precompute energies)
     * 
     * @param seqA sequence A
     * @param seqB sequence B
     * @param special_open_state whether open state is treated as special state 
     * (that is never merged with other states)
     * @param maxsitesize maximum size of a sequence in an interaction site
     * @param maxsitesize maximum size difference of sequences in an interaction site
     * @param max_recover_energy maximum energy, where neighbors of a state are still recovered
     * @param consider_double_sites whether double interaction sites are allowed
     * @param gradient whether to combine states into basins due to gradient walks
     * @param verbose turns on verbose output 
     * @param debug_out turns on debugging output 
     */
    RRIBarrierGraph(const std::string &seqA, 
		    const std::string &seqB,
		    bool special_open_state,
		    size_t maxsitesize,
		    size_t maxsitesize_diff,
		    double max_recover_energy,
		    bool consider_double_sites,
		    bool gradient,
		    bool verbose,
		    bool debug_out
		    );

    //! @brief destructor
    ~RRIBarrierGraph() {
    }

    /**
     * @brief read states from input stream
     *
     * @param in input stream of energy sorted states
     * @param binary assume binary encoding of states in input
     *
     * Reads an energy-sorted list of states from input stream and
     * constructs the barrier graph for these states. Merges basins
     * that are separated by a barrier less than min_barrier_height
     * above of one of the minima.  Note that the merging can be done
     * during the construction of the barrier graph, since barriers
     * are maxima of origin, target and transition state
     * energies. This causes barrier energies to be sorted.
     *
     */
    void
    read_states(std::istream &in,
		bool binary);

    
    /**
     * @brief access to the hybrid ensemble model of the object
     * @return model
     */
    const HybEnsModel &model() const {
	return model_;
    }

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
    
protected:
    /** 
     * @brief Register a new transition to a target basin in the source basin 
     *
     * 1) Determines the barrier energy and neighbor basin where the
     * basin connects to with lowest barrier.
     *
     * 2) Determines the partition functions of transitions in the barrier graph 
     *
     * @param tr transition
     */
    void
    add_transition( const transition_t &tr );


private:
    /** 
     * @brief Process a move during processing of a state
     * 
     * @param move move to be applied
     * @param source_state the source state
     * @param source_energy energy of the source state
     * @param[out] min_transition_energy minimum of input value and transition energy of move
     * @param[out] min_transition_index set to basin index of target basin if min_transition_energy is changed; set to maximum of size_t if target is not in hash
     * @param[out] trans, vector of transitions; if 0L, dont't record transitions 
     *
     * The method applies the move and looks up the neighbor in the
     * hash.  The neighbor state with lowest transition state over all
     * neighor states is determined, where only neighbors are
     * considered that are in the hash or larger than maxsitesize.
     */
    void
    process_move(const HybEnsModel::Move *move,
		 const HybEnsModel::StateDescription &source_state,
		 const HybEnsModel::energy_t &source_energy,
		 HybEnsModel::energy_t &min_transition_energy,
		 size_t &min_transition_index,
		 std::vector<transition_t> &trans
		 );


    /**
     * @brief Assign state to a basin
     *
     * @param state state to be assigned to a basin
     * @param energy energy of source state
     * @param basin_index index of the basin
     *
     * Assigns the state to the basin, increasing its partition function;
     * registers assignment in state hash
     */
    void
    assign_to_basin(const HybEnsModel::StateDescription &state,
		    const HybEnsModel::energy_t &energy,
		    size_t basin_index);
    
    /**
     * @brief Create new basin with one state 
     *
     * @param state state to be assigned to a basin
     * @param energy energy of source state
     *
     * @return index of the basin
     *
     * Creates new basin, and assigns the state;
     * registers assignment in state hash
     */
    size_t
    create_new_basin(const HybEnsModel::StateDescription &state,
		     const HybEnsModel::energy_t &energy);
    
    /** 
     * @brief register transitions
     * @param source_basin_index index of the source basin
     * @param trans list of transitions from source to neighbor states, which are in the hash
     *
     * For all transitions from the source basin to a different basin,
     * the transition state pf is increased; registers the forward and
     * the backward transition.
     */
    void
    register_transitions(size_t source_basin_index,
			 std::vector<transition_t> &trans); 

    /**
     * @brief Process a single state in the construction of the barrier graph
     * 
     * Assumes that states are processed in ascending order of their energy.
     * Generates the neighbors of source_state and in this way determines,
     * whether the state is a local minimum or belongs to a known basin.
     *
     * @param source_state     state to be processed
     * @param energy           energy of source_state
     * @param force_new_basin  force creation of new basin, if true
     */
    void
    process_state(const HybEnsModel::StateDescription &source_state, 
		  double source_energy);

    

    /**
     * @brief Perform gradient walk and create missing states on the walk
     *
     * @param state source state of the walk
     * @param energy energy of source state
     * @return target basin index of the walk
     *
     * Iteratively walks to the lowest transition energy neighbor of
     * the source state, until the transition energy is larger than
     * the source energy. Furthermor, terminates immediately if the
     * neighbor is in the hash.
     * 
     * The source state and states on the walk are registered in the
     * hash, and assigned to the target basin.
     */
    size_t create_gradient_walk(const HybEnsModel::StateDescription &state, 
				double energy);

};


std::ostream &
operator <<(std::ostream &out,const RRIBarrierGraph::transition_t &t);


#endif // RRI_BARRIER_GRAPH_HH
