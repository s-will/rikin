#ifndef RRI_BARRIER_GRAPH_HH
#define RRI_BARRIER_GRAPH_HH

#include "hybrid_ensemble_model.hh"
#include "barrier_graph.hh"
#include "basin_info.hh"
#include <zlib.h>

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
    typedef HybEnsModel::energy_t energy_t;
    typedef HybEnsModel::StateDescription::code_t code_t;
    
private:
    //! type of state hash, key=encoded state, value=index of assigned basin
    typedef std::unordered_map< code_t,
        size_t,
	HybEnsModel::StateDescription::code_t_hash
	> state_hash_t;
    
    //! hybrid ensemble model
    const HybEnsModel &model_;

    //! maximum recover energy
    double max_recover_energy_;

    //! start of region of interaction sites in A
    size_t region_startA_;
    
    //! end of region of interaction sites in A
    size_t region_endA_;

    //! start of region of interaction sites in B
    size_t region_startB_;
    
    //! end of region of interaction sites in B
    size_t region_endB_;

    //! whether double sites are supported
    bool consider_double_sites_;

    //! whether gradient walks are used to combine states into basins
    bool gradient_;

    //! hash of the state basin assignment indexed by encoded states
    state_hash_t state_hash_;

    //! whether gradient basin information is tracked
    bool track_basins_;
    
    //! vector of additional infos for each basin (if track_) 
    std::vector<BasinInfo> basin_infos_;

    bool verbose_;

    bool debug_out_;

public:


    /** 
     * @brief Construct and initialize (e.g. precompute energies)
     * 
     * @param special_open_state whether open state is treated as special state 
     * (that is never disolved into other states)
     * @param max_recover_energy maximum energy, where neighbors of a state are still recovered
     * @param region_startA 5' end of region in A
     * @param region_endA   3' end of region in A
     * @param region_startA 3' end of region in B
     * @param region_endA   5' end of region in B
     * @param consider_double_sites whether double interaction sites are allowed
     * @param gradient whether to combine states into basins due to gradient walks
     * @param verbose turns on verbose output 
     * @param debug_out turns on debugging output 
     *
     * @todo consider interface change, such that the model is passed to the constructor! This
     * avoids passing through of model specific parameters
     */
    RRIBarrierGraph(const HybEnsModel &model,
		    bool special_open_state,
		    double max_recover_energy,
		    size_t region_startA,
		    size_t region_endA,
		    size_t region_startB,
		    size_t region_endB,
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
    const HybEnsModel &
    model() const {
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
     * @return whether move is valid (within region)
     * The method applies the move and looks up the neighbor in the
     * hash.  The neighbor state with lowest transition state over all
     * neighor states is determined, where only neighbors are
     * considered that are in the hash or larger than maxsitesize.
     */
    bool
    process_move(const HybEnsModel::Move *move,
		 const HybEnsModel::StateDescription &source_state,
		 energy_t source_energy,
		 energy_t &min_transition_energy,
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
		    energy_t energy,
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
		     energy_t energy);
    
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
    
    /** @brief check for state in region
     * @param state the state to check
     * @return whether state is outside of region
     */
    bool
    state_outside_region(const HybEnsModel::StateDescription &state) const;

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

private:

    /**
     * @brief Append new basin to basins list
     *
     * @param basin_index basin index
     * @param state_code  code of the basin minimum state
     * @param energy      energy of the basin minimum
     *
     * @note keeps tracking information up to date
     */
    void 
    push_back_basin(size_t basin_index, const code_t &state_code, energy_t energy);
    
    /**
     * @brief add a state to a basin
     *
     * @param basin_index basin index
     * @param state_code  code of the added state
     * @param energy      energy of the added state
     *
     * @note keeps tracking information up to date
     */
    void
    add_state_to_basin(size_t basin_index, const code_t &state_code, double energy);

public:
    
    /**
     * @brief turn on tracking of basin information
     */
    void
    track_basins();
    
    /**
     * @brief write basin information
     * @param out output stream
     * @return stream
     *
     * write lines 
     * basin_idx no_states ensemble_energy min_mesostate min_mesostate_energy list_of_mesostates,
     * where list_of_mesostates enumerates all mesostates in the basin
     * 
     * the last three entries are available only if basin tracking was turned on (@see track_basins())
     */
    std::ostream &
    write_basin_track(std::ostream &out) const;
    
    /**
     * @brief write basin information
     * @param file gz file handle opened for write 
     * @return file handle; NULL on error
     *
     * @see write_basin_track()
     */
    gzFile
    gzwrite_basin_track(gzFile file) const;
    
private:

    /**
     * @brief write interactions pair probabilities of a single basin
     * @param out output stream
     * @param min_prob minimum probability 
     * @param basin_idx index of basin
     * @return stream
     */
    std::ostream &
    write_basin_ipps_single_basin(std::ostream &out, 
				  double min_prob,
				  size_t basin_idx) const;

public:

    /**
     * @brief write interactions pair probabilities of basins
     * @param out output stream
     * @param min_prob minimum probability
     * @return stream
     *
     * Write list with lines for each basin (@see
     * write_basin_ipps_single_basin()):
     * 
     * basin_index basin_partition_function { i j p_ij }
     *
     * @note the information about the partition function of the basin is as well given in the
     * basin track (@see write_basin_track()); there, in the form of basin ensemble energy
     */
    std::ostream &
    write_basin_ipps(std::ostream &out,
		     double min_prob) const;

    /**
     * @brief write interactions pair probabilities of basins
     * @param file gz file handle opened for write 
     * @param min_prob minimum probability
     * @return file handle; NULL on error
     *
     * @see write_basin_ipps()
     */
    gzFile
    gzwrite_basin_ipps(gzFile file,
		       double min_prob) const;

};


std::ostream &
operator << (std::ostream &out,
	     const RRIBarrierGraph::transition_t &t);


#endif // RRI_BARRIER_GRAPH_HH
