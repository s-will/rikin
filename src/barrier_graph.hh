#ifndef BARRIER_GRAPH_HH
#define BARRIER_GRAPH_HH

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
    transitions_map_t transitions_;
    
    //! type of state hash, key=encoded state, value=index of assigned basin
    typedef std::tr1::unordered_map< HybEnsModel::StateDescription::code_t,
				     size_t > state_hash_t;
    
    //! hybrid ensemble model
    HybEnsModel model_;
    
    //! whether open state is special or treated like all other input states
    bool special_open_state_;
    
    //! whether double sites are supported
    bool consider_double_sites_;

    //! whether gradient walks are used to combine states into basins
    bool gradient_;

    /* control output */
    bool verbose_;
    bool debug_out_;


    //! vector storing basins of all local minima
    std::vector<Basin> basins_;
    
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
     * @param binary assume binary encoding of states in input
     * @param special_open_state whether open state is treated as special state (that is never merged with other states)
     * @param consider_double_sites whether double interaction sites are allowed
     * @param gradient whether to combine states into basins due to gradient walks
     * @param verbose turns on verbose output 
     * @param debug_out turns on debugging output 
     *
     */
    BarrierGraph(const std::string &seqA, 
		 const std::string &seqB,
		 bool binary,
		 bool special_open_state,
		 bool consider_double_sites,
		 bool gradient,
		 bool verbose,
		 bool debug_out
		 );

    /** @brief total outflow of basin as partition function
	@param x basin
	@return  sum over transition state partition functions
     */
    double
    outflow_pf(const Basin &x) const {
	double total_out=0;
	
	transitions_map_t::const_iterator trs_x_it=transitions_.find(x.idx());
	
	if (transitions_.end() == trs_x_it) return total_out;
	
	const transitions_map_row_t &trs_x = trs_x_it->second;
	
	for (transitions_map_row_t::const_iterator it=trs_x.begin();
	     trs_x.end()!=it; ++it) {
	    
	    if (basins_[it->first].merged()) continue;
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
	double max_outflow = - std::numeric_limits<double>::infinity();

	transitions_map_t::const_iterator trs_x_it=transitions_.find(x.idx());
	
	if (transitions_.end() == trs_x_it) return 0.0;
	
	const transitions_map_row_t &trs_x = trs_x_it->second;

	
	for (transitions_map_row_t::const_iterator it=trs_x.begin();
	     trs_x.end()!=it; ++it) {
	    
	    if (basins_[it->first].merged()) continue;
	    if (it->first==x.idx()) continue;
	    
	    max_outflow = std::max(max_outflow, it->second.get_Z());
        }
	
	return max_outflow/x.get_Z();
    }

    /** 
     * @brief compute partition function of all basins 
     * 
     * @return partition function 
     */
    double
    compute_Z() const;
    
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

private:
    class compBasinIdxs {
	const BarrierGraph &bg_;
    public:
	compBasinIdxs(const BarrierGraph &bg): bg_(bg) {}
	
	bool
	operator() (size_t a,size_t b) const {
	    return bg_.basins_[a].get_Z() < bg_.basins_[b].get_Z();
	}
    };
public:
    
    /** @brief merge basins with large outflow to their neighbors
     *
     * @param max_outflow maximal outflow (sum of rates) where basin
     * is retained; otherwise it is distributed to its neighbors
     *
     * @param min_p_equ minimum equilibrium probability; basins with
     * lower probability are distributed to their neighbors

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
    merge_basins(double max_outflow, double min_p_equ, double min_rate);

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
	for(size_t i=0; i<basins_.size(); ++i) {
	    if (!basins_[i].merged()) {
		total_outflow += outflow(basins_[i]);
	    }
	}
	out << "mean outflow:      "<<total_outflow/num_basins()<<std::endl;
    }

    /** 
     * @brief print treekin-compatible barriers list to stream 
     * 
     * @param out output stream 
     */
    void
    print_barriers(std::ostream &out) const;    

    
    /** 
     * @brief check validity of rates
     * @note assume there are no merged basins, e.g. after reindex()
     * @returns maximum difference between Zs of transitions_[i][j] and transitions_[j][i] 
     */ 
    double
    check_rates() const;
    
    /** 
     * @brief determine connected components
     * @param[out] components mapping of each basin index to a (1-based) component index
     * @return vector of component sizes
     *
     * @note the current code requires that every state has an entry in the hash transitions
     */
    std::vector<size_t>
    connected_components(std::vector<size_t> &components) const;
    
    /** 
     * @brief keep only a single component
     * @param c index of component to keep 
     * @param components mapping of each basin index to its component index
     *
     * Removes all basins i where components[i]!=c
     */
    void
    keep_single_component(size_t c,const std::vector<size_t> &components);

    /** 
     * @brief print treekin-compatible rates matrix to stream 
     * 
     * @param out output stream 
     */
    void
    print_treekin_ratesmatrix(std::ostream &out) const;    

    /**
     * @brief convert state description to dot bracket notation
     * @param sd state description
     */
    std::string
    to_dotbracket(const HybEnsModel::StateDescription &sd) const;

    /** @brief reindex basins, remove merged basins
     */
    void
    reindex();

    /** @brief reindex basins according to vector keep
     */
    void
    reindex(const std::vector<bool> &keep);


    /** 
     * Swap two basin indices
     * 
     * @param x first index
     * @param y second index
     */
    void swap_indices(size_t x, size_t y);

private:
    /** @brief print a double suited as rate in the ratematrix file for treekin
     */
    static 
    std::string
    format_rate_for_treekin(double x);

    
    /** @brief dump barrier graph to stream
	
	such that it can be read in again by read_graph()
    */
    void dump_graph(std::ostream &out) const;
    
    /** @brief read barrier graph from stream
	
	in the format written by dump_graph()
    */
    void read_graph(std::istream &in) const;
    
};

#endif
