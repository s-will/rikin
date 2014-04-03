#ifndef BARRIER_GRAPH_HH
#define BARRIER_GRAPH_HH

#include <vector>
#include <set>
#include <iostream>
#include <limits>
#include <tr1/unordered_map>

#include "basin_transition.hh"

/**
 * @brief Graph of basins and barriers
 * 
 * Basins and barriers are characterized by their respective partition
 * functions.  Supports merge of basins by distributing merged basins
 * to their neighbors proportionally to the flow into each neighbor
 * basin.
 *
 * @todo binary file format of barrier graph and pf-file format are equivalent => merge?!
 */
class BarrierGraph
{
protected:

    //! a row of a transition map "matrix"
    typedef 
    std::tr1::unordered_map<size_t,BasinTransition> 
    transitions_map_row_t;
    
    //! map from basin indices to basin transition define as map of
    //! maps, such that we can iterate over first index (equivalently,
    //! second index, in the case of a symmetric matrix; see below!)
    typedef
    std::tr1::unordered_map<size_t,transitions_map_row_t>
    transitions_map_t;

    //! whether first state is special or treated like all other input
    //! states.  This is used to keep the open state basin.
    bool special_first_state_;
    
    /* control output */

    //! verbose output
    bool verbose_;

    //! debugging output
    bool debug_out_;
    
    /**
     * @brief map of all transitions from one basin to another basins,
     * keys=basin indices of source and target 
     * 
     * @note this matrix is
     * symmetric, we maintain the equivalence of symmetric pairs
     * transitions[i][j] and transitions[j][i].
     */
    transitions_map_t transitions_;        
    
    //! vector storing basins of all local minima
    std::vector<Basin> basins_;

protected:
    
    /**
     * @brief minimal constructor for derived classes
     */
    BarrierGraph( bool special_first_state,
		  bool verbose,
		  bool debug_out
		  )
	: special_first_state_(special_first_state),
	  verbose_(verbose),
	  debug_out_(debug_out)
    {
    }
    
public:

    /**
     * @brief Read barrier graph from binary stream
     * @param in input stream
     *
     * Read in format written by write_binary()
     */
    BarrierGraph(std::istream &in,
		 bool special_first_state,
		 bool verbose,
		 bool debug_out
		 );

    /**
     * @brief write barrier graph to stream such that it can be read
     * in again
     * 
     * @param out output stream
    */
    std::ostream &
    write_binary(std::ostream &out) const;


    /** @brief total outflow of basin as partition function
	@param x basin
	@return sum over transition state partition functions
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
    outflow_rate(const Basin &x) const {
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
    
    //! compare basin indices by their associated pfs
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
    
    /** 
     * @brief select basins to be merged bases on their outflow and minimum probability 
     * 
     * @param x0 basin
     * @param max_outflow maximal outflow (sum of rates) where basin
     * is retained; otherwise it is distributed to its neighbors
     *
     * @param min_p_equ minimum equilibrium probability; basins with
     * lower probability are distributed to their neighbors
     * 
     * @returns whether basin should be merged
     */
    bool
    is_to_be_merged(const Basin &x0, double max_outflow, double min_p_equ) const;


    /**
     * @brief merge a single basins to its neighbors
     * @param x0 Basin
     */
    void
    dissolve_basin(Basin &x0);

    /** @brief Prunes the barrier graph
     *
     * Dissolve small or volatile, high-outflow basins. Remove small
     * transitions.
     *
     * @param max_outflow maximal outflow (sum of rates) where basin
     * is retained; otherwise it is distributed to its neighbors
     *
     * @param min_p_equ minimum equilibrium probability; each basin with
     * lower probability is distributed to its neighbors
     *
     * @param min_rate minimum rate; transitions with lower rate (in
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
    prune(double max_outflow, double min_p_equ, double min_rate);
    

    /** 
     * Reduce the set of basins to those in keep set
     * 
     * @param to_keep set of indices of basins to keep
     * @param min_rate minimal rate (for merging)
     *
     * all basins that are not in the set are merged (in decreasing order of their size)
     */
    void
    reduce_basin_set(const std::set<size_t> &to_keep, double min_rate);
    

    /**
     * @brief filter transitions of a single basin by minimum rate, also remove self-transitions
     * @param x0 basin
     * @param min_rate minimum rate; transition with lower rate (in
     * both directions!) are removed
     */
    void
    filter_basin_transitions(const Basin &x0, double min_rate);
    
    /**
     * @brief filter transitions by minimum rate, also remove self-transitions
     * @param min_rate minimum rate; transition with lower rate (in
     * both directions!) are removed
     *
     */
    void filter_transitions(double min_rate);
    

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

    
    //! @brief destructor
    ~BarrierGraph() {
    }

    /** @brief print some statistics of the barrier graph
	@param ostream output stream
     */ 
    std::ostream &
    print_stats(std::ostream &out) const;
        
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
     * @brief print partition functions of basins and transitions to stream 
     * 
     * Prints full matrix of transition pfs and sets diagonal to basin pfs 
     *
     * @param out output stream 
     * @param binary write binary data
     */
    void
    print_pfs(std::ostream &out,bool binary=false) const;

    /** 
     * @brief print reactions to stream
     * 
     * Prints reactions in format 'rxns' (compatible to Stochastinator)
     *
     * @param out output stream
     * @param nameA name of homomer A
     * @param nameB name of homomer B
     * @param nameAB name of dimer
     */
    void
    print_rxns(std::ostream &out,
	       const std::string &nameA,
	       const std::string &nameB,
	       const std::string &nameAB) const;

    /** 
     * @brief print species to stream
     * 
     * Prints reactions in format 'spcs' (compatible to Stochastinator)
     *
     * @param out output stream
     * @param nameA name of homomer A
     * @param nameB name of homomer B
     * @param nameAB name of dimer
     */
    void
    print_spcs(std::ostream &out,
	       const std::string &nameA,
	       const std::string &nameB,
	       const std::string &nameAB) const;
    
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


    /** 
     * @brief print treekin-compatible barriers list to stream 
     * 
     * @param out output stream 
     */
    void
    print_treekin_barriers(std::ostream &out, std::string header, double RT) const;


private:
    /** @brief print a double suited as rate in the ratematrix file for treekin
     */
    static 
    std::string
    format_rate_for_treekin(double x);

    
    /** 
     * @brief determine outflow distribution
     * 
     * @param outflows[out] vector of sorted outflows
     */
    void
    outflow_distribution(std::vector<double> &outflows) const;    
    
    /** 
     * @brief determine equilibrium probability distribution
     * 
     * @param pequs[out] vector of sorted p_equs
     */
    void
    pequ_distribution(std::vector<double> &pequs) const;

public:
    /** 
     * @brief maximum outflow by quantile
     * 
     * @param q quantile of basins to keep
     * 
     * @return max outflow
     */
    double
    max_outflow_by_quantile(double q) const;
    
    /** 
     * @brief maximum outflow by number
     * 
     * @param n number of basins to keep
     * 
     * @return max outflow
     */
    double
    max_outflow_by_number(size_t n) const;
    
    /** 
     * @brief mimum p_equ by quantile
     * 
     * @param q quantile of basins to keep
     * 
     * @return min pequ
     */
    double
    min_pequ_by_quantile(double q) const; 
    
    /** 
     * @brief mimum p_equ by number
     * 
     * @param q number of basins to keep
     * 
     * @return min pequ
     */
    double
    min_pequ_by_number(size_t n) const;
    
};


#endif
