/**
 * @file rrikin_barriers.cc
 *
 * Defines main() of the program rrikin_barriers
 *
 * rrikin_barriers constructs the barrier tree/graph and rate matrix for the
 * macro state process that moves between basins in the energy
 * landscape.
 *
 * @note Due to our non-standard definition of transition energies, it
 * becomes necessary to rethink the definition of local minima and and
 * basins in the energy landscape!  Usually these are defined in terms
 * of (non-transition) state energies, comparing the energies of
 * source and target state in each move/transition.  In general, these
 * notions should be defined in terms of transition state energies,
 * since these energies govern the speed of transitions.  We define a
 * neighborship relation x ->_t y with the semantics that there is a
 * move from x to y with transition state t.  For a state (set of
 * states/ensemble) s, E_s denotes it's energy (ensemble energy) and
 * Z_s denotes it's Boltzmann weight (partition function).  A state x
 * is a local minimum, iff forall x ->_t y: E_x <= E_t. Note that the
 * common definition puts E_y in place of E_t. Analogously, the
 * criterion for adaptive walks and steepest descent walk (i.e. basin
 * assignment) is modified to compare to the transition state energy
 * in place of the target state energy.
 *
 * IS THIS TRUE?: Note that both (i.e., old and new) local minimum
 * definitions are equivalent for Kawasaki and Metropolis rates. The
 * steepest descent definition is equivalent for the case of Kawasaki
 * rates but not *equivalent* for Metropolis. This could indicate a
 * problem of the new (or the old) definition!?: For x->_t y, where
 * E_t<=E_x, which implies E_y<=E_x, the comparison to E_t does not
 * differentiate in the case of Metropolis rates, since there E_t=E_x!
 *
 * @note In our setting, for all x->_t y the property holds that E_t
 * >= min(E_x,E_y). [Note that this is a consequence of either the
 * source or the target being a relaxation of the transition state and
 * Prop.: For energies of ensembles x and y, there holds that x
 * subset y implies E_x >= E_y.] This important property allows us to
 * use a barrier-like algorithm, where states are processed in the
 * order of increasing energy, for identifying local minima and
 * assigning basin membership to all states.
 *
 * @note For the computation of macro-state rates for the transitions
 * between basins it suffices to accumulate the partition functions of
 * basins and transition states between basins. Note that the
 * definition e.g. given by Flamm, Hofacker, 2007 is equivalent to the
 * Arrhenius rate derived from these partition functions, i.e.
 * k(a->b) = sum_x in a,y in b Pr[x|a]k(x->y) = Z^trans_ab / Z_a,
 * where Z^trans_ab = sum_x in a,y in b: Z^trans_xy.
 *
 * @note Although the model is degenerate, steepest descent walks are 
 * uniquely defined, since the order on moves is fixed
 *
 * @note Barriers/transition states between basins: the energy of the
 * barrier has to be the maximum of source state energy, transition
 * state energy, and target state energy
 * (cf. transition_t::barrier_energy()). Otherwise, the barrier can be
 * lower than one of the explicit states. In this case the barrier
 * would not reflect that the energy has to be raised to this explicit
 * state for the transition between the basins.
 *
 * @todo CHECK: is the process irreducibel?
 * @todo accumulate transition state pf, calculate rates, print rate matrix
 */

#include  <iostream>
#include  <sstream>
#include  <limits>

#include  <stdlib.h>
#include  <string.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/fold_vars.h" // defines global variables
}

// #ifdef _OPENMP
// #include <omp.h>
// #endif

#include "rrikin_barriers_cmdline.h"

#include <LocARNA/matrices.hh>
#include "hybrid_ensemble_model.hh"

#include <tr1/unordered_map>


const bool debug_out=false;

class BarrierGraph;


/**
 * @brief A transition between two basins
 * 
 */
class BasinTransition {
    double barrier; //! minimal energy of a transition state
    double Z;       //! partition function of macro transition state
public:
    BasinTransition(): barrier(std::numeric_limits<double>::max()),Z(0) {} 
    
    /** 
     * Construct with transition energy and model
     * 
     * @param transition_energy 
     * @param model 
     */
    BasinTransition(double transition_energy,const HybEnsModel &model)
	:barrier(transition_energy),Z(model.boltzmann_weight(transition_energy))
    {}
    
    /** 
     * Update transition with energy and model
     * 
     * @param transition_energy transition/barrier energy
     * @param model hybrid ensemble model
     */
    void
    update(double transition_energy,const HybEnsModel &model) {
	barrier = std::min(barrier,transition_energy);
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
	barrier = std::min(barrier,transition_energy);
	Z+=Z_;
    }
    
    /** 
     * Update transition adding new BasinTransition
     * 
     * @param transition addition trransition 
     */
    void update(const BasinTransition &t) {
	update(t.barrier,t.Z);
    }
    
    /** 
     * Get barrier
     *  
     * @return barrier of transition
     */
    double
    get_barrier() const {return barrier;}
    
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
    size_t basin_index; //<! index of basin 
    HybEnsModel::StateDescription local_minimum; //!< local minimum
    double minimum_energy; //!< energy of the local minimum
    double Z;      //!< partition function
    size_t states; //!< size in number of states   
	
    double barrier;    //!< energy of minimal energy transition out of the basin
    size_t connect_to; //!< basin index of deepest basin accessible under barrier

public:
    //! if the basin is merged to a basin, index of the other basin. Otherwise, equal to this.basin_index  
    size_t merged_to;

private:
    //! map from basin index to basin transition 
    typedef std::tr1::unordered_map<size_t,BasinTransition> transitions_map_t;
	
    //! map of all transitions from this basin to other basins, key=target basin index
    transitions_map_t transitions;
public:
    /**
     * A transition between two states
     * 
     */
    struct transition_t {
	double source_state_energy; //!< energy of source state
	size_t target_basin_index;  //!< index of target basin
	double target_state_energy; //!< energy of target state
	double transition_energy;   //!< energy of transition
	    
	transition_t(double source_state_energy_,
		     size_t target_basin_index_,
		     double target_state_energy_,
		     double transition_energy_)
	    :source_state_energy(source_state_energy_),
	     target_basin_index(target_basin_index_),
	     target_state_energy(target_state_energy_),
	     transition_energy(transition_energy_)
	{}
	    
	transition_t
	reverse(size_t source_basin_index) const {
	    return transition_t(target_state_energy,source_basin_index,source_state_energy,transition_energy);
	}
	
	double
	barrier_energy() const {
	    return std::max(source_state_energy,std::max(target_state_energy,transition_energy));
	    //return std::max(source_state_energy,target_state_energy);
	}
    };
    
    
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
	  barrier(std::numeric_limits<double>::infinity()),
	  connect_to(index),
	  merged_to(index)
    {
	add_state(energy_of_minimum, model);
    }
    
    void
    add_state(double energy, const HybEnsModel &model) {
	states++;
	Z += model.boltzmann_weight(energy);
    }

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
    add_transition( const transition_t &tr, const HybEnsModel &model) {
	// keep track of minimal energy transition from this basin to some smaller basin
	if ( (basin_index > tr.target_basin_index) && tr.barrier_energy() < barrier ) {
	    if (debug_out) std::cerr << "  New barrier for basin "<<basin_index<<" to basin "<< tr.target_basin_index  <<": "<< tr.barrier_energy()<<std::endl;
	    barrier = tr.barrier_energy();
	    connect_to = tr.target_basin_index;
	}
    
	// add transition to partition function for the transition
	// between the source and target basin
	if (transitions.find(tr.target_basin_index) == transitions.end()) {
	    transitions[tr.target_basin_index]=
		BasinTransition(tr.barrier_energy(),model);
	} else {
	    transitions[tr.target_basin_index].update(tr.barrier_energy(),model);
	}
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
    
    size_t get_merged_to() const {
	return merged_to;
    }

    /** 
     * @brief Merge in a second basin
     *
     * Transfer all states and transitions from the given basin to this basin.
     * Don't transfer transition from from_basin to this basin!
     *
     * @todo Update transitions in all other basins!
     *
     * @param from_basin basin that should be merged in
     */
    void
    merge_in_basin(Basin &from_basin, const BarrierGraph &bg) {
	assert(minimum_energy <= from_basin.minimum_energy);
		
	// mark from basin as merged
	from_basin.merged_to = basin_index;
	
	// transfer partition function and number of states from from_basin to to_basin
	states += from_basin.states;
	Z += from_basin.Z;
	
	/*
	// transfer transitions
	for(transitions_map_t::const_iterator it=from_basin.transitions.begin();
	    it != from_basin.transitions.end();
	    ++it) {
	    
	    size_t target_basin_index=it->first;
	    
	    // delete transitions from from_basin to this basin
	    if (target_basin_index==this->basin_index) continue;
	    
	    const BasinTransition &transition=it->second;
	    
	    if ( transitions.find(target_basin_index) == transitions.end() ) {
		transitions[target_basin_index]=transition;
	    } else {
		transitions[target_basin_index].update(transition.get_barrier(),transition.get_Z());
	    }
	}
	*/
    }
    
    size_t
    number_of_states() const {
	return states; 
    }

    
    void print_header(std::ostream &out) const {
	printf("%5s %-30s %4s %6s %5s %6s %6s\n",
	       "idx",
	       "description",
	       "n_s",
	       "ensE",
	       "ct",
	       "minE",
	       "h"
	       );
    }

    void
    print(std::ostream &out, const HybEnsModel &model) const {
	printf("%5lu %-30s %4lu %6.2f %5lu %6.2f %6.2f\n",
	       basin_index,
	       local_minimum.toString().c_str(),
	       states,
	       - model.RT() * log(Z),
	       connect_to,
	       minimum_energy,
	       barrier-minimum_energy
	       );
    }

    void
    print_edges(std::ostream &out) const {
	out << basin_index << "  -> ";
	
	for(transitions_map_t::const_iterator it=transitions.begin();
	    it != transitions.end();
	    ++it) {
	    out << " " << it->first << " " << (it->second.get_Z() / get_Z());
	}
	out << std::endl;
    }
};


std::ostream & operator <<(std::ostream &out,const Basin::transition_t &t) {
    out << "Transition to "<<t.target_basin_index
	<< " (sourceE="<<t.source_state_energy
	<< " tgtE="<<t.target_state_energy
	<< " transE="<<t.transition_energy
	<< " transE-sourceE="<<(t.transition_energy-t.source_state_energy)<<")";
    return out;
}


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
    //! type of state hash, key=encoded state, value=index of assigned basin
    typedef std::tr1::unordered_map< std::string, size_t > state_hash_t;
    
    //! hybrid ensemble model
    HybEnsModel model;
    
    /** minimal barrier height between basins. Basins that are
     *  seperated by only lower barrier height are merged.
     */
    double min_barrier_height; 

    
    std::vector<Basin> basins;
    
    //! hash of the state basin assignment indexed by encoded states
    state_hash_t state_hash;
    
    /** 
     * Read a state description line from input stream
     * 
     * @param in input stream 
     * @param[out] state state description 
     * @param[out] energy state energy
     * 
     * @return whether valid state description line could be read
     */
    bool 
    read_state(std::istream &in, HybEnsModel::StateDescription &state, double &energy) const {
	std::string line;
	
	if (!getline(in,line)) return false;
	
	std::istringstream linein(line);
	linein >> energy;
	
	std::vector<size_t> state_vec;
	for (size_t i; linein >> i;) {
	    state_vec.push_back(i);
	}
	
	switch ( state_vec.size() ) {
	case 0:
	    state = HybEnsModel::StateDescription();
	    break;

        case 4:
            state = HybEnsModel::StateDescription(state_vec[0],
                                                  state_vec[1],
                                                  state_vec[2],
                                                  state_vec[3]);
            break;
        case 8:
            state = HybEnsModel::StateDescription(state_vec[0],
                                                  state_vec[1],
                                                  state_vec[2],
                                                  state_vec[3],
                                                  state_vec[4],
                                                  state_vec[5],
                                                  state_vec[6],
                                                  state_vec[7]);
            break;
        default:
            std::cerr << "ERROR: ignored invalid input line" << std::endl;
            return false;
        }
        
        if (not state.is_valid(model)) {
            std::cerr << "Input error: read state "<<state<<" is not valid in model."<<std::endl;
            return false;
        }
        
        return true;
    }

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
     * @param max_adaptive_walk_delta_ initializes BarrierGraph::max_adaptive_walk_delta
     *
     */
    BarrierGraph(const std::string &seqA, const std::string &seqB,
		 double min_barrier_height_
             )
        : model(seqA,seqB),
          min_barrier_height(min_barrier_height_)
    {
        
        HybEnsModel::StateDescription orig_state;
        double energy;
        
        // counter for states that are read for construction of graph
        size_t state_counter;

        state_counter=0;
        
        while (read_state(std::cin,orig_state,energy)) {
            
            if (debug_out) std::cerr << "read " << state_counter << " " << energy << " " << " "  << orig_state << std::endl;
            
            process_state(orig_state,energy);
            
            state_counter++;
        }
    }



    /** 
     * @brief Generate newick representation of barrier tree
     *
     * Generates the barrier tree in newick format from already
     * computed information in array basins
     *
     * @return barrier tree in newick format
     * @todo implement
     */
    std::string barrier_tree() {
	std::string tree="";
	std::cerr << "Generating newick representation of barrier tree no implemented yet."
		  << "Thus, return empty string." << std::endl;
	
	return tree;
    }

    /** 
     * Resolve basin index taking potential basin merge into account
     * 
     * @param index basin index 
     * 
     * @return true basin index, resolving potential merge(s)
     *
     * @todo this could short-circuit to save time in previous resolvings 
     */
    size_t resolve_basin_index(size_t index) const {
	while (basins[index].get_merged_to() != index) {
	    index = basins[index].get_merged_to();
	}
	return index;
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
    process_state(const HybEnsModel::StateDescription &source_state, double source_energy) {
	
	// minimum transition state energy from source state source_state to a neighbor
	double min_transition_energy = std::numeric_limits<double>::infinity();
	
	// index of the basin of the neighbor state with minimum transition energy
	size_t min_tE_target_basin_index = std::numeric_limits<size_t>::max();

	std::vector<Basin::transition_t> transitions;

	size_t moves_counter=0;

	// ----------------------------------------
	// Enumerate neighbors of source_state
	//
	// Register all transitions from source_state to previously read neighbor states.
	// Keep track of neighbor state with smallest transition energy 
	//
	HybEnsModel::MoveIterator mi(source_state,model);
	for (HybEnsModel::Move *move = mi.firstMove(); move != NULL; move = mi.nextMove(move)) {
	    
	    //std::cout << " move "; move->print(std::cout); std::cout<<std::endl;
	    
	    moves_counter++;
	    
	    HybEnsModel::energy_t tE=move->transitionEnergy();
	    
	    HybEnsModel::StateDescription neigh_state=source_state;
	    move->apply(neigh_state);
	    
// #ifndef NDEBUG
// 	    // check whether generated neighbor state is a valid state in the model
// 	    if( not neigh_state.is_valid(model) ) {
		
// 		std::cerr
// 		    << "ERROR: Move produces invalid state" << std::endl 
// 		    << source_state << " == ";
// 		move->print(std::cerr); 
// 		std::cerr << " ==> " << neigh_state << std::endl;
// 		abort();
// 	    }
// #endif
	    
	    // encode neighbor and search neighbor code in hash
	    
	    std::string neigh_code; // string for holding code
	    neigh_state.encode(neigh_code);
	    
	    state_hash_t::const_iterator it = state_hash.find(neigh_code);
	    if (state_hash.end() != it) {
		//found => belongs to basin of already seen local minimum
		
		// NOTE: in the transition from source_state to
		// neigh_state, neigh_state has lower or equal energy
		// than the source_state! The transition state has higher
		// or equal energy than neigh_state (and possibly source_state)
		
		// NOTE: the following differs from the standard
		// barriers algorithm, since there transition energies
		// are sorted in the same order as target state
		// energies
		
		size_t target_basin_index=resolve_basin_index(it->second);
		
		// if transition energy to neigh_state is smaller than the former minimum
		if ( tE < min_transition_energy ) {
		    // record new minimal energy transition
		    min_transition_energy = tE;
		    min_tE_target_basin_index = target_basin_index;
		}
		
		// compute energy of neighbor state
		double neigh_energy = model.energy(neigh_state);
		
		// record new transition.
		// In this way, we collect all transitions from
		// the source state to energetically lower target states
		transitions.push_back(Basin::transition_t(source_energy,      // energy of source state
							  target_basin_index, // index of target basin
							  neigh_energy,       // energy of target state
							  tE)                 // energy of the transition state
				      );

		if (debug_out) {
		    std::cerr <<"\t"<< transitions[transitions.size()-1]
			      << " by ";
		    move->print(std::cerr);
		    std::cerr << std::endl;
		}
	    }
	}

	if (debug_out) std::cerr << "  " << transitions.size() << " transitions, "<< moves_counter << " moves" <<std::endl;
	


	if (moves_counter==0) {
	    if (debug_out) std::cerr << "Ignore frozen state"<<std::endl;
	    return;
	}

	// index of the basin of the processed source state
	size_t source_basin_index;
	
	// ------------------------------------------------------------
	// Perform basin assignment of source_state.
	//
	// Either construct a new basin with local minimum source_state
	// or assign source_state to an existing basin
	//
	if ( min_transition_energy > source_energy ) {
	    // no transition state to an energetically lower target
	    // state is energetically lower than the source state
	    //
	    // Consequently, source_state is a new local minimum
	    
	    source_basin_index = basins.size();
	    
	    if (debug_out) std::cerr << "  New basin "<<source_basin_index<<std::endl;

	    // put state into hash
	    state_hash[source_state.encode()] = source_basin_index;
	    
	    // generate new basin and put into object's basin list
	    Basin new_basin(source_basin_index,source_state,source_energy,model);
	    basins.push_back(new_basin);
	    
	} else {
	    // source_state is not a local minimum but belongs to basin
	    // min_tE_target_basin_index, which is the basin that is reached
	    // with the lowest transition energy.
	    
	    // handle case where basin with index source_basin_index
	    // was merged before
	    source_basin_index = resolve_basin_index(min_tE_target_basin_index);

	    if (debug_out) std::cerr << "  Assign to basin "<<source_basin_index<<std::endl;
	    
	    // assign basin index source_basin_index to source_state and register
	    // source_state as new member of the basin
	    state_hash[source_state.encode()] = source_basin_index;
	    basins[source_basin_index].add_state(source_energy,model);
	}
	
	// ----------------------------------------
	// Register all transitions from source_basin_index to other basins.
	// DON'T Register reverse transitions from target to source
	//
	for(std::vector<Basin::transition_t>::const_iterator it=transitions.begin(); 
	    transitions.end()!=it; ++it) {
	    if (it->target_basin_index != source_basin_index) {
		
		if (debug_out)
		std::cerr << "add transition "<<source_basin_index<<" -> "
			  << it->target_basin_index << " "
			  << it->transition_energy << " "
			  << it->barrier_energy() << " "
			  << std::endl;
		
		basins[source_basin_index].add_transition(*it,model);
		
		//basins[it->target_basin_index].add_transition(it->reverse( source_basin_index ),model);
		
	    }
	}

	Basin &source_basin = basins[source_basin_index];

	// ------------------------------------------------------------
	// Merge basins.  First, merge neighbor basins with barrier
	// height smaller than min_barrier_height to the source basin.
	// Second, merge the source basin to the neighbor basin where
	// transition has the smallest barrier height below the
	// threshold.
	//
	for(std::vector<Basin::transition_t>::const_iterator it=transitions.begin(); 
	    transitions.end()!=it; ++it) {
	    Basin &target_basin=basins[it->target_basin_index];
	    
	    if ((it->target_basin_index > source_basin_index)
		&& ((it->barrier_energy() 
		     - target_basin.get_minimum_energy()) <  min_barrier_height)
		&& (target_basin.get_merged_to()!=source_basin_index)
		) {
		// merge target basin to source basin
		
		if (debug_out) 
		    std::cerr << "  Merge target basin "
			      <<it->target_basin_index
			      <<" to source "<< source_basin_index
			      << " via barrier height "
			      << (it->barrier_energy() 
				  - target_basin.get_minimum_energy())
			      <<std::endl;
		basins[source_basin_index].
		    merge_in_basin(basins[it->target_basin_index],*this);
	    }
	}
	
	// determine the lowest transition barrier from the source basin
	// to some other basin and record basin
	size_t source_merge_to=source_basin_index;
	double source_min_barrier=std::numeric_limits<double>::max();
	
	for(std::vector<Basin::transition_t>::const_iterator
		it=transitions.begin(); transitions.end()!=it; ++it) {
	    Basin &target_basin=basins[it->target_basin_index];
	    if ((it->target_basin_index != source_basin_index) 
		&& (target_basin.get_merged_to()!=source_basin_index)) {
		double barrier = it->barrier_energy();
		if (barrier < source_min_barrier) {
		    source_min_barrier = barrier;
		    source_merge_to = it->target_basin_index;
		}
	    }
	}
	
	if (source_min_barrier - source_basin.get_minimum_energy() 
	    < min_barrier_height) {
	    
	    if (debug_out) 
		std::cerr << "  Merge source basin "
			  << source_basin_index
			  <<" to target "<<source_merge_to
			  << " via barrier height "
			  << (source_min_barrier - source_basin.get_minimum_energy())
			  << std::endl;
	    basins[source_merge_to].
		merge_in_basin(basins[source_basin_index],*this);
	}
	
    }

    /** 
     * @brief Print the list of all basins of the barriers graph 
     * 
     * Output is sorted by basin index
     *
     * @param out output stream
     */
    void
    print_basins(std::ostream &out) const {
	basins[0].print_header(out);
	for(size_t i=0; i<basins.size(); ++i) {
	    if (basins[i].get_merged_to()==i) {
		basins[i].print(out,model);
	    } else {
		// printf("merged to %4lu",basins[i].merged_to); basins[i].print(std::cout,model);
	    }
	}
    }

    /**
       @brief print represenation of barrier graph with weights of basins and transitions
    */
    void print_barrier_graph(std::ostream &out) const {
	for(size_t i=0; i<basins.size(); ++i) {
	    if (basins[i].get_merged_to()==i) {
		basins[i].print_edges(out);
	    }
	}	
    }
    
    ~BarrierGraph() {
    }
    
};

int
main(int argc, char **argv)
{
    gengetopt_args_info args_info;
    
    // get options (call gengetopt command line parser)
    if (cmdline_parser (argc, argv, &args_info) != 0)
	exit(1) ;
    
    if ( args_info.inputs_num != 2 ) {
	std::cerr << "Expect two sequences as input on command line"<<std::endl;
	cmdline_parser_print_help();
	cmdline_parser_free(&args_info);
	exit(1);
    }
    
    std::string seqA=args_info.inputs[0];
    std::string seqB=args_info.inputs[1];

    double min_barrier_height=args_info.minh_arg;

    cmdline_parser_free(&args_info);
    
    // set some global variables for Vienna libRNA
    dangles=2;
    
    // construct barrier graph
    BarrierGraph barriers(seqA,seqB,min_barrier_height);
    
    // print basins of barrier graph
    barriers.print_basins(std::cout);
    
    std::cout << std::endl
	      << std::endl;
    barriers.print_barrier_graph(std::cout);
    
    exit(0);
}
