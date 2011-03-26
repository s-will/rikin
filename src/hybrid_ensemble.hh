#ifndef HYBRID_ENSEMBLE_HH
#define HYBRID_ENSEMBLE_HH

#include <ell/State.hh>
#include <biu/Alphabet.hh>

#include "unpaired_pf.hh"
#include "hybrid_pf.hh"

/*
\file hybrid_ensemble.hh

Idea:

In our model microstates are ensembles of hybridization structures,
which share the same hybridization sites.

An object of HybridEnsembleModel knows the RNA sequences and
computes/stores all ensemble energies, it knows how to describe a
state of the model (zero, one or two 4-tuples, possible extension:
several 4-tuples) and compute the energy of a state.  It defines neighborship
by knowing how to traverse the neighbors of a state. (Maybe
neighborship test?)

In our model, we don't explicitely represent all states, i.e. O(n^8),
but sparsify the 'microstate' space. The HybridEnsembleModel object
needs to know this kind of sparsification and how to enumerate all
states (maybe sorted!?, what do we need for the barrier tree
construction in ELL?).

The model class is used to implement the HybridEnsembleState class,
wich interfaces to the ELL.  An object of this class knows its
description (HybridEnsembleModel::StateDescription) and its model
(HybridEnsembleModel) in order to compute its energy.

*/




/**
 * \brief Class describing the two interacting RNAs
 * that can assign an energy two the description of
 * an hybridization complex
 * ...
 *
 * This class computes pf/probability tables at construction, holds
 * the tables and for energy computation, looks up single energy terms
 * from tables and combines them.
 *
 * @todo Specification, Implementation
 *
 */
class HybridEnsembleModel {
    
public:

    typedef double energy_t;
    
    /**
     * @brief Class for the description of an HybridizationEnsemble in our
     * model
     *
     * An object describes an ensemble of RNA-RNA interaction structures
     * by zero, one ore two hybridization sites. The ensemble comprises all
     * RNA-RNA-interaction structures that are composed from
     * hybridization at each hybridizaztion site (only interaction base
     * pairs, no intramolecular base pairs in each single RNA) and nested
     * intramolecular base pairs outside of hybridization sites.  Each
     * hybridization site is described by a 4-tuple of positions
     * i1,j1,i2,j2 such that, for RNAs A and B, the site consists of
     * subsequence i1..j1 of A and subsequence i2..j2 of B.
     * @todo Implementation
     */
    class StateDescription {
    public:
	
	//! \brief an interaction site
	//!
	//! an interaction site is a four tuple of sequence positions
	//! <pre>
	//! ----\        /-----
	//!     i1------j1     
	//!     i2------j2     
	//!  ---/        \--------
	//! <\pre>
	struct ISite {
	    size_t i1; //!< start position in first sequence 
	    size_t j1; //!< end position in first sequence   
	    size_t i2; //!< start position in second sequence
	    size_t j2; //!< end position in second sequence  
	    
	    // /**
	    //  * Default constructor
	    //  * Initialize with 0
	    //  */
	    // ISite()
	    // 	: i1(0),j1(0),i2(0),j2(0)
	    // {}
	    
	    
	    /** 
	     * \brief construct with four sequence positions
	     * @param i1_ start position in first sequence 
	     * @param j1_ end position in first sequence   
	     * @param i2_ start position in second sequence
	     * @param j2_ end position in second sequence  
	     */
	    ISite(size_t i1_,size_t j1_,size_t i2_,size_t j2_)
		: i1(i1_),j1(j1_),i2(i2_),j2(j2_)
	    {}
	};

	//! type of encoded class representation
	//typedef std::vector<char> code_t;
	typedef std::vector<unsigned char> code_t;   
	
	/** 
	 * @brief Construct without interaction
	 */
	StateDescription();
	
	/** 
	 * @brief construct with single site and ensemble energy
	 * 
	 * @param i1 start position in sequence 1
	 * @param j1 end position in sequence 1
	 * @param i2 start position in sequence 2
	 * @param j2 end position in sequence 2
	 */
	StateDescription(size_t i1, size_t j1, size_t i2, size_t j2);
    
	/** 
	 * @brief construct with two sites and ensemble energy
	 * 
	 * @param i1 start position first site in sequence 1
	 * @param j1 end position first site in sequence 1
	 * @param i2 start position first site in sequence 2
	 * @param j2 end position first site in sequence 2
	 * @param k1 start position second site in sequence 1
	 * @param l1 end position second site in sequence 1
	 * @param k2 start position second site in sequence 2
	 * @param l2 end position second site in sequence 2
	 */
	StateDescription(size_t i1, size_t j1, size_t i2, size_t j2,
			 size_t k1, size_t l1, size_t k2, size_t l2);
	
	//! @brief Get number of hybridization sites
	//! @return number of hybridization sites
	size_t
	num_sites() const;
	
	/** 
	 * Encode class to compressed representation
	 * 
	 * @param[out] code Encoded, compressed representation.
	 */
	void
	encode(code_t &code) const;

	/** 
	 * @brief Decode class from compressed representation
	 * 
	 * @param code Encoded, compressed representation. 
	 */
	void
	decode(const code_t &code);
	

    private:
	//! positions of hybridization sites
	std::vector<ISite> isites;
    };
    
    /**
     * @brief Defines the move set for hybridization ensembles.
     *
     * A move between two hybridization ensembles is either the growth, shrinkage, split,
     * or merge of one or two hybridization sites.
     *
     * The class supports iteration over all neighbors of an hybridization ensemble.
     * @todo Implementation
     *
     * @note How to define neighborship: Two states with the same
     * number of hybridisation sites are \Delta-neighbors if the sum
     * of position differences is smaller or equal \Delta.  Two states
     * with different number of hyb. sites are \Delta,\mu-neighbors,
     * iff 1.) the state with more sites has a site with a fixed
     * minimal size \mu (and the states without that site are
     * \Delta-neighbors) 2.) a \Delta neighbor of the state with more
     * sites has touching site ends, after merging the sites one
     * obtains the second state.
     * 
     */
    class NeighborIterator {
    public:
	
	class ItState {
	    
	    //! @brief Mode of iterator
	    //! 
	    
	    //! Distinguish three modes for the iterator.  In the
	    //! first mode, enumerate shifts of the left and right end
	    //! positions of interaction sites.  Merge of sites is
	    //! implicit, when sites touch each other due to shifts.
	    //! In the second mode, enumerate possible splits of a
	    //! site. In the third mode, enumerate possible new
	    //! sites. (Modes 2 and 3 exist only if this does not
	    //! exceed the maximum number of sites.)
	    enum {SHIFT,SPLIT,NEW_SITE} mode;
	    
	    //! absolute values of deltas for each position of a site
	    //!
	    //! the sign of a delta shift is controlled by sign_permutation
	    std::vector<size_t> xs;
	    
	    //! permutation of signs
	    //!
	    //! The variable encodes signs of shifts by bits, such
	    //! that we can enumerate permutations by counting.  The
	    //! sign of the i-th non-null shift is controlled by the
	    //! i-th bit. We fix that 1 means minus and 0 means plus.
	    size_t sign_permutation;
	    
	    ItState();
	}
	
    private:
	//! origin state
	const StateDescription & origin;
	
    public:
	
	/** 
	 * @brief Constructs with origin 
	 * 
	 * @param origin Description of origin state 
	 */
	NeighborIterator(const StateDescription &origin);
	
	/** 
	 * @brief First iterator state
	 * 
	 * @return first iterator state
	 */
	const ItState &
	firstItState() const;
	
	/** 
	 * @brief Next iterator state
	 * 
	 * @param[in,out] itstate iterator state
	 * @post itstate is next iterator state
	 */
	void
	nextItState(ItState &itstate) const;
	
	/** 
	 * @brief Test for iteration end
	 * 
	 * @return whether iterator state is beyond the last one
	 */
	bool
	isEndItState(const ItState &itstate) const;

	/** 
	 * Apply move corresponding to iterator state to origin
	 * 
	 * @param itstate iterator state 
	 * @param[out] sd state description
	 * @post sd contains description of state after applying itstate to StateDescription origin
	 */
	void
	applyItState(const ItState &itstate, StateDescription &sd)  const;
		
	/** 
	 * Apply move corresponding to iterator state to origin
	 * 
	 * @param itstate iterator state 
	 * 
	 * @return description of state after applying itstate to StateDescription origin
	 */
	StateDescription
	applyItState(const ItState &itstate) const;

    };
    
    //! construct from sequences
    //! @param seqA sequence A
    //! @param seqB sequence B
    HybridEnsembleModel(std::string seqA, std::string seqB);
    
    /** 
     * @brief energy of a state
     * 
     * @param desc Description of a state
     * 
     * @return energy of the state described by desc in the model *this
     *
     * @note O(1) time due to table lookup
     */
    energy_t
    energy(const StateDescription &desc) const;
    
};


/**
 * Represents a ELL-compatible state of an hybrid ensemble model
 *
 */
class HybridEnsembleState : public ell::State {
    
protected:
    
    class NeighborList : public State::NeighborList
    {
    protected:
	const HybridEnsembleModel::NeighborIterator nit;
	
    public:
	
	//! State of an iterator. The information in ItState is used
	//! for computing the next list entry
	class ItState : 
	    public State::NeighborList::ItState,
	    public HybridEnsembleModel::NeighborIterator::ItState
	{
	public:
	    
	    ItState() 
		: State::NeighborList::ItState(),
		  HybridEnsembleModel::NeighborIterator::ItState()
	    {}
	    
	    virtual 
	    ~ItState()
	    {}
	    
	    ItState(const ItState &itstate)
		: State::NeighborList::ItState(itstate),
		  HybridEnsembleModel::NeighborIterator::ItState(itstate)
	    {}

	    ItState(const HybridEnsembleModel::NeighborIterator::ItState &itstate)
		: HybridEnsembleModel::NeighborIterator::ItState(itstate)
	    {}
	    
	    virtual ItState* 
	    clone() const {
		return new ItState(*this);
	    }
	};

	/** 
	 * Construct with origin
	 * 
	 * @param origin State of origin. The neighbor list consists
	 * of the neighbors of this state.
	 * 
	 */
	NeighborList(const HybridEnsembleState* origin) 
	    : nit(origin->description()) 
	{}
	
	virtual ~NeighborList();
	
	/*! Returns a pointer to the first element of the
	 *  virtual list of neighbors of the state of origin.
	 *  Returns NULL if no neighbors exist.
	 */
	virtual HybridEnsembleState* 
	first(State::NeighborList::ItState** itstate) const;
	
	/*! Returns a pointer to the next element of the
	 *  virtual list of neighbors of the state of origin.
	 *  Returns NULL if no more neighbors exist.
	 */
	virtual HybridEnsembleState*
	next(State::NeighborList::ItState* itstate_,
	     State* elem_) const;

    };
    
    // class RandomNeighborList : public State::NeighborList
    // {
    // public:
    // 	class ItState;
	
    // protected:
    // 	//! the state of origin of the neighbors
    // 	const HybridEnsembleState* const origin;
    // public:
    // 	class ItState : public State::NeighborList::ItState
    // 	{
    // 	public:	
    // 	    ItState()
    // 		: lastMove(NULL)
    // 	    {}
	    
    // 	    ItState(const ItState &itstate)
    // 		: lastMove(new MoveSet::Move(*itstate.lastMove))
    // 	    {}
	    
    // 	    virtual 
    // 	    ~ItState() { 
    // 		if (lastMove!=NULL) delete lastMove;
    // 	    }
	    
    // 	    virtual ItState* 
    // 	    clone() const {
    // 		return new ItState(*this);
    // 	    }
    // 	};
	
    // 	RandomNeighborList(const HybridEnsembleState *origin_);
	
    // 	virtual ~RandomNeighborList();
	
    // 	/*! 
    // 	 * @brief  First neighbor state
    // 	 * @return a pointer to a random element of the
    // 	 *  virtual list of neighbors of the state of origin.
    // 	 *  Returns NULL if no neighbors exist.
    // 	 */
    // 	virtual State*
    // 	first(State::NeighborList::ItState** itstate) const;
	
    // 	/*!
    // 	 * @brief  Next neighbor state
    // 	 * @return a pointer to the next random element of the
    // 	 *  virtual list of neighbors of the state of origin.
    // 	 *  Returns NULL if no more neighbors exist.
    // 	 */
    // 	virtual State*
    // 	next(State::NeighborList::ItState* itstate,
    // 	     State* elem) const;
    // };
    
    
private:
    HybridEnsembleModel::StateDescription sd;

public:
    
    const HybridEnsembleModel::StateDescription &
    description() const {
	return sd;
    }
    
    HybridEnsembleState(const HybridEnsembleModel::StateDescription &sd_)
	:sd(sd_)
    {}
    
    HybridEnsembleState()
    {}
    
    virtual ~HybridEnsembleState()
    {}
    
    /*!
     * Access to a State subclass specific ID string to identify 
     * instances of this class.
     * @return the subclass specific ID string
     */
    virtual const std::string& getID( void ) const;
		
    virtual bool
    operator== (const State& state2) const;
    
    virtual bool
    operator!= (const State& state2) const;
		
    /*! Implements a unique order on states based on their energy and
     * string representation. A state is smaller than another one iff
     * it has smaller energy or it has equal energy and a 
     * lexicographically smaller string representation (tie breaker).
     * 
     * This function will be overwritten to achive a better performance
     * than calling getEnergy() and getString().
     * 
     * If the energy function is non-degenerate the string comparison
     * is obsolete.
     * 
     * @param state2 the State object to compare to
     * @return true if this state is smaller than state2 according to 
     *         the unique order of the states
     */
    virtual bool
    operator< (const State& state2) const;

    /*! Comparison function that compares two State pointer based on the
     * less than operator '<' of the first State.
     * The function can be used in STL algorithms, e.g. if only State
     * pointer are stored but the ordering should be based on the 
     * object order. 
     * 
     * @param s1 the State pointer of the object that is asked to be 
     *           smaller (!=NULL)
     * @param s2 the State pointer of the object that is asked to be 
     *           bigger (!=NULL)
     * @return s1->operator <(*s1);
     */
    static
    bool less(const State* s1, const State* s2);
		
    /*! Returns the state specific energy. */
    virtual double
    getEnergy() const;
		
    /*! Returns the minimal number of steps via valid
     *  neighbored states from this to another valid State.
     *  @param state2 the State to reach
     */
    virtual unsigned int
    getMinimalDistance(const State& state2) const;
		
    /*! Returns a pointer to a clone of the current state. 
     * @param toFill a State to make a copy of this, or NULL if a new 
     *        State should be created
     * @return a pointer to a copy of this state, if toFill is != NULL
     *         the return value corresponds to the updated state pointed
     *         to by toFill, otherwise a new object is created !!!
     * */
    virtual State*
    clone( State* toFill = NULL) const;
		
    /*! Returns a new State based on the current state. The new state
     *  differs only by the information given by stringRep. */
    virtual State*
    fromString(const std::string& stringRep) const;
		
    /*! Returns a specific std::string representation of this
     *  State.
     */
    virtual std::string
    toString() const;

    /*! Fill the given string with a specific std::string 
     * representation of this State.
     * @param toFill the string to overwrite
     * @return the changed toFill in parameter
     */
    virtual std::string& toString( std::string & toFill ) const;

    // neighborhood
    /*! Returns a virtual list of all VALID neighbored states in
     *  the energy landscape in a specific order.
     */
    virtual NeighborListPtr
    getNeighborList() const;

    /*! Returns a virtual list of all VALID neighbored states in
     *  the energy landscape in a random order.
     *  If a state is given, this is changed to a neighbor.
     */
    virtual NeighborListPtr
    getRandomNeighborList() const;

    /*! Returns a VALID random neighbored state.
     *  If a state is given, this is changed to a neighbor.
     */
    virtual State* getRandomNeighbor(State* inPlaceNeigh = NULL) const;
		
    //////////////  COMPRESSED REPRESENTATION  ////////////////////////////
		
    /*! Access to a compressed sequence representation of the state.
     * @return the compressed sequence representation
     */
    virtual biu::Alphabet::CSequence
    compress(void) const;
		
    /*! Access to a compressed sequence representation of the state.
     * @param toFill a data structure to write the compressed 
     *               representation too
     * @return the compressed sequence representation
     */
    virtual biu::Alphabet::CSequence&
    compress(biu::Alphabet::CSequence& toFill) const;
		
    /*! Uncompresses a compressed sequencce representation into a new
     * State object.
     * @param cseq the compressed sequence representation of a state
     * @param toFill a state object to uncompress too or NULL if a new 
     *               object has to be created
     * @return new State object that is encoded in cseq or NULL in error
     *         case.
     */
    virtual State* uncompress(const biu::Alphabet::CSequence& cseq, State* toFill ) const;

    /*! Uncompresses a compressed sequencce representation into a this
     * State object.
     * @param cseq the compressed sequence representation of a state
     * @return this or NULL in error case
     */
    virtual State* uncompress(const biu::Alphabet::CSequence& cseq);

};

#endif
