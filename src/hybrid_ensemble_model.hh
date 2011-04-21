#ifndef HYBRID_ENSEMBLE_MODEL_HH
#define HYBRID_ENSEMBLE_MODEL_HH

#include "unpaired_pf.hh"
#include "hybrid_pf.hh"

/**
 * \file hybrid_ensemble_model.hh
 *
 * Idea:
 *
 * In our model microstates are ensembles of hybridization structures,
 * which share the same hybridization sites.
 * 
 * An object of HybridEnsembleModel knows the RNA sequences and
 * computes/stores all ensemble energies, it knows how to describe a
 * state of the model (zero, one or two 4-tuples, possible extension:
 * several 4-tuples) and compute the energy of a state.  It defines
 * neighborship by knowing how to traverse the neighbors of a
 * state. (Maybe neighborship test?)
 *
 * In our model, we don't explicitely represent all O(n^8) many states
 * but sparsify the 'microstate' space. The HybridEnsembleModel
 * object needs to know this kind of sparsification and how to
 * enumerate all states (maybe sorted!?, what do we need for the
 * barrier tree construction in ELL?).
 *
 * The model class is used to implement the HybridEnsembleState class,
 * wich interfaces to the ELL.  An object of this class knows its
 * description (HybridEnsembleModel::StateDescription) and its model
 * (HybridEnsembleModel) in order to compute its energy.
 *
 */


/**
 * \brief Describes model of two interacting RNAs
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
     * @brief Describes an hybridization ensemble state in our model
     *
     * An object describes an ensemble of RNA-RNA interaction
     * structures by zero, one ore two hybridization sites. The
     * ensemble comprises all RNA-RNA-interaction structures that are
     * composed from hybridization at each hybridizaztion site (only
     * interaction base pairs, no intramolecular base pairs in each
     * single RNA) and nested intramolecular base pairs outside of
     * hybridization sites.  Each hybridization site is described by a
     * 4-tuple of positionsClass for the description of i1,j1,i2,j2
     * such that, for RNAs A and B, the site consists of subsequence
     * i1..j1 of A and subsequence i2..j2 of B.
     *
     * In order to ensure that two-site states
     * i1..j1,k1..l1,i2..j2,k2..l2 are always disjoint from one-site
     * states, the sites of two site states have to be separated by a
     * minimal number $L$ of bases (k1-j1-1)+(k2-j2-1) and all
     * interior loops within an interaction site are smaller.
     * 
     * In an extension of the idea, we allow smaller separation, but
     * then require some pairing that separates the sites. This is
     * done by using the partition function where i1..j1 and k1..l1
     * are unpaired but the whole site i1..l1 is not unpaired (and the
     * same for the second molecule).
     * 
     * @todo Implementation
     */
    class StateDescription {
    public:
	
	//! \brief an interaction site
	//!
	//! an interaction site is a four tuple of sequence positions
	//! 
	//! <pre>
	//! ----\        /-----
	//!     i1------j1     
	//!     i2------j2     
	//!  ---/        \--------
	//! </pre>
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
    
    /** @brief A move by adding/removing an interior loop
     *
     * A move between two hybridization ensembles is either the growth, shrinkage, split,
     * or merge of one or two hybridization sites.
     *
     * @note The implementation of moves and iteration over moves
     * relies on polymorphism. This makes allows defining methods for
     * iteration, move application and computation of transition
     * energy for each type of moves well-arranged in separate
     * methods.  Compared to non-polymorph implementations this causes
     * some overhead, which could, however, be largely reduced by
     * redefining the new and delete operators for moves.
     */
    class Move {
	//! a move knows its iterator
	const MoveIterator & const mi; 
	
    public:
	Move(const MoveIterator &mi_): mi(mi_) {}
	
	virtual ~Move();
	
	/** 
	 * @brief Return first move
	 * 
	 * @return Pointer to first move.
	 */
	virtual
	Move *
	first() = 0;
	
	/** 
	 * @brief return next move
	 * 
	 * @return Pointer to next move.
	 */
	virtual
	Move *
	next() = 0;
	
	/** 
	 * @brief return pointer to new object of next move type
	 * 
	 * @return Pointer to new object of next move type.
	 *
	 * @note The instantiations of this method determine the order of
	 * move types, when iterating through all moves.
	 */
	virtual
	Move *
	nextMoveType() = 0;
	
	/** 
	 * Energy of transition state
	 * 
	 * @param sd origin state
	 * @param m  model
	 * 
	 * @return energy of the transition state when applying move to state s
	 */
	virtual
	energy_t
	transitionEnergy(const StateDescription &sd) const = 0;
	
	/** 
	 * Apply move to a state
	 * 
	 * @param sd origin state
	 */
	virtual
	void
	apply(const StateDescription &sd) const = 0;
    };

    /**
     * @brief Move that grows or shrinks left ends of first site 
     * 
     */
    class GrowShrinkMove0 : public Move {
    public:
	GrowShrinkMove0(MoveIterator &mi_): mi(mi_) {}
    
	virtual
	~GrowShrinkMove0();
	
	Move *
	nextMoveType() const {return new GrowShrinkMove1(mi);}
	
	Move *
	first();
	
	Move *
	next();

	energy_t
	transitionEnergy(const StateDescription &sd) const;

	void
	apply(const StateDescription &sd) const;
    };
    
    //! @brief stopper class for the chain of move classes
    //!
    //! on invokation of first()
    //! immediately deletes itself and returns NULL
    class StopMove : public Move {
    public:
	StopMove(MoveIterator &mi_): mi(mi_) {}
	
	virtual
	~StopMove();

	Move *
	nextMoveType() const {return NULL;}
	
	Move *
	first();
	
	Move *
	next();

	energy_t
	transitionEnergy(const StateDescription &sd) const;

	void
	apply(const StateDescription &sd) const;
    };
    
    
    /**
     * @brief Defines iteration over moves (together with Move).
     *
     * The class supports iteration over all neighbors of an hybridization ensemble.
     * @todo Implementation
     *
     * The transition state is defined by the ensemble that contains
     * the added/removed internal loop. 
     *
     * @todo: Update description according to new ideas from Vancouver meeting
     *
     * @note Each move (but new site introduction) grows or shrinks an
     * interaction site by an internal loop.
     *
     * Introducing new sites is limited to new sites of minimal size.
     *
     * The number of unpaired pairs at each side of an interior loop is restricted.
     *
     * Each state can be described by a vector of coordinates of the
     * removed or added loop or the newly introduced site.
     *
     * For new site creation, the transition state is the interior loop
     * closed by the site ends.
     *
     * <b>Implementation note:</b> In contrast to the ELL virtual list
     * this class does not support polymorphism for the 'iterator
     * state' Move.
     *
     * @see Move
     */
    class MoveIterator {
    public:

    private:
	//! origin state
	const StateDescription &origin;
	
	const HybridEnsembleModel &model;
	
	//! maximal number of unpaired bases at one side of an interior/hybridization loop
	const size_t maxunpinloop;
	
    public:
	
	/** 
	 * @brief Constructs with origin 
	 * 
	 * @param origin Description of origin state 
	 */
	MoveIterator(const StateDescription &origin, const HybridEnsembleModel &model, size_t maxunpinloop);
	
	/** 
	 * @brief First move
	 * 
	 * @return Pointer to move to the first neighbor.
	 *	 
	 * @note firstMove() and nextMove() are designed such that
	 * they can be used in a for loop like
	 * for(MoveIterator::Move *m=mi.firstMove(); m!=NULL; m=mi.nextMove(m)) {}
	 *
	 * @note The method creates a new move object on the
	 * heap. This object is deleted automatically by a call to
	 * nextMove() that returns a null pointer, which happens at
	 * the end of a full traversal through the list of
	 * neighbors. Otherwise it should be deleted explicitely by
	 * disposeMove().
	 */
	Move *
	firstMove() const;
	
	/** 
	 * @brief Next move
	 * 
	 * @param[in,out] m pointer to move
	 * 
	 * @return pointer to next move
	 *
	 * @note destructive on *m
	 *
	 * @see firstMove()
	 */
	Move *
	nextMove(Move *m) const;

	/** 
	 * @brief Dispose a move
	 * 
	 * @param move pointer to move
	 *
	 * Deletes the object referenced by move, unless move is the
	 * null pointer.
	 */
	void
	disposeMove(Move *move) const;		
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


#endif
