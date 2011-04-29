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
 * An object of HybEnsModel knows the RNA sequences and
 * computes/stores all ensemble energies, it knows how to describe a
 * state of the model (zero, one or two 4-tuples, possible extension:
 * several 4-tuples) and compute the energy of a state.  It defines
 * neighborship by knowing how to traverse the neighbors of a
 * state. (Maybe neighborship test?)
 *
 * In our model, we don't explicitely represent all O(n^8) many states
 * but sparsify the 'microstate' space. The HybEnsModel
 * object needs to know this kind of sparsification and how to
 * enumerate all states (maybe sorted!?, what do we need for the
 * barrier tree construction in ELL?).
 *
 * The model class is used to implement the HybridEnsembleState class,
 * wich interfaces to the ELL.  An object of this class knows its
 * description (HybEnsModel::StateDescription) and its model
 * (HybEnsModel) in order to compute its energy.
 *
 * @todo Splitting and merging is too restrictive, (and probably also
 * new site creation)! Currently, splitting introduces separation
 * always of exactly the same size and merge is only possible from
 * exactly one fixed distance of sites.  Can splitting and merging be
 * assymmetric?  The current strategy may be ok for splitting but
 * seems to be worse for merging.  Requiring an exactly fixed size for
 * new sites, will in practice almost limit new sites to perfect
 * stems. Do we have to relax size of merge/split sites. Fixing the
 * distance is not justified.
 */


/**
 * \brief Describes model of two interacting RNAs
 * 
 * This class computes pf/probability tables at construction, holds
 * the tables and for energy computation, looks up single energy terms
 * from tables and combines them.
 */
class HybEnsModel {
    
public:
    
    // forward ref
    class Move;
    
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
     */
    class StateDescription {
    public:
	
	/**
	 * \brief an interaction site
	 *
	 * an interaction site is a four tuple of sequence positions
	 * 
	 * <pre>
	 * ----\        /-----
	 *     i1------j1     
	 *     i2------j2     
	 *  ---/        \--------
	 * </pre>
	 */
	struct ISite {
	    size_t i1; //!< start position in first sequence 
	    size_t i2; //!< start position in second sequence
	    size_t j1; //!< end position in first sequence   
	    size_t j2; //!< end position in second sequence  
	    
	    /**
	     * @brief Default constructor
	     *
	     * Initialize with 0
	     */
	    ISite()
	     	: i1(0),i2(0),j1(0),j2(0)
	    {}
	    
	    
	    /** 
	     * \brief construct with four sequence positions
	     * @param i1_ start position in first sequence 
	     * @param i2_ start position in second sequence
	     * @param j1_ end position in first sequence   
	     * @param j2_ end position in second sequence  
	     */
	    ISite(size_t i1_,size_t i2_,size_t j1_,size_t j2_)
		: i1(i1_),i2(i2_),j1(j1_),j2(j2_)
	    {
		assert(i1<=j1);
		assert(i2<=j2);
	    }
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
	 * @param i2 start position in sequence 2
	 * @param j1 end position in sequence 1
	 * @param j2 end position in sequence 2
	 */
	StateDescription(size_t i1, size_t i2, size_t j1, size_t j2);
    
	/** 
	 * @brief construct with two sites and ensemble energy
	 * 
	 * @param i1 start position first site in sequence 1
	 * @param i2 start position first site in sequence 2
	 * @param j1 end position first site in sequence 1
	 * @param j2 end position first site in sequence 2
	 * @param k1 start position second site in sequence 1
	 * @param k2 start position second site in sequence 2
	 * @param l1 end position second site in sequence 1
	 * @param l2 end position second site in sequence 2
	 */
	StateDescription(size_t i1, size_t i2, size_t j1, size_t j2,
			 size_t k1, size_t k2, size_t l1, size_t l2);
	
	/**
	 * @brief Get number of hybridization sites
	 * @return number of hybridization sites
	 */
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
	
	
	/** 
	 * @brief Read access to interaction sites
	 * 
	 * @param i index of site
	 * 
	 * @return interation site number i
	 */
	const ISite & operator [](size_t i) const {return isites[i];}
	
	/** 
	 * @brief Write access to interaction sites
	 * 
	 * @param i index of site
	 * 
	 * @return non-const reference to interaction site number i
	 
	 * @note We actually need write access for the moves only
	 * (method apply()). Consider restricting this.
	 */
	ISite & operator [](size_t i) {return isites[i];}
	
	/** 
	 * @brief Set number of interaction sites
	 * 
	 * @param s new number of interaction sites
	 */
	void
	resize(size_t s) {isites.resize(s);}
	
    private:
	/**
	 * @brief positions of hybridization sites
	 *
	 */
	std::vector<ISite> isites;
    };
    
    
    // forward reference
    class MoveIterator;
    
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
    protected:
	//! a move knows its iterator
	const MoveIterator & mi; 
	
    public:
	Move(const MoveIterator &mi_): mi(mi_) {
	}
	
	virtual ~Move();
	
	/** 
	 * @brief Set to first move of object type
	 * 
	 * @return whether there is a first move
	 */
	virtual
	bool
	first() = 0;
	
	/** 
	 * @brief Set to next move of object type
	 * 
	 * @return whether there is a next move
	 */
	virtual
	bool
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
	nextMoveType() const = 0;
	
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
	transitionEnergy() const = 0;
	
	/** 
	 * Apply move to a state
	 * 
	 * @param sd origin state
	 */
	virtual
	void
	apply(StateDescription &sd) const = 0;

	/** 
	 * Print move to stream (for debugging)
	 *  
	 * @param out output stream
	 *
	 * @return output stream
	 *
	 * @note since we want polymorphism, this cannot be realized
	 * via operator <<
	 */
	virtual
	std::ostream &
	print(std::ostream &out) const;
    };
    
    /**
     * @brief Move that grows or shrinks left ends
     *
     * Abstract class implementing common functionality for grow
     * shrink moves
     * 
     */
    class GrowShrinkMove : public Move {
    protected:
	size_t k1; //!< new position for site end in seq 1
	size_t k2; //!< new position for site end in seq 2
	
	size_t i1; //!< original site end position in seq 1
	size_t i2; //!< original site end position in seq 2
	
	size_t min1; //!< min position for site end in seq 1
	size_t min2; //!< min position for site end in seq 2
	size_t max1; //!< max position for site end in seq 1
	size_t max2; //!< max position for site end in seq 2
	
	GrowShrinkMove(const MoveIterator &mi);

	virtual
	~GrowShrinkMove();
	

	/** 
	 * @brief generic first grow shrink move
	 * 
	 * @param i1_ position of moved site end in seq 1
	 * @param i2_ position of moved site end in seq 2
	 * @param min1_ minimal position of new site end seq 1
	 * @param min2_ minimal position of new site end seq 2
	 * @param max1_ maximal position of new site end seq 1
	 * @param max2_ maximal position of new site end seq 2
	 * 
	 * @return whether there is a first move
	 */
	bool
	first(size_t i1_,size_t i2_,
	      size_t min1_,size_t min2_,
	      size_t max1_,size_t max2_);
	
	/** 
	 * @brief generic next grow shrink move
	 * 
	 * @return whether there is a next move
	 */
	bool
	next();
	
	/** 
	 * @brief generic first grow shrink move to left
	 * 
	 * @return whether there is a first move
	 */
	bool
	firstLeft();
	
	/** 
	 * @brief generic next grow shrink move to left
	 * 
	 * @return whether there is a next move
	 */
	bool
	nextLeft();
	
	/** 
	 * @brief generic first grow shrink move to right
	 * 
	 * @return whether there is a first move
	 */
	bool
	firstRight();

	/** 
	 * @brief generic next grow shrink move to right
	 * 
	 * @return whether there is a next move
	 */
	bool
	nextRight();
	
	/** 
	 * @brief generic transition state energy calculation
	 * 
	 * @param sd_small state description
	 * @param sd_large state description
	 * @param loop_i1 left position of loop sequence 1
	 * @param loop_i2 left position of loop sequence 2
	 * @param loop_j1 right tposition of loop sequence 1
	 * @param loop_j2 right position of loop sequence 1
	 * 
	 * @return energy of transition state where the loop is attached to
	 * the smaller state in order to yield the larger state
	 */	
	energy_t transitionEnergy(const StateDescription &sd_small,
				  const StateDescription &sd_large,
				  size_t loop_i1, size_t loop_i2, size_t loop_j1, size_t loop_j2
				  ) const;
	
	/** 
	 * @brief generic transition state energy calculation
	 * 
	 * @param site modified interaction site
	 * @param left whether left end is modified by move
	 * 
	 * @return transition energy where the loop by i1.i2, k1.k2 is attached to 
	 * the specified end of specified site
	 */
	energy_t
	transitionEnergy(size_t site, size_t left) const;	

	/** 
	 * Print move
	 * 
	 * @param out output stream
	 * 
	 * @return output stream
	 */
	virtual
	std::ostream &
	print(std::ostream &out) const;
    };
    
    /**
     * @brief Move that grows or shrinks left ends of first site 
     * 
     */
    class GrowShrinkMoveFL : public GrowShrinkMove {
    public:
	GrowShrinkMoveFL(const MoveIterator &mi);
    
	virtual
	~GrowShrinkMoveFL();
	
	Move *
	nextMoveType() const;
	
	bool
	first();
	
	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };

    /**
     * @brief Move that grows or shrinks right ends of first site 
     * 
     */
    class GrowShrinkMoveFR : public GrowShrinkMove {
    public:
	GrowShrinkMoveFR(const MoveIterator &mi);
    
	virtual
	~GrowShrinkMoveFR();
	
	Move *
	nextMoveType() const;
	
	bool
	first();
	
	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };

    /**
     * @brief Move that grows or shrinks left ends of second site 
     * 
     */
    class GrowShrinkMoveSL : public GrowShrinkMove {
    public:
	GrowShrinkMoveSL(const MoveIterator &mi);
    
	virtual
	~GrowShrinkMoveSL();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };
 
    /**
     * @brief Move that grows or shrinks right ends of second site 
     * 
     */
    class GrowShrinkMoveSR : public GrowShrinkMove {
    public:
	GrowShrinkMoveSR(const MoveIterator &mi);
    
	virtual
	~GrowShrinkMoveSR();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };


    /**
     * @brief Move that removes an interaction site
     *
     * @note A site can be removed only if the site is not larger than
     * one interaction loop (checks HybEnsModel::maxunpinloop()) 
     *
     * @note the tranition state for these moves is defined by
     * changing the to be removed site into one interaction loop
     * closed by the end base pairs of the site.
     */
    class RemoveSiteMove : public Move {
	size_t site; //!< number of site that is removed
    public:
	RemoveSiteMove(const MoveIterator &mi);
    
	virtual
	~RemoveSiteMove();

	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;

	std::ostream &
	print(std::ostream &out) const; 

    private:
	bool
	thisornext();
    };
    
    class NewSiteMove : public Move {
    protected:
	size_t i1; //!< left end of new site in seq 1
	size_t i2; //!< left end of new site in seq 2
    public:
	NewSiteMove(const MoveIterator &mi);
	
	virtual
	~NewSiteMove();
	
	std::ostream &
	print(std::ostream &out) const;
    };


    /**
     * @brief Move that introduces the first interaction site
     *
     * @note New sites have exactly the minimal size of an interaction site (HybEnsModel::minsitesize()) 
     *
     * @note the transition state for these moves is defined by
     * changing the new site into one interaction loop
     * closed by the end base pairs of the site.
     */
    class NewSiteMoveF : public NewSiteMove {
    public:
	NewSiteMoveF(const MoveIterator &mi);
    
	virtual
	~NewSiteMoveF();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };
    

    /**
     * @brief Move that introduces a new interaction site at the left
     *
     * @note New sites have exactly the minimal size of an interaction site (HybEnsModel::minsitesize()) 
     *
     * @note the transition state for these moves is defined by
     * changing the new site into one interaction loop
     * closed by the end base pairs of the site.
     */
    class NewSiteMoveL : public NewSiteMove {
    public:
	NewSiteMoveL(const MoveIterator &mi);
    
	virtual
	~NewSiteMoveL();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };
    
    /**
     * @brief Move that introduces a new interaction site at the right
     *
     * @note New sites have exactly the minimal size of an interaction
     * site (HybEnsModel::minsitesize())
     *
     * @note the transition state for these moves is defined by
     * changing the new site into one interaction loop
     * closed by the end base pairs of the site.
     */
    class NewSiteMoveR : public NewSiteMove {
    public:
	NewSiteMoveR(const MoveIterator &mi);
    
	virtual
	~NewSiteMoveR();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;
    };
    

    /**
     * @brief Move that merges two interaction sites
     *
     * Merging of two sites is allowed iff the site distance is model.minsitedist()
     *
     * @note the transition state for these moves is defined by
     * conditioning the ensemble energy of the new single site by the
     * interaction base pairs at the old site ends
     *
     */
    class MergeMove : public Move {
    public:
	MergeMove(const MoveIterator &mi);
    
	virtual
	~MergeMove();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;

	/** 
	 * Print move
	 * 
	 * @param out output stream
	 * 
	 * @return output stream
	 */
	virtual
	std::ostream &
	print(std::ostream &out) const;

    };
    
    /**
     * @brief Move that splits two interaction sites
     *
     * Splitting of a site separates the new sites by exactly model.minsitedist()
     *
     * @note the transition state for these moves is defined by
     * conditioning the ensemble energy of the single site by the
     * interaction base pairs at the new site ends
     *
     */
    class SplitMove : public Move {
	size_t i1; //!< left end of split site in seq 1
	size_t i2; //!< left end of split site in seq 2
    public:
	SplitMove(const MoveIterator &mi);
    
	virtual
	~SplitMove();
	
	Move *
	nextMoveType() const;
	
	bool
	first();

	bool
	next();

	energy_t
	transitionEnergy() const;

	void
	apply(StateDescription &sd) const;

	/** 
	 * Print move
	 * 
	 * @param out output stream
	 * 
	 * @return output stream
	 */
	virtual
	std::ostream &
	print(std::ostream &out) const;
    };
        
    
    /**
     * @brief Defines iteration over moves (together with Move).
     *
     * The class supports iteration over all neighbors of an hybridization ensemble.
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
	const StateDescription &origin_;
	
	const HybEnsModel &model_;
		
    public:
	
	/** 
	 * @brief Constructs with origin 
	 * 
	 * @param origin Description of origin state 
	 */
	MoveIterator(const StateDescription &origin, const HybEnsModel &model);
	
	const StateDescription &
	origin() const {
	    return origin_;
	}

	const HybEnsModel &
	model() const {
	    return model_;
	}
		
	
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
    
    
    /**
     * construct from sequences
     * @param seqA sequence A
     * @param seqB sequence B
     */
    HybEnsModel(std::string seqA, std::string seqB);
    
    /** 
     * @brief Energy of a state
     * 
     * @param desc Description of a state
     * 
     * @return energy of the state described by desc in the model *this
     *
     * @note O(1) time due to table lookup
     */
    energy_t
    energy(const StateDescription &desc) const;
  
    /** 
     * @brief Read sequence A
     * 
     * @return sequence A
     */
    const std::string &
    seqA() const {
	return seqA_;
    }
    
    /** 
     * @brief Read sequence B
     * 
     * @return sequence B
     */
    const std::string &
    seqB() const {
	return seqB_;
    }

    /** 
     * Energy of a hybridization at one interaction site
     * 
     * @param i1 left end of interaction site in seq 1
     * @param i2 left end of interaction site in seq 2
     * @param j1 right end of interaction site in seq 1
     * @param j2 right end of interaction site in seq 2
     * 
     * @return ensemble energy of the hybridization
     */
    energy_t
    energy_hybrid(size_t i1,size_t i2,size_t j1,size_t j2) const;

    /** 
     * Energy of a hybridization at one interaction site
     * 
     * @param is interaction site 
     * 
     * @return ensemble energy of the hybridization
     */
    energy_t
    energy_hybrid(const StateDescription::ISite &is) const;
    
    /** 
     * Energy of unpairing of one interaction site
     * 
     * @param is interaction site
     *
     * @return energy difference of the unpairing
     */
    energy_t
    energy_unpair(const StateDescription::ISite &is) const;
    
    /** 
     * Energy of unpairing of two interaction sites
     * 
     * @param is1 interaction site 1
     * @param is2 interaction site 1
     *
     * @return energy difference of the unpairing
     */
    energy_t
    energy_unpair(const StateDescription::ISite &is1,
		  const StateDescription::ISite &is2) const;
    
    /** 
     * Energy of one hybridisation loop
     * 
     * @param i1 left base pair, position in seq 1
     * @param i2 left base pair, position in seq 2
     * @param j1 right base pair, position in seq 1
     * @param j2 right base pair, position in seq 2
     * 
     * @return energy of hybridisation loop closed by i1.i2 and enclosing j1.j2
     */
    energy_t
    energy_hybrid_loop(size_t i1,size_t i2,size_t j1,size_t j2) const {
	return hybridpf_.ILoopE(i1,i2,j1,j2);
    }

    energy_t
    energy_duplex_init() const {return hybridpf_.DuplexInit();}


public:
    // give access to some constants of the model
    
    /** 
     * Size of interaction loops
     * 
     * @return Maximal number of unpaired bases at one side of an interaction loop
     */
    const size_t
    maxunpinloop() const {
	return maxunpinloop_;
    }
    
    /** 
     * Minimal size of interaction sites
     * 
     * @return Minimal number of bases of a site in each sequence
     */
    const size_t
    minsitesize() const {
	return minsitesize_;
    }
    
    /** 
     * Minimal separation of two sites
     *  
     * @return Minimal number of bases that seperate two sites in each sequence 
     *
     * @note in our model, we explicitely set such a speparation
     * distance.  In order to make ensembles with two sites disjoint
     * from ensembles with one site, it will suffice to make this
     * distance larger than maxunpinloop_
     */
    const size_t
    minsitedist() const {
	return minsitedist_;
    }


    
private:
    const std::string &seqA_; //!< sequence A
    const std::string &seqB_; //!< sequence B

    const UnpairedPF uppfA_; //!< unpaired pf sequence A
    const UnpairedPF uppfB_; //!< unpaired pf sequence A
    const HybridPF hybridpf_; //!< hybrid pf
    
    const size_t maxunpinloop_;
    
    const size_t minsitesize_;
    
    const size_t minsitedist_;
    
};

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription::ISite &isite);

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription &sd);


#endif
