#ifndef MOVES_HH
#define MOVES_HH

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
     * @return energy of the transition state when applying move
     */
    virtual
    energy_t
    transitionEnergy() const = 0;
	
    /** 
     * Apply move to a state
     * 
     * @param[in,out] sd origin state
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

/**
 * @brief Move that introduces a new site
 *
 * Defines the common features and methods for all new site moves
 * 
 * @note New sites have at least the size minsitesize and at most
 * loopsize
 */    
class NewSiteMove : public Move {
protected:
    StateDescription::ISite is; //! new site
    
    size_t mini1;
    size_t mini2;
    size_t maxj1;
    size_t maxj2;
public:
    NewSiteMove(const MoveIterator &mi);
    
    virtual
    ~NewSiteMove();
    
    bool
    first();
    
    bool
    next();
    
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


#endif //MOVES_HH
