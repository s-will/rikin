#ifndef HYBRID_ENSEMBLE_STATE_HH
#define HYBRID_ENSEMBLE_STATE_HH

#include <ell/State.hh>
#include <biu/Alphabet.hh>

#include "unpaired_pf.hh"
#include "hybrid_pf.hh"


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
	:sd()
    {}
    
    virtual ~HybridEnsembleState();
    
    /*!
     * Access to a State subclass specific ID string to identify 
     * instances of this class.
     * @return the subclass specific ID string
     */
    virtual const std::string& getID() const;
		
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
