#ifndef HYBRID_ENSEMBLE_MODEL_HH
#define HYBRID_ENSEMBLE_MODEL_HH

#include <sstream>

#include "unpaired_pf.hh"
#include "hybrid_pf.hh"

#include <math.h>

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
     * structures with zero, one ore two hybridization sites. The
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
	    
	    bool
	    operator == (const ISite &is) const {
		return
		    ( (*this).i1 == is.i1 ) &&
		    ( (*this).i2 == is.i2 ) &&
		    ( (*this).j1 == is.j1 ) &&
		    ( (*this).j2 == is.j2 )    
		    ;
	    }
	    
	    /**
	     * \brief provide parametrized view on the site ends (read only)
	     */
	    size_t end(size_t seq, bool left) const {
		if (seq==0) {
		    return left?i1:j1; 
		} else {
		    return left?i2:j2; 
		}
	    }

	    /**
	     * \brief provide parametrized view on the site ends
	     */
	    size_t &end(size_t seq, bool left) {
		if (seq==0) {
		    return left?i1:j1; 
		} else {
		    return left?i2:j2; 
		}
	    }

	};
	
	//! type of encoded class representation
	typedef std::string code_t;   
	
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
	 * @brief Encode to compressed representation
	 * 
	 * @return code
	 */
	code_t
	encode() const;
	
	/** 
	 * @brief Encode to compressed representation
	 * 
	 * @param[out] code Encoded, compressed representation.
	 *
	 * @return reference to code
	 */
	code_t &
	encode(code_t &code) const;

	/** 
	 * @brief Decode from compressed representation
	 * 
	 * @param code Encoded, compressed representation.
	 *
	 * @return *this
	 */
	HybEnsModel::StateDescription &
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
	
	
	/** 
	 * @brief Equality operator
	 * 
	 * @param sd state description
	 * 
	 * @return whether state description sd and *this are equal
	 */
	bool
	operator == (const StateDescription &sd) {
	    bool equal = this->num_sites() == sd.num_sites();
	    for (size_t i=0; i<this->num_sites() && equal; i++) {
		equal = equal && (*this)[i] == sd[i];
	    }
	    return equal;
	}
	
	/** 
	 * @brief Convert object to string
	 * 
	 * @return String describing the state 
	 */
	std::string
	toString() const;
	
	/** 
	 * Test validity
	 * 
	 * @param model
	 * 
	 * @return whether object describes a valid state of the model
	 */
	bool
	is_valid(const HybEnsModel &model) const;
	
    private:
	/**
	 * @brief positions of hybridization sites
	 *
	 */
	std::vector<ISite> isites;
	
    }; // end class StateDescription
	    
    
#   include "moves.hh"
    
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

    
    /**
     * @brief pair type for pair of interacting positions
     * @param i1 position in sequence A
     * @param i2 position in sequence B
     * @return pair type of i1.i2
     */
    int
    pair_type(size_t i1, size_t i2) const {
	return hybridpf_.pair_type(i1,i2);
    }


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
    
    double RT() const {
	return hybridpf_.RT();
    }

    double
    boltzmann_weight(energy_t energy) const {
	return exp(-energy/RT());
    }

    
private:
    const std::string seqA_; //!< sequence A
    const std::string seqB_; //!< sequence B

    const UnpairedPF uppfA_; //!< unpaired pf sequence A
    const UnpairedPF uppfB_; //!< unpaired pf sequence A
    const HybridPF hybridpf_; //!< hybrid pf
    
    //! @brief maximal number of unpaired bases in a loop (for each sequence)
    const size_t maxunpinloop_;
    
    //! @brief minimal number of bases in a site (for each sequence)
    const size_t minsitesize_;
    
    //! @brief minimal number of bases between two sites (for each sequence)
    //! 
    //! assume the first sites ends at j and the second starts at k,
    //! we require j+minsitedist_+1 <= k (e.g. 1 2 3 4 5 6 7)
    //! @note in order to make the states represnting disjoint sets of structures,
    //! minsitedist_ should be larger than maxunpinloop_
    const size_t minsitedist_;
   
};

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription::ISite &isite);

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription &sd);


#endif
