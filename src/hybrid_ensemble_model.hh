#ifndef HYBRID_ENSEMBLE_MODEL_HH
#define HYBRID_ENSEMBLE_MODEL_HH

#include <sstream>

#include "unpaired_pf.hh"
#include "hybrid_pf.hh"

#include <algorithm>

#include <math.h>
#include <inttypes.h>

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
 * enumerate all states.
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
    typedef double pf_t;
    
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
	
	typedef std::vector<ISite> isites_t;
	

	//! type of encoded class representation
	// typedef uint64_t code_t;
	typedef std::pair<uint64_t,uint64_t> code_t;   
	
	class code_t_hash {
	public:
	    size_t
	    operator ()( const code_t &x ) const {
		return x.first + x.second*3571;
	    }
	};
	
	/** 
	 * @brief Construct without interaction
	 */
	StateDescription();
	
	/** 
	 * @brief Construct from code
	 */
	StateDescription(const code_t &code);
	
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
	size() const;
	
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
	 * @brief write binary-encoded state to stream
	 */
	std::ostream &
	write_binary(std::ostream &out) const;
	
	/**
	 * @brief read binary-encoded state from istream
	 */
	bool
	read_binary(std::istream &in);
	
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
	
	// support iteration over isites
	isites_t::const_iterator begin() const {return isites.begin();}
	isites_t::const_iterator end() const {return isites.end();}
	isites_t::iterator begin() {return isites.begin();}
	isites_t::iterator end() {return isites.end();}

	
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
	    bool equal = this->size() == sd.size();
	    for (size_t i=0; i<this->size() && equal; i++) {
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
	to_string() const;
	
	/** 
	 * Test validity
	 * 
	 * @param model
	 * 
	 * @return whether object describes a valid state of the model
	 */
	bool
	is_valid(const HybEnsModel &model) const;

	
	/** 
	 * Maximum site length
	 * 
	 * @return largest size of any subsequence in interaction sites 
	 */
	size_t
	max_site_size() const;
	
	/** 
	 * @brief symmetric state for homodimers
	 * @return symmetric state
	 * @note no checks are performed, works only for homodimers
	 */
	StateDescription
	symmetric_state(size_t n) const {
	    StateDescription sd;
	    
	    sd.resize(this->size());
	    
	    for (size_t i=0; i<this->size(); i++) {
		const ISite &os = (*this)[i];
		ISite &ns = sd[this->size()-i-1];
		ns.i1 = n-os.j2+1;
		ns.i2 = n-os.j1+1;
		ns.j1 = n-os.i2+1;
		ns.j2 = n-os.i1+1;
	    }
	    return sd;
	}
		
    private:
	/**
	 * @brief positions of hybridization sites
	 *
	 */
	isites_t isites;
	
    }; // end class StateDescription
	    
    /**
     * @brief convert state description to dot bracket notation
     * @param sd State description
     * @return state as dot bracket string
     */
    std::string
    to_dotbracket(const StateDescription &sd) const;

    
#   include "moves.hh"
    
    /**
     * @brief construct from sequences
     *
     * @param seqA sequence A in 5'->3'
     * @param seqB sequence B in 3'->5'
     * @param region_startA 5' end of region in A
     * @param region_endA   3' end of region in A
     * @param region_startA 3' end of region in B
     * @param region_endA   5' end of region in B
     * @param maxsitesize maximum length of a unpaired site
     * @param model_double_sites whether to model double sites
     */
    HybEnsModel(std::string seqA, std::string seqB,
		size_t maxsitesize,
		size_t maxsitesize_diff,
		size_t region_startA,
		size_t region_endA,
		size_t region_startB,
		size_t region_endB,
		size_t span,
		size_t window,
		bool model_double_sites
		);
    
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
     * @brief Partition function of a state
     * 
     * @param desc Description of a state
     * 
     * @return partition function of the state described by desc in the model *this
     *
     * @note O(1) time due to table lookup
     *
     * @note the partition function equals the Boltzmann weight of
     * energy()
     */
    pf_t
    partition_function(const StateDescription &desc) const {
	return boltzmann_weight(energy(desc));
    }
    
    /**
     * @brief Probability of an interaction in a state
     *
     * @param k1 interaction position in sequence 1
     * @param k2 interaction position in sequence 2
     *
     * @note works only for states with at most one interaction site
     * @todo implement for two interaction sites
     */
    double
    interaction_probability(size_t k1, size_t k2,const StateDescription &sd) const;
    
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
	return hybrid_pf_.ILoopE(i1,i2,j1,j2);
    }

    energy_t
    energy_duplex_init() const {return hybrid_pf_.DuplexInit();}

    
    /**
     * @brief pair type for pair of interacting positions
     * @param i1 position in sequence A
     * @param i2 position in sequence B
     * @return pair type of i1.i2
     */
    int
    pair_type(size_t i1, size_t i2) const {
	return hybrid_pf_.pair_type(i1,i2);
    }
    
    /**
     * @brief Read access to hybrid pf object
     * @return hybrid pf object
     */
    const HybridPF &
    hybrid_pf() const {return hybrid_pf_;}

    // give access to some constants of the model
    
    /** 
     * Size of interaction loops
     * 
     * @return Maximal number of unpaired bases at one side of an interaction loop
     */
    size_t
    maxunpinloop() const {
	return maxunpinloop_;
    }
    
    /** 
     * Minimal size of interaction sites
     * 
     * @return Minimal number of bases that the model allows for sites
     * in each sequence
     */
    size_t
    minsitesize() const {
	return minsitesize_;
    }

    /** 
     * Maximal size of interaction sites
     * 
     * @return Maximal number of bases that the model allows for sites
     * in each sequence
     */
    size_t
    maxsitesize() const {
	return maxsitesize_;
    }
    
    
    /** 
     * Maximum site length difference
     * 
     * @return largest size of any subsequence in interaction sites 
     */
    size_t
    maxsitesize_diff() const {
	return maxsitesize_diff_;
    };

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
    size_t
    minsitedist() const {
	return minsitedist_;
    }
    
    /**
     * @brief Start of region in sequence A
     * @return position
     */
    size_t
    region_startA() const {
	return region_startA_;
    }

    /**
     * @brief End of region in sequence A
     * @return position
     */
    size_t
    region_endA() const {
	return region_endA_;
    }
    
    /**
     * @brief Start of region in sequence B
     * @return position
     */
    size_t
    region_startB() const {
	return region_startB_;
    }

    /**
     * @brief End of region in sequence B
     * @return position
     */
    size_t
    region_endB() const {
	return region_endB_;
    }
    
    double 
    RT() const {
	return hybrid_pf_.RT();
    }

    pf_t
    boltzmann_weight(energy_t energy) const {
	return exp(-energy/RT());
    }

    // static sequence manipulation
    
    //! \brief in place complement of sequences over bases ACGTU
    static
    void 
    complement(std::string &s) {
	static char tab[]={'A','C','G','T','U'};
	static char ctab[]={'U','G','C','A','A'};
    
	for (size_t i=0; i<s.length(); ++i) {
	    for (size_t j=0; j<5; ++j) {
		if (s[i]==tab[j]) {
		    s[i]=ctab[j];
		    break;
		}
	    }
	}
    }

    //! \brief in place reverse of strings
    static
    void
    reverse(std::string &s) {
	std::reverse(s.begin(),s.end());
    }

    //! \brief in place reverse complement of strings
    static
    void
    reverse_complement(std::string &s) {
	reverse(s);
	complement(s);
    }

    //! \brief in place normalization of RNA sequence (upcase, T->U)
    static
    void
    normalize_RNA_sequence(std::string &s) {
	std::transform(s.begin(),s.end(),s.begin(),toupper);
	std::replace(s.begin(),s.end(),'T','U');
    }

    bool
    is_homodimer() const {return homodimer_;}

private:
    const std::string seqA_; //!< sequence A
    const std::string seqB_; //!< sequence B

    const UnpairedPF uppfA_; //!< unpaired pf sequence A
    const UnpairedPF uppfB_; //!< unpaired pf sequence A
    const HybridPF hybrid_pf_; //!< hybrid pf
    
    //! @brief maximal number of unpaired bases in a loop (for each sequence)
    const size_t maxunpinloop_;
    
    //! @brief minimal number of bases in one site (for each sequence)
    const size_t minsitesize_;
    
    //! @brief minimal number of bases between two sites (for each sequence)
    //! 
    //! assume the first sites ends at j and the second starts at k,
    //! we require j+minsitedist_+1 <= k (e.g. 1 2 3 4 5 6 7)
    //! @note in order to make the states represnting disjoint sets of structures,
    //! minsitedist_ should be larger than maxunpinloop_
    const size_t minsitedist_;


    //! @brief maximal number of bases in one site (for each sequence)
    const size_t maxsitesize_;

    //! @brief maximal difference between numbers of bases of sequences of one site
    const size_t maxsitesize_diff_;
        
    const size_t region_startA_;
    const size_t region_endA_;
    const size_t region_startB_;
    const size_t region_endB_;


    //! @brief whether seqA and seqB form homodimers
    //!
    //! true iff seqB is reverse complement of seqA
    const bool homodimer_;


}; // end class HybEnsModel

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription::ISite &isite);

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription &sd);
    


    
#endif
