#ifndef UNPAIRED_PF_HH
#define UNPAIRED_PF_HH

extern "C" {
#include "ViennaRNA/data_structures.h" // defined FLT_OR_DBL
}

#include <LocARNA/matrices.hh>

/**
 * \brief Partition functions of RNA ensembles with unpaired sites
 *
 * of single RNA with up to two unpaired ranges i..j and k..j
 *
 * @note Sequence can be given and indexed 5'->3' (forward) or 3'->5'
 * (reverse) direction. In the forward case, the sequence is indexed
 * 1..length(seq) from 5'->3' as usual. In the latter case, the
 * sequence positions start with 1 at the 3' end and run to
 * length(seq) at the 5' end.
 */
class UnpairedPF {
    
public:
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    //! type of probabilities
    typedef FLT_OR_DBL prob_t;
    
    /**
     * @brief construct with sequence
     *
     * construct with given sequence and compute partition
     * functions where a range i..j is unpaired for all ranges
     * i..j, 1<=i<=j<=len
     * @param seq The RNA sequence
     * @param orientation +1 if seq is in forward orientation 5'->3', -1 otherwise.
     * @param maxsitesize maximum length of unpaired regions
     * @param span base pair span for local folding
     * @param window window size for local folding
     * @param model_double_sites whether to support dependant unpaired regions by conditional probabilities
     * (if false, model sites as independent)
     */
    UnpairedPF(const std::string &seq,
	       int orientation, 
	       size_t maxsitesize,
	       size_t span,
	       size_t window,
	       bool model_double_sites);
    
    /**
     * @brief get probability of unpaired range (i.e. pf divided by Z)
     *
     * @param i left end of range
     * @param j right end of range
     * @return probability p of range i..j is unpaired in object's sequence
     * according to Turner energy model
     * @note sequence positions in 1..len
     * @note -RT ln p equals -RT ln Z_unpaired/Z_total equals the
     * energy difference E_unpaired - E_total, where E_unpaired = -RT
     * ln Z_unpaired and E_total=-RT ln Z_total
     * @see unpaired_prob_conditional()
     */
    prob_t
    unpaired_prob_single(size_t i, size_t j) const;

    /**
     * get probability of unpaired range under condition that
     * another range is unpaired
     * @param i left end of first range
     * @param j right end of first range
     * @param k left end of second range
     * @param l right end of second range
     * @return probability p that the first range is unpaired under the condition
     * that the second range is unpaired
     * @see unpaired_prob_single()
     */
    prob_t
    unpaired_prob_conditional(size_t i, size_t j,size_t k, size_t l) const;

    /**
     * get joint probability that two ranges are unpaired at the same time
     * @param i left end of first range
     * @param j right end of first range
     * @param k left end of second range
     * @param l right end of second range
     * @return joint probability p that both ranges are unpaired
     * @see unpaired_prob_single()
     */
    prob_t
    unpaired_prob_joint(size_t i, size_t j,size_t k, size_t l) const {
	return
	    unpaired_prob_single(k,l)
	    *
	    unpaired_prob_conditional(i,j,k,l);
    }
        
    /** 
     * Scaling factor for partition functions
     * 
     * @return R * temperature
     */
    double
    RT() const {return RT_;}
    
private:
    
    std::string seq_; //!< the RNA sequence

    /**
     * @brief scaling factor for partition functions/boltzmann weights
     * @note has to be compatible to Vienna lib
     */
    double RT_;

    const size_t maxsitesize_;
    
    const size_t span_;
    const size_t window_;

    const bool model_double_sites_;
    

    LocARNA::Matrix<prob_t> Psingle; //!< matrix to hold unpaired probabilities for single ranges
    
    /**
     * matrix to hold conditional unpaired probabilities.
     * Pcond(i,j)(k,l) is the probability of range k..l unpaired under condition range i..j unpaired  
     */
    LocARNA::Matrix<LocARNA::Matrix<prob_t> > Pcond;


    /** 
     * \brief revert an upper triangular matrix
     * 
     * @param size matrix size
     * @param m matrix with element access (i,j)
     */
    template<class M>
    void
    revert_upper_triangle_matrix(size_t size,M &m);
    
    /** 
     * \brief revert the single probabilities matrix Psingle
     * 
     */
    void
    revert_single_probs();

    /** 
     * \brief revert the conditional probabilities matrix Pcond
     * 
     */
    void
    revert_cond_probs();
    
    /**
     * \brief Calculates single range unpaired probabilities.
     *
     * Calculates all unpaired probabilities calling plfold of
     * libRNA
     * \note Sequence seq has to be upper case and must not
     * contain Ts
     */
    void
    computeSingleProbs();
    
    /**
     * \brief Calculates conditional unpaired probabilities for one condition.
     *
     * Calculates unpaired probabilities under condition i..j unpaired calling
     * plfold of libRNA, which was modified to support unpaired constraints
     * \note Sequence seq has to be upper case and must not contain Ts
     * \post matrix Pcond(i,j) contains conditional probablities
     */
    void
    computeCondProbs(size_t i, size_t j);

    /**
     * \brief calculate all conditional unpaired probabilities.
     * \post matrix of matrices Pcond contains all conditional probablities
     */
    void
    computeCondProbs();

    /**
     * @brief Generic calculation of unpaired probabilities
     * @param[out] P matrix to hold probabilities
     * @param structure constraint structure string or NULL
     * @note Given a constraint string the probabilities will be conditioned.
     * Setting structure to NULL will compute unconditional probabilities.
     */
    void
    computeProbsGeneric(LocARNA::Matrix<prob_t> &P, const char *structure);
    
    
    /**
     * total partition function of the sequence seq
     * \note  currently unused
     */
    pf_t
    total_pf() const;
};

#endif
