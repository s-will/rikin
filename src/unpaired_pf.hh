#ifndef UNPAIRED_PF_HH
#define UNPAIRED_PF_HH

extern "C" {
#include "ViennaRNA/data_structures.h" // defined FLT_OR_DBL
}

#include <LocARNA/matrices.hh>

//! \brief Partition functions of RNA ensembles with unpaired sites

//! of single RNA with unpaired
//! ranges i..j
class UnpairedPF {
    
public:
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    //! type of probabilities
    typedef FLT_OR_DBL prob_t;
    
    //! construct with sequence
    //!
    //! construct with given sequence and compute partition
    //! functions where a range i..j is unpaired for all ranges
    //! i..j, 1<=i<=j<=len
    //! @param seq The RNA sequence
    UnpairedPF(const std::string seq);

    //! get probability of unpaired range (i.e. pf divided by Z)
    //! @param i left end of range, position in sequence
    //! @param j right end of range, position in sequence
    //! @return probability p of range i..j is unpaired in object's sequence
    //! according to Turner energy model
    //! @note sequence positions in 1..len
    //! @note -kT ln p equals -kT ln Z_unpaired/Z_total equals the
    //! energy difference E_unpaired - E_total, where E_unpaired = -kT
    //! ln Z_unpaired and E_total=-kT ln Z_total
    //! @see get_unpaired_prob_conditional
    prob_t
    get_unpaired_prob_single(size_t i, size_t j) const;

    //! get probability of unpaired range under condition that
    //! another range is unpaired
    //! @param i left end of first range, position in sequence
    //! @param j right end of first range, position in sequence
    //! @param k left end of second range, position in sequence
    //! @param l right end of second range, position in sequence
    //! @return probability p that the first range is unpaired under the condition
    //! that the second range is unpaired
    //! @see get_unpaired_prob_single
    prob_t
    get_unpaired_prob_conditional(size_t i, size_t j,size_t k, size_t l) const;
    
    // for convenience we consider to add methods that return joint
    // probability and ensemble energy differences corresponding to
    // the various probabilities
    
private:
    
    const std::string seq; //!< the RNA sequence

    LocARNA::Matrix<prob_t> Psingle; //!< matrix to hold unpaired probabilities for single ranges
    
    //! matrix to hold conditional unpaired probabilities.
    //! Pcond(i,j)(k,l) is the probability of range k..l unpaired under condition range i..j unpaired  
    LocARNA::Matrix<LocARNA::Matrix<prob_t> > Pcond;
    
    //! \brief Calculates single range unpaired probabilities.
    //!
    //! Calculates all unpaired probabilities calling plfold of
    //! libRNA
    //! \note Sequence seq has to be upper case and must not
    //! contain Ts
    void
    computeSingleProbs();
    
    //! \brief Calculates conditional unpaired probabilities for one condition.
    //!
    //! Calculates unpaired probabilities under condition i..j unpaired calling
    //! plfold of libRNA, which was modified to support unpaired constraints
    //! \note Sequence seq has to be upper case and must not contain Ts
    //! \post matrix Pcond(i,j) contains conditional probablities
    void
    computeCondProbs(size_t i, size_t j);

    //! \brief calculate all conditional unpaired probabilities.
    //! \post matrix of matrices Pcond contains all conditional probablities
    void
    computeCondProbs();

    //! @brief Generic calculation of unpaired probabilities
    //! @param[out] P matrix to hold probabilities
    //! @param structure constraint structure string or NULL
    //! @note Given a constraint string the probabilities will be conditioned.
    //! Setting structure to NULL will compute unconditional probabilities.
    void
    computeProbsGeneric(LocARNA::Matrix<prob_t> &P, const char *structure);
    
    
    //! total partition function of the sequence seq
    //! \note not used currently
    pf_t
    total_pf() const;
};

#endif
