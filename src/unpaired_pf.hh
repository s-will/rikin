#ifndef UNPAIRED_PF_HH
#define UNPAIRED_PF_HH

extern "C" {
#include "ViennaRNA/data_structures.h" // defined FLT_OR_DBL
}

#include <LocARNA/matrices.hh>

//! class for partition functions of single RNA with unpaired
//! ranges i..j
class UnpairedPF {
    
public:
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    
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
    pf_t
    get_unpaired_prob_single(size_t i, size_t j) const;

    //! get probability of unpaired range under condition that
    //! another range is unpaired
    //! @param i left end of first range, position in sequence
    //! @param j right end of first range, position in sequence
    //! @param k left end of second range, position in sequence
    //! @param l right end of second range, position in sequence
    //! @return probability p of first range unpaired under condition
    //! that second range is unpaired
    //! @see get_unpaired_prob_single
    pf_t
    get_unpaired_prob_conditional(size_t i, size_t j,size_t k, size_t l) const;
    
    // for convenience we consider to add methods that return joint
    // probability and ensemble energy differences corresponding to
    // the various probabilities
    
private:
    
    const std::string seq; //!< the RNA sequence
    const size_t length; //!< length of the RNA sequence

    LocARNA::Matrix<pf_t> Q; //!< matrix to hold partition functions
    
    //! calculate all partition functions calling plfold subrouting of
    //! libRNA
    //! \note Sequence seq has to be upper case and must not contain Ts
    void
    computePFs();
    
    //! total partition function of the sequence seq
    //! \note UNUSED
    pf_t
    total_pf() const;
};

#endif
