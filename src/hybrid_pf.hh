#ifndef HYBRID_PF_HH
#define HYBRID_PF_HH

extern "C" {
#include "ViennaRNA/data_structures.h" // defines FLT_OR_DBL
}

#include <LocARNA/matrices.hh>

/**
 * \brief Partition functions of hybridizations
 * Computes and maintaines pf for all subsequences
 * of two RNA sequences seqA and seqB.
 *
 * @note sequence B is expected in 3'->5' orientation and is indexed
 * 1..length(seqB) ascending from 3' to 5'
 */ 
class HybridPF {
    
public:
    
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    
    //! type of energy
    typedef double energy_t;
    
    /**
     * @brief construct with sequences A and B
     * @param seqA sequence A (5' -> 3')
     * @param seqB sequence B (3' -> 5')
     * @param maxsitesize maximum sequence length in a site
     * @param maxsitesize_diff maximum sequence length difference in a site
     *
     * @note sequence B is in reverse orientation
     */
    HybridPF(const std::string &seqA,
	     const std::string &seqB,
	     size_t maxsitesize,
	     size_t maxsitesize_diff
	     );
    
    //! Destructor frees librna arrays
    ~HybridPF();
    
    /**
     * @brief get partition funtion of hybridization of subsequences
     * @param i1 start position in sequence A
     * @param j1 end position in sequence A
     * @param i2 start position in sequence B
     * @param j2 end position in sequence B
     * @return partition function of hybridization of seqA[i1..j1] and seqB[i2..j2],
     *          where i1.i2 and j1.j2 interact
     */
    pf_t
    partition_function(size_t i1, size_t j1,size_t i2, size_t j2) const;

    /**
     * @brief get partition function of hybridization of subsequences with unpaired condition
     * @param i1 start position in sequence A
     * @param j1 end position in sequence A
     * @param i2 start position in sequence B
     * @param j2 end position in sequence B
     * @param u1 size of unpaired region in A
     * @param u2 size of unpaired region in B
     * @param left whether left end (true) or right end (false)
     *
     * @note: this could be used to compute the energy of a shift move
     * transition state. Howeverm there is a simpler definition of
     * this transition state that does not require conditional hybrid pfs.
     *
     * @return partition function of hybridization of seqA[i1..j1] and seqB[i2..j2],
     *          where i1.i2 and j1.j2 interact and u bases are unpaired at given position pos
     */
    pf_t
    partition_function_cond(size_t i1, size_t j1,size_t i2, size_t j2, size_t u1, size_t u2, bool left) const;
    
    /** 
     * Scaling factor for partition functions
     * 
     * @return R * temperature
     */
    double
    RT() const {return RT_;}
    
private:
    
    
    //! 4D matrix for holding partition functions.
    typedef LocARNA::Matrix<LocARNA::Matrix<FLT_OR_DBL> > PFMatrix4D;


    /**
     * \brief structure containing Boltzmann weights for energy parameter
     * 
     * scaled to current temperature
     */
    pf_paramT *pf_params; 

    /**
     * \brief structure containing energy parameter
     * 
     * scaled to current temperature
     */
    paramT *params; 

    const std::string &seqA; //!< sequence A
    const std::string &seqB; //!< sequence B
    
    const size_t maxsitesize_;
    const size_t maxsitesize_diff_;

    const size_t lenA; //!< length of seqA
    const size_t lenB; //!< length of seqB
    
    /**
     * @brief scaling factor for partition functions/boltzmann weights
     * @note has to be compatible to Vienna lib
     */
    double RT_;

    short *SA; //!< encoded sequence A (type 0)
    short *SA1; //!< encoded sequence A (type 1)
    short *SB; //!< encoded reverse sequence B (type 0)
    short *SB1; //!< encoded reverse sequence B (type 1)
    
    /**
     * 4-dim matrix that holds all partition functions 
     * Q[i1][i2][j1][j2] is the hybrid partition function 
     * for subsequences seqA[i1..j1] and reverse(seqB)[j1..j2].
     */
    PFMatrix4D Q;
    
    
    //! \brief Compute table Q.
    
    /**
     * Compute hybrid partition functions for
     * all subsequences of seqA and seqB
     */
    void 
    compute_hybrid_pf();
    
    //! intitalize matrix Q
    void
    initialize_hybrid_pf();
    
    /**
     * @brief compute all hybrid partition functions for common left ends
     * @param i1 left end position in seqA
     * @param i2 left end position in seqB
     */
    void
    compute_hybrid_pf_common_start(size_t i1, size_t i2);
        
    /**
     * \brief create temporary data structures
     *
     * that are required in the computation of Q
     */
    void 
    create_temporary();
    
    //! free space of temporary data structures
    void
    free_temporary();


public:
    
    /**
     * @brief pair type for pair of interacting positions
     * @param i1 position in sequence A
     * @param i2 position in sequence B
     * @return pair type of i1.i2
     */
    int
    pair_type(size_t i1, size_t i2) const;

    /**
     * @brief Energy of interaction loop
     * @param i1 3' position in sequence A
     * @param j1 5' position in sequence A
     * @param i2 5' position in sequence B
     * @param j2 3' position in sequence B
     *
     * @return Energy of interaction loop closed by i1.i2 
     *         base pair with interior base pair j1.j2 in kcal/mol
     */
    energy_t
    ILoopE(size_t i1, size_t j1, size_t i2, size_t j2) const;
    
   
    /**
     * @brief Boltzmann weight of interaction loop
     * @param i1 3' position in sequence A
     * @param j1 5' position in sequence A
     * @param i2 5' position in sequence B
     * @param j2 3' position in sequence B
     * @return Boltzmann weight of interaction loop closed by i1.i2 
     *          base pair with interior base pair j1.j2
     */
    pf_t
    exp_ILoopE(size_t i1, size_t j1, size_t i2, size_t j2) const;
    
    /** 
     * Duplex initiation energy
     * 
     * @return duplex init energy
     */
    energy_t
    DuplexInit() const {
	return params->DuplexInit;
    }

};


#endif
