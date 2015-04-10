#ifndef HYBRID_PF_HH
#define HYBRID_PF_HH

extern "C" {
#include "ViennaRNA/data_structures.h" // defines FLT_OR_DBL
}

#include <LocARNA/matrices.hh>
#include <LocARNA/sparse_matrix.hh>

#include "const_sparse_matrix.hh"

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
     * @param region_start smallest left end of interaction site
     * @param region_end largest right end of interaction site
     *
     * @note sequence B is in reverse orientation
     */
    HybridPF(const std::string &seqA,
	     const std::string &seqB,
	     size_t maxsitesize,
	     size_t maxsitesize_diff,
	     size_t region_startA,
	     size_t region_endA,
	     size_t region_startB,
	     size_t region_endB
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
     * transition state. However there is a simpler definition of
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
    
    //! type of pf matrix slice for common left ends
    typedef ConstSparseMatrix<unsigned short int,FLT_OR_DBL> pf_matrix_slice_t;
    //typedef LocARNA::SparseMatrix<FLT_OR_DBL> pf_matrix_slice_t;

    /**
     * @brief Hybrid partition functions for all sites with given left ends
     *
     * @param i1 left end of sites in sequence 1 
     * @param i2 left end of sites in sequence 1 
     *
     * @return 2D matrix slice
     */
    HybridPF::pf_matrix_slice_t
    hybrid_pfs_common_left_ends(size_t i1,size_t i2) const { return Q_(i1,i2); }

private:
    
   
    //typedef LocARNA::OMatrix<LocARNA::OMatrix<FLT_OR_DBL> > PFOffsetMatrix4D;
    typedef std::pair<size_t,size_t> idx_pair_t;

    //! 4D matrix for holding partition functions.
    typedef LocARNA::OMatrix< pf_matrix_slice_t > PFOffsetMatrix4D;
    
    typedef LocARNA::OMatrix<FLT_OR_DBL> PFOffsetMatrix2D;
    
    /**
     * \brief structure containing Boltzmann weights for energy parameter
     * 
     * scaled to current temperature
     */
    pf_paramT *pf_params_; 

    /**
     * \brief structure containing energy parameter
     * 
     * scaled to current temperature
     */
    paramT *params_;

    const std::string &seqA; //!< sequence A
    const std::string &seqB; //!< sequence B
    
    const size_t maxsitesize_;
    const size_t maxsitesize_diff_;

    size_t region_startA_;
    size_t region_endA_;
    size_t region_startB_;
    size_t region_endB_;
    
    const size_t lenA_; //!< length of seqA
    const size_t lenB_; //!< length of seqB
    
    /**
     * @brief scaling factor for partition functions/boltzmann weights
     * @note has to be compatible to Vienna lib
     */
    double RT_;

    short *SA_; //!< encoded sequence A (type 0)
    short *SA1_; //!< encoded sequence A (type 1)
    short *SB_; //!< encoded reverse sequence B (type 0)
    short *SB1_; //!< encoded reverse sequence B (type 1)
    
    /**
     * 4-dim matrix that holds all partition functions 
     * Q_[i1][i2][j1][j2] is the hybrid partition function 
     * for subsequences seqA[i1..j1] and reverse(seqB)[i2..j2].
     */
    PFOffsetMatrix4D Q_;
    
    
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
	return params_->DuplexInit;
    }

};


#endif
