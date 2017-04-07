#ifndef PAIR_PFS_HH
#define PAIR_PFS_HH

#include "const_sparse_matrix.hh"
#include <LocARNA/sparse_matrix.hh>
#include <algorithm>
#include <vector>
#include <string>
#include <zlib.h>
#include <math.h>

/**
 * @brief Represent probabilities / partition functions of pairs for a set of states
 *
 * Read and write pair partition functions for a sequence of
 * states resp. from and to files; provide access
 *
 * Supported file format: one line per state; each line consists of
 * the state index in the first column, the total partition function
 * in the second column, and the sequences of triples i j p_ij, where
 * p_ij is the probability of pair (i,j) in the state; blank-separated
 *
 * @note We store pair probabilities (in files and in memory); one
 * reason is that probabilities can be reasonably sparsified due to a
 * minimal probability cutoff.
 */
class PairPfs {
public:
    typedef double prob_t; //!< probability type
    typedef double pf_t;   //!< partition function type
    typedef ConstSparseMatrix<size_t,prob_t> const_pp_matrix_t; //!<matrix of pair probabilties
    typedef LocARNA::SparseMatrix<prob_t> pp_matrix_t; //!<matrix of pair probabilties

    //! vector of pairs with probabilities
    typedef std::vector< std::pair<std::pair<size_t,size_t>, prob_t> > pp_vec_t;
private:
    //! vector of matrices of pair probabilities of each state
    std::vector<const_pp_matrix_t> pps_;
    
    //! vector of total partition function of each state
    std::vector<pf_t> Z_;
    
public:
    /**
     * @brief construct from file
     *
     * Read pair probabilities for a series of states from
     * an "pp" flat file. The input file may be gzip-compressed or
     * uncompressed, which is auto-detected.
     */
    PairPfs(const std::string &filename);

    /**
     * @brief read access to pair probabilities of a state
     * @param state_idx index of the state
     * @return matrix of the state's base pair partition functions
     *
     * @note the state's pair partition functions are products of the
     * pair probabilities and the total partition function of the
     * state; @see Z()
     */
    const const_pp_matrix_t &
    pair_probs(size_t state_idx) const {return pps_[state_idx];}

    /**
     * @brief total partition function of a state
     * @param state_idx index of the state
     * @return the state's partition functions
     */    
    const pf_t &
    Z(size_t state_idx) const {return Z_[state_idx];}
    
protected:
    /**
     * @brief read pp matrices from compressed or uncompressed file
     */
    gzFile
    gzread_pps(gzFile fh);
    
public:
    /**
     * @brief write pp matrices from file
     * @param out output stream
     * @param min_prob minimal pair probability for output 
     */
    std::ostream &
    write_pps(std::ostream &out,prob_t min_prob) const;

    /**
     * @brief write pp matrices from compressed file
     * @param fh zlib file handle
     * @param min_prob minimal pair probability for output 
     */
    gzFile
    gzwrite_pps(gzFile fh,prob_t min_prob) const;

    /**
     * @brief Write information of one state to line of pair probability file
     *
     * @param out output stream
     * @param state_idx index of state
     * @param state_pf partition function of state
     * @param pfs iteratable association of partition functions associated with pairs 
     * @param min_prob minimal probability in output
     * @return stream
     */
    template<class PFs_t>
    static
    std::ostream &
    write_state(std::ostream &out,
		size_t state_idx, 
		pf_t state_pf,
		PFs_t pfs, 
		prob_t min_prob,
		bool is_pf=true);

    /**
     * @brief Read information of one state in line of pair probability file
     *
     * @param in input stream
     * @param[out] state_idx index of state
     * @param[out] state_pf partition function of state
     * @param[out] pps pair 
     * @return stream
     */
    static
    std::istream &
    read_state(std::istream &in,
	       size_t &state_idx,
	       pf_t &state_pf,
	       pp_matrix_t &pps);


};

template<class PFs_t>
std::ostream &
PairPfs::write_state(std::ostream &out,
		     size_t state_idx,
		     pf_t state_pf,
		     PFs_t pfs,
		     prob_t min_prob,
		     bool is_pf) {
    
    pp_vec_t vec;
    for (auto &x: pfs) {
	prob_t prob = is_pf ? x.second/state_pf : x.second;
	if ( prob  > min_prob ) {
	    vec.push_back(std::make_pair(x.first,prob));
	}
    }    

    sort(vec.begin(),vec.end());
    
    out << state_idx << "\t";

    out << state_pf << "\t";
    
    for (auto &x: vec) {
	prob_t prob = x.second;
	prob = round( prob / min_prob ) * min_prob ;
	    
	out << x.first.first << " "
	    << x.first.second << " " 
	    << prob << "\t";
    }
    out << std::endl;
    return out;
}

#endif //PAIR_PFS_HH
