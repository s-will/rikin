#ifndef BASIN_PRUNING_INFO_HH
#define BASIN_PRUNING_INFO_HH

#include <iosfwd>
#include <tr1/unordered_map>
#include <limits>

//! Additional information about a basin
class BasinPruningInfo {
    
    size_t orig_index_;
    
    typedef std::tr1::unordered_map<size_t,double> contributions_t;
    contributions_t contributions_;
    
public:

    /**
     * @brief "empty" constructor
     *
     * @param index original basin index ("input index")
     *
     * No other, dissolved basins contribute to this basin yet
     */
    BasinPruningInfo(size_t index=std::numeric_limits<size_t>::max());

    /** 
     * @brief update info after merging in a basin 
     * 
     * @param orig_index original index of merged in basin
     * @param contribution proportion of the merged in basin's contribution
     * @param bpi pruning info of merged in basin
     * @param min_contribution minimum contribution that is recorded (as single term in summations)
     *
     * @note the merged in basin itself is going to be dissolved
     */    
    void
    merge_in_basin(size_t orig_index, double contribution, 
		   const BasinPruningInfo &bpi,
		   double min_contribution);
    
    // /**
    //  * @brief sparsify
    //  * @param min_contrib minimum contribution
    //  *
    //  * Remove all stored contributions less than min_contribution 
    //  */
    // void
    // sparsify(double min_contrib);
    
    /**
     * @brief access original index
     */
    size_t 
    orig_index() const {return orig_index_;} 
    
    /** @brief write to stream
     * @param out output stream
     * @param min_contribution smallest contribution to be written
     * @param num_input_basins number of input basins
     * @param sparse if true, write sparse
     */
    std::ostream &
    write(std::ostream &out, double min_contribution, size_t num_input_basins, bool sparse) const;

    
};


#endif // BASIN_PRUNING_INFO_HH
