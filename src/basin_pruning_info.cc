#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>
#include "basin_pruning_info.hh"

std::ostream &
BasinPruningInfo::write(std::ostream &out,
			double min_contribution, 
			size_t num_input_basins,
			bool sparse) const {
    if (sparse) {
	bool is_first=true;
	for (BasinPruningInfo::contributions_t::const_iterator it=contributions_.begin(); 
	     contributions_.end()!=it; ++it) {
	    if (it->second >= min_contribution) { 
		if (!is_first) {out<< '\t';}
		is_first=false;
		std::ostringstream entry;
		entry << it->first << ":";
		entry.precision(3);
		entry << it->second;
		out << entry.str();
	    }
	}
    } else {

	std::vector<double> row;
	row.resize(num_input_basins);
	
	for (BasinPruningInfo::contributions_t::const_iterator it=contributions_.begin(); 
	     contributions_.end()!=it; ++it) {
	    if (it->second >= min_contribution) { 
		row[it->first] = it->second;
	    }
	}
	
	for (size_t i=0; i<row.size(); i++) {
	    std::ostringstream entry;
	    entry.precision(3);
	    entry << row[i];
	    
	    if (i!=0) out<< '\t';
	    out << entry.str();
	}
    }

    return out;
}

BasinPruningInfo::BasinPruningInfo(size_t index)
    : orig_index_(index),
      contributions_() // no contributions to new basin, other than from this basin itself
{
    contributions_[orig_index_]=1.0;
}

void
BasinPruningInfo::merge_in_basin(size_t orig_index, double contribution, const BasinPruningInfo &bpi, double min_contribution) {
    
    assert(contribution<=1.0);
    assert(orig_index_ != orig_index);

    if (contribution<min_contribution) return;

    // index must not contribute to the basin already
    assert(contributions_.find(orig_index) == contributions_.end());
    
    for ( contributions_t::const_iterator it = bpi.contributions_.begin();
	  bpi.contributions_.end() != it; ++it ) {
	
	if ( contribution * it->second >= min_contribution ) {
	    contributions_[it->first] += contribution * it->second;
	}
    }
}

void BasinPruningInfo::sparsify(double min_contrib) {
    for ( contributions_t::const_iterator it = contributions_.begin();
	  contributions_.end() != it; ) {
        
        if ( it->second < min_contrib ) {
            it=contributions_.erase(it);
        } else {
            ++it;
        }
    }
}

void BasinPruningInfo::clear() {
    contributions_t().swap(contributions_);   // clear x reallocating
}
