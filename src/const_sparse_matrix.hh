#ifndef SIMPLE_MAP_HH
#define SIMPLE_MAP_HH

#include <LocARNA/sparse_matrix.hh>
#include "ordered_sparse_matrix.hh"
#include <vector>
#include <algorithm>
//#include <cassert>

/**
 * @brief Space saving replacement for constant sparse matrix of value_t
 * 
 * Initialized with a vector of key value pairs, 
 * provides log time access (after sorting the input.)
 */

template<class _index_t, class _value_t>
class ConstSparseMatrix {
public:
    typedef _value_t value_t;
    typedef _index_t index_t;
    typedef std::pair<index_t,index_t> key_t;
    typedef std::pair<key_t, value_t> key_value_t;
    typedef std::vector<key_value_t> key_value_vec_t;
    typedef typename key_value_vec_t::const_iterator const_iterator;
protected:

    key_value_vec_t vec_;
    value_t def_;
    
    static
    bool
    key_compare(const key_t &x,const key_t &y) {
	return x.first<y.first || 
		       (x.first==y.first && x.second<y.second);
    }
    
    bool
    comp_entry_val(const key_value_t &x,
		   const key_t &y) const {
	return key_compare(x.first,y);
    }

    const_iterator
    binsearch (const_iterator first, const_iterator last, const key_t& key) const
    {
	first = std::lower_bound(first,last,key,
				 [] (const key_value_t &x,const key_t &y) {
				     return key_compare(x.first,y);
				 }
				 );
	
	if (first!=last && key < first->first) {
	    return last;
	}
	return first;
    }

    void
    sort() {
	std::sort(vec_.begin(),vec_.end(),
	     [this] (const key_value_t &x,const key_value_t &y) {
		 return key_compare(x.first,y.first);
	     }
	     );
    }

    const_iterator 
    find(const key_t &key) const {
	return binsearch(vec_.begin(), vec_.end(), key);
    };
    
public:
    ConstSparseMatrix():def_() {}
    
    ConstSparseMatrix(const key_value_vec_t &x) : vec_(x), def_() {
	sort();
    }
    
    /** 
     * @brief assign contents of sparse matrix
     * 
     * @param sm sparse matrix
     * @return this
     */
    ConstSparseMatrix &
    operator = (const LocARNA::SparseMatrix<value_t> &sm) {
	vec_.clear();
	def_ = sm.def();
	for (auto &x : sm) {
	    vec_.push_back(x);
	}
	sort();
	return *this;
    }

    /** 
     * @brief copy contents of ordered sparse matrix
     * 
     * @param sm sparse matrix
     * @return this
     */
    ConstSparseMatrix &
    operator =(const OrderedSparseMatrix<value_t> &sm) {
	vec_.clear();
	def_ = sm.def();
	for (auto &x : sm) {
	    vec_.push_back(x);
	}
    }
    
    bool
    exists(size_t x,size_t y) const {
	return find(key_t(x,y)) != vec_.end();
    }

    const value_t &
    operator () (size_t x,size_t y) const {
	auto it = find(key_t(x,y));
	if (it==vec_.end()) {
	    return def_;
	} else {
	    return it->second;
	}
    }

    const_iterator begin() const {
	return vec_.begin();
    }
    
    const_iterator end() const {
	return vec_.end();
    }
        
    size_t
    size() const {
	return vec_.size();
    }

    bool
    empty() const {
	return vec_.size()==0;
    }

    size_t
    capacity() const {
	return vec_.capacity();
    }
};

#endif //SIMPLE_MAP_HH
