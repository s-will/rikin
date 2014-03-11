#ifndef SYMMETRIC_SPARSE_MATRIX
#define SYMMETRIC_SPARSE_MATRIX

#include <tr1/unordered_map>
#include <map>

template <class index_t, class value_t>
class SymmetricSparseMatrix {
    
    typedef std::map<index_t,value_t> vector_t;
    typedef std::map<index_t,vector_t> matrix_t;
        
    //! default value
    const value_t def_;
    
    //! matrix
    matrix_t matrix_;
    
public:
    
    /** @brief Proxy class for elements
     */
    class Element {
	SymmetricSparseMatrix &m_;
	const index_t i_;
	const index_t j_;
	
    public:
	
	Element(SymmetricSparseMatrix &m, const index_t &i,const index_t &j)
	    : m_(m), i_(i),j_(j)
	{}

	operator value_t () {
	    typename matrix_t::iterator rowit=m_.matrix_.find(i_);
	    if (m_.matrix_.end()==rowit) return m_.def_;
	    vector_t &row=rowit->second;
	    typename vector_t::iterator it=row.find(j_);
	    if (row.end()==it) return m_.def_;
	    return it->second;
	}
	
	Element &
	operator =(const value_t &v) {
	    m_.matrix_[i_][j_]=v;
	    m_.matrix_[j_][i_]=v;
	    return *this;
	}
	
	Element &
	operator +=(const value_t &v) {
	    value_t res=m_[i_][j_] + v;
	    
	    m_.matrix_[i_][j_] = res;
	    m_.matrix_[j_][i_] = res;
	    
	    return *this;
	}

	Element &
	operator -=(const value_t &v) {
	    value_t res=m_[i_][j_] - v;
	    
	    m_.matrix_[i_][j_] = res;
	    m_.matrix_[j_][i_] = res;
	    
	    return *this;
	}

	Element &
	operator *=(const value_t &v) {
	    value_t res=m_[i_][j_] * v;
	    
	    m_.matrix_[i_][j_] = res;
	    m_.matrix_[j_][i_] = res;
	    
	    return *this;
	}

	Element &
	operator /=(const value_t &v) {
	    value_t res=m_[i_][j_] / v;
	    
	    m_.matrix_[i_][j_] = res;
	    m_.matrix_[j_][i_] = res;
	    
	    return *this;
	}

    };

    typedef Element elem_t;

    class RowIterator;
    
    /** @brief Proxy class for rows/cols
     */
    class Row {
	friend class RowIterator;

	SymmetricSparseMatrix &m_;
	index_t i_;
	
    public:
	Row(SymmetricSparseMatrix &m, const index_t &i)
	    : m_(m), i_(i)
	{}

	const elem_t&
	operator [] (const index_t &j) const { 
	    return Element(m_,i_,j);
	}
	
	elem_t
	operator [] (const index_t &j) { 
	    return Element(m_,i_,j);
	}
	
	const index_t & i() const {return i_;}

	typedef typename vector_t::iterator iterator;
	
	iterator begin() { return m_.matrix_[i_].begin(); }
	
	iterator end() { return m_.matrix_[i_].end(); }
	
    };
    
    class RowIterator : public matrix_t::iterator {
	Row *row_;
    public:
	
	RowIterator(typename matrix_t::iterator it,SymmetricSparseMatrix &m) :
	    matrix_t::iterator(it),row_(new Row(m,it->first)) {
	}
	
	~RowIterator() {delete row_;}
	
	Row &
	operator*() const {
	    const typename matrix_t::iterator &p=*this;
	    row_->i_=p->first;
	    return *row_;
	}

	Row *
	operator->() const {
	    const typename matrix_t::iterator &p=*this;
	    row_->i_=p->first;
	    return row_;
	}

    };
    

    typedef Row row_t;
    
    SymmetricSparseMatrix(const value_t &def=0): 
	def_(def) {
    }
    
    const row_t
    operator [](const index_t &i) const {
	return Row(*this,i);
    }
    
    row_t
    operator [](const index_t &i) {
	return Row(*this,i);
    } 


    typedef RowIterator iterator;	
    
    iterator 
    begin() { 
	return(RowIterator(matrix_.begin(),*this));
    }
    
    iterator
    end() { 
	return RowIterator(matrix_.end(),*this);
    }
    
};

#endif
