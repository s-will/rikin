#include <map>
#include <iostream>


template <class index_t, class value_t>
class SparseVector {

    typedef std::map<index_t,value_t> vector_t;
    
    //! vector
    vector_t vec_;

    //! default
    static const value_t def_;
    
public:

    /*
      Noch zu implementieren ist ein generischer Mechanismus um
      mehrere Dimensionen zu unterstützen.
      
      Optimal scheint mir, alles über templates zu lösen. D.h. wir
      benutzten z.B. pair<pair<size_t,size_t>,size_t> als index tuple etc.
      
      Die gesamte Element-Klasse ist dann ein template aus dem
      Instanzen für k Dimensionen erzeugt werden.
     */
    

    /** @brief Base of proxy classes (const and non-const) for elements
     */
    template<class SparseVectorRef>
    class ElementBase {
    protected:
	SparseVectorRef v_;
	vector<index_t> is_;
	
	template<class inner_value_t>
	const inner_value_t &
	read_lookup(const SparseVector<index_t,inner_value_t> &v, size_t dimension) const {
	    --dimension;
	    const inner_value_t iv=v[is_[dimension]];
	    
	    if (dimension==0) {
		return iv;
	    } 
	    
	    return read_lookup(iv,dimension);
	}
	
    public:
	ElementBase(SparseVectorRef v, const index_t &i)
	    : v_(v)
	{
	    is_.push_back(i);
	}
	
	
	// never allocate space in lookup
	operator const value_t& () const {
	    typename vector_t::const_iterator it=v_.vec_.find(i_);
	    if (v_.vec_.end()==it) return def_;
	    return it->second;
	}

	// multiple indices require special handling
	//...
	
    };
    

    /** @brief Proxy class for non-constant elements
     */
    class Element: public ElementBase<SparseVector &> {
	typedef ElementBase<SparseVector &> base;

    public:

	Element(SparseVector &v, const index_t &i)
	    :  base(v,i)
	{}

	/** Assigment */
	Element &
	operator =(const value_t &val) {
	    base::v_.vec_[base::i_]=val;
	    return *this;
	}

	// for convenience implement combinations of operation and assignment
	Element &
	operator +=(const value_t &val) {
	    base::v_.vec_[base::i_]+=val;
	    return *this;
	}

	Element &
	operator -=(const value_t &val) {
	    base::v_.vec_[base::i_]-=val;
	    return *this;
	}

	Element &
	operator *=(const value_t &val) {
	    base::v_.vec_[base::i_]*=val;
	    return *this;
	}

	Element &
	operator /=(const value_t &val) {
	    base::v_.vec_[base::i_]/=val;
	    return *this;
	}

	Element &
	operator %=(const value_t &val) {
	    base::v_.vec_[base::i_]%=val;
	    return *this;
	}

    };
    
    /** @brief Proxy class for constant elements
     */
    typedef ElementBase<const SparseVector &> ConstElement;
    
    
    /** @brief Generic iterator class
     */
    template <class SparseVectorRef,class vectoriterator,class Element>
    class Iterator {
	SparseVectorRef v_;
	vectoriterator it_;
	
    public:
	Iterator(SparseVectorRef v,vectoriterator it)
	    : v_(v), it_(it)
	{}
	
	Iterator(const Iterator &it)
	    : v_(it.v_), it_(it.it_)
	{}
	
	Element
	operator *() const {
	    return Element(v_,it_->first);
	}
	
	Iterator
	operator ++() {
	    Iterator tmp(*this);
	    it_++;
	    return tmp;
	}
	
	Iterator &
	operator ++(int) {
	    ++it_;
	    return *this;
	}

	bool
	operator !=(const Iterator &it) const {
	    return it_ != it.it_;
	}
	
	index_t index() const {
	    return it_->first;
	}
    };

    typedef Iterator<SparseVector &,
		      typename vector_t::iterator,
		      Element
		      > iterator;
    typedef Iterator<const SparseVector &,
		      typename vector_t::const_iterator,
		      ConstElement
		      > const_iterator;
    
    SparseVector()
    {}
    
    
    ConstElement
    operator [](const index_t &i) const {
	return ConstElement(*this,i);
    }

    Element
    operator [](const index_t &i) {
	return Element(*this,i);
    }
    
    iterator begin() {
	return iterator(*this,vec_.begin());
    }

    iterator end() {
	return iterator(*this,vec_.end());
    }

    const_iterator begin() const {
	return const_iterator(*this,vec_.begin());
    }

    const_iterator end() const {
	return const_iterator(*this,vec_.end());
    }

};


template <class index_t, class value_t>
const value_t SparseVector<index_t,value_t>::def_ = value_t();
