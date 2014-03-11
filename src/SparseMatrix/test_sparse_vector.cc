#include "sparse_vector.hh"

#include <iostream>


template<class index_t,class value_t>
std::ostream &
print(std::ostream &out, const SparseVector<index_t,value_t> &v) {
    for (typename SparseVector<index_t,value_t>::const_iterator it=v.begin();
	 it!=v.end(); ++it)
	out <<it.index()<<":"<<*it <<" ";

    return out;
}

template<class index_t,class value_t, class Function>
SparseVector<index_t,value_t> &
map(SparseVector<index_t,value_t> &v, Function f) {
    for (typename SparseVector<index_t,value_t>::iterator it=v.begin();
	 v.end()!=it; ++it)
	*it = f(*it);
    return v;
}

template<class T>
class mal {
    T x_;
public:
    mal(const T&x):x_(x) {};
    
    T
    operator ()(const T &x) {return x*x_;} 
};

int
main() {
    typedef SparseVector<size_t,double> sparsevector_t;
        
    sparsevector_t v;

    std::cout << v[1] << std::endl;

    v[2]=1.0;
    v[3]=1.5;
    v[5]=2.5;
    v[7]=3.5;

    map(v,mal<double>(4.0));
    
    v[7]*=2;

    print(std::cout,v);
    std::cout << std::endl;
    std::cout << std::endl;
    
    
    typedef SparseVector<size_t,SparseVector<size_t,double> > sparsematrix_t;
    sparsematrix_t m;
    sparsematrix_t::Element e=m[1];
    SparseVector<size_t,double> v2=m[1];
    
    std::cout << v2[2] << std::endl;
    
    
    //std::cout << ((SparseVector<size_t,double>)m[1])[2];
    //std::cout << (m[1])[2];


    // typedef SparseVector<size_t,SparseVector<size_t,double> > sparsematrix_t;
    // sparsematrix_t m;
    
    // sparsematrix_t::Element row=m[2];

    // std::cout << row[3]<<std::endl;    
    // std::cout << m[2][3]<<std::endl;
    

    return 0;

}
