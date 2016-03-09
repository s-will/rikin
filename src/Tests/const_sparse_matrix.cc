#include "../const_sparse_matrix.hh"
#include <string>
#include <cstdlib>

/** @file some unit tests for the class ConstSparseMatrix 
*/

void
test(bool ok, const std::string &text) {
   std::cerr << (ok?"ok":"fail")
	     << " -- "<<text<<std::endl;
   if (!ok) std::exit(-1);
}

int
main(int argc, char **argv) {
    const size_t SIZE=100;

    LocARNA::SparseMatrix<size_t> m;
    
    size_t c=1;
    for (size_t i=1; i<SIZE; i++) {
	m(i,3*i) = c++;
    }
    
    ConstSparseMatrix<size_t,size_t> csm;
    csm=m;

    bool ok=true;
    c=1;
    for (size_t i=1; i<SIZE; i++) {
	ok &= csm(i,3*i) == c++;
	ok &= csm.exists(i,3*i);
	ok &= csm(i,3*i+1) == 0;
	ok &= csm(i,3*i+2) == 0;
	ok &= !csm.exists(i,3*i+2);
    }
    
    test(ok,"copy and reread");

    test(csm.size()==m.size()?"ok":"fail"
	 ,"matrix size after copy");
    
}

