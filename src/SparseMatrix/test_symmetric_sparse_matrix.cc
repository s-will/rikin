#include "symmetric_sparse_matrix.hh"

#include <iostream>


int
main() {


    typedef SymmetricSparseMatrix<size_t,double> DoubleMat;
    
    DoubleMat m(1);
    
    m[10][20]=0.5;
    
    std::cout << "0.5=="<<m[20][10]<<std::endl;
    std::cout << "1=="<<m[10][21]<<std::endl;
    
    m[20][10] = m[20][10] * 4;
    m[20][10] += m[10][20];
    
    std::cout << "4.0=="<<m[10][20]<<std::endl;

    m[1][2] += 1;
    
    std::cout << "2=="<<m[2][1]<<std::endl;
    
    m[20][10]/=m[1][2]*2;
    std::cout << "1=="<<m[10][20]<<std::endl;
    
    m[10][1]=7;
    
    m[10][5]=13;
    
    std::cout << "Row 10:";
    for (DoubleMat::row_t::iterator it=m[10].begin();
	 m[10].end()!=it; ++it) {
	std::cout  << " " << it->first<<":"<<it->second;
    }
    std::cout << std::endl;

    std::cout << "Complete Matrix:"<<std::endl;
    for (DoubleMat::iterator it=m.begin();
	 m.end()!=it; ++it) {
	for (DoubleMat::row_t::iterator it2=it->begin();
	     it->end()!=it2; ++it2) {
	    std::cout  << " " << it->i()<< " " << it2->first<<":"<<it2->second;
	}
	std::cout << std::endl;
    }
    std::cout << std::endl;



    typedef SymmetricSparseMatrix<std::string,std::string> StringMat;
    StringMat m2("");
    
    m2["Welt"]["Hallo"]="!";
    
    std::string idx="Hallo";
    
    std::cout << idx;
    for (StringMat::row_t::iterator 
	     it=m2[idx].begin();
	 m2[idx].end()!=it; ++it)
	{
	    std::cout  << " " << it->first<<" "<<it->second;
	}
    std::cout << std::endl;
    
    return 0;
}
