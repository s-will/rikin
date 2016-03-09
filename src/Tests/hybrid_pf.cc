#include "../hybrid_pf.hh"
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
    
    const std::string seqA="AAUGGGGGAAACCCCAUACUCCUCACACACCAAAUCGCCCGAUUUAUCGGGC";
    const std::string seqB="GGGUGUGUAGGAGAGGUUUUACUGAAACAGUGGAAACAAGGAAACACUUGAUUUAGUUAAACCAAGGAAAGACCUCUAUUAGGGGUC";
    size_t maxsitesize=std::max(seqA.length(),seqB.length());
    size_t maxsitesize_diff=10;
    size_t region_startA=1;
    size_t region_endA=seqA.length();
    size_t region_startB=1;
    size_t region_endB=seqB.length();
    
    HybridPF hpf(seqA,seqB,
		 maxsitesize,
		 maxsitesize_diff,
		 region_startA,
		 region_endA,
		 region_startB,
		 region_endB);


    bool ok=true;
    
    for (size_t i1=region_startA; i1<=region_endA; i1++) {
	for (size_t i2=region_startB; i2<=region_endB; i2++) {
	    if (hpf.pair_type(i1,i2)<=0) continue;
	    for (size_t j1=i1+1; j1<=region_endA; j1++) {
		for (size_t j2=i2+1; j2<=region_endB; j2++) {
		    if (hpf.pair_type(j1,j2)<=0) continue;
		    if(abs((int)(j1-i1)-(int)(j2-i2))>maxsitesize_diff) continue;

		    for (size_t k1=i1; k1<=j1; k1++) {
			for (size_t k2=i2; k2<=j2; k2++) {
			    if (hpf.pair_type(k1,k2)<=0) continue;
			    if(abs((int)(k1-i1)-(int)(k2-i2))>maxsitesize_diff) continue;
			    if(abs((int)(j1-k1)-(int)(j2-k2))>maxsitesize_diff) continue;

			    HybridPF::pf_matrix_slice_t
				Qi1i2 = hpf.hybrid_pfs_common_left_ends(i1,i2);

			    HybridPF::pf_matrix_slice_t
				Qk1k2 = hpf.hybrid_pfs_common_left_ends(k1,k2);
			    
			    if(!(
				 Qi1i2(k1,k2)
				 * Qk1k2(j1,j2)
				 <=
				 Qi1i2(j1,j2)
				 )) {
				ok=false;
				std::cerr
				    << i1 << " " << i2 << " "
				    << k1 << " " << k2 << " "
				    << j1 << " " << j2 << " : "
				    << Qi1i2(k1,k2)
				    << "*" << Qk1k2(j1,j2) << "="
				    << ( Qi1i2(k1,k2)
					 * Qk1k2(j1,j2) )
				    << " > "
				    << Qi1i2(j1,j2)
				    << std::endl;
			    }
			}
		    }
		}
	    }
	}
    }
    test(ok,"submultiplicativity of Q");
}
