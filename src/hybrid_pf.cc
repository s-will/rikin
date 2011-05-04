#include  <string.h> // strchr needed in pair_mat.h ! dependency BUG
#include <assert.h>

#include "hybrid_pf.hh"

extern "C" {
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loop_energies.h"
}

HybridPF::HybridPF(const std::string &seqA_,const std::string &seqB_):
    seqA(seqA_),
    seqB(seqB_),
    lenA(seqA_.length()),
    lenB(seqB_.length()),
    RT_( (temperature+K0)*GASCONST/1000.0 )
{
    // std::cout << "Create HybridPF from sequences "<<seqA<<" and "<<seqB<<std::endl; 
    
    create_temporary();

    compute_hybrid_pf();
    
}

HybridPF::~HybridPF() {
    free_temporary();
}

void
HybridPF::create_temporary() {
    SA = encode_sequence(seqA.c_str(),0);
    SA1 = encode_sequence(seqA.c_str(),1);
    
    SB = encode_sequence(seqB.c_str(),0);
    SB1 = encode_sequence(seqB.c_str(),1);
    
    params =  scale_parameters();
    pf_params = get_scaled_pf_parameters();

    make_pair_matrix();
}

void
HybridPF::free_temporary() {
    free(pf_params);
    free(params);
    free(SA);
    free(SA1);
    free(SB);
    free(SB1);
}

int
HybridPF::pair_type(size_t i1, size_t i2) const {
    assert(1 <= i1);
    assert(i1 <= lenA);
    assert(1 <= i2);
    
    if (! (i2 <= lenB) ) {
	std::cerr <<"i2: "<<i2<<" lenB: "<<lenB<<std::endl;
	assert(i2 <= lenB);
    }

    return pair[SA[i1]][SB[i2]];
}

HybridPF::pf_t
HybridPF::exp_ILoopE(size_t i1, size_t i2, size_t k1,  size_t k2) const {
        
    int ptype_closing = pair_type(i1,i2);
    // note: the enclosed bp type is 'turned around' for lib call
    int ptype_enclosed = rtype[pair_type(k1,k2)];

    if (ptype_closing==0 || ptype_enclosed==0) return 0;

    return
	exp_E_IntLoop(k1-i1-1,k2-i2-1, 
		      ptype_closing,
		      ptype_enclosed,
		      SA[i1+1],
		      SB[i2+1],
		      SA[k1-1],
		      SB[k2-1],
		      pf_params);
}


HybridPF::energy_t
HybridPF::ILoopE(size_t i1, size_t i2, size_t k1,  size_t k2) const {
    
    int ptype_closing = pair_type(i1,i2);
    
    // note: enclosed bp type 'turned around' for lib call
    int ptype_enclosed = rtype[pair_type(k1,k2)]; 
    
    if (ptype_closing==0 || ptype_enclosed==0) return INF/100.0;

    return
	E_IntLoop(k1-i1-1,k2-i2-1, 
		  ptype_closing,
		  ptype_enclosed,
		  SA[i1+1],
		  SB[i2+1],
		  SA[k1-1],
		  SB[k2-1],
		  params)/100.0;
}

void
HybridPF::initialize_hybrid_pf() {
    // resize Q
    Q.resize(lenA+1,lenB+1);
    for (size_t i1=1; i1<=lenA; i1++) {
	for (size_t i2=1; i2<=lenB; i2++) {
	    Q(i1,i2).resize(lenA+1,lenB+1);
	}
    }
 
    // initialisation
    
    // fill with 0
    for (size_t i1=1; i1<=lenA; i1++) {
	for (size_t i2=1; i2<=lenB; i2++) {
	    Q(i1,i2).fill(0);
	} 
    }

    
    // set Q[i][i2][i][i2] to 1
    // if (i,i2) is a possible interaction base pair
    // and 0 otherwise
    for (size_t i1=1; i1<=lenA; i1++) {
	for (size_t i2=1; i2<=lenB; i2++) {    
	    int ptype = pair_type(i1,i2);
	    
	    // Do we want the Duplex Init Energy added to each hybridisation????
	    // This makes energies more similar to cofold results.
	    // Currently, we add no duplex energy (see below)
	    
	    Q(i1,i2)(i1,i2) = 
		(ptype > 0) // is i.k pairing canonical ?
		// ?pf_params->expDuplexInit // instead of 1.0 ?
		?1.0
		:0.0;
	}
    }
}


void
HybridPF::compute_hybrid_pf_common_start(size_t i1, size_t i2) {
    
    // iterate over interaction loops closed by i1.i2 with innner
    // base pair k1.k2 such that k1-i1 + k2-i2 <= MAXLOOP
    
    size_t max_k1 = std::min(lenA, i1+MAXLOOP+1);
    for (size_t k1=i1+1; k1<=max_k1; k1++) {
	
	// here compute min_l1 such that maxlooplength constraint holds
	size_t max_k2 = std::min(lenB, i2+MAXLOOP+1-(k1-i1-1));
	
	for (size_t k2=i2+1; k2<=max_k2; k2++) {
	    
	    // factor for energy contribution of loop
	    // closed by k1.l1 with inner bp k.l
		    
	    pf_t exp_loopE = 
		exp_ILoopE(i1,i2,k1,k2);
	    
	    for (size_t j1=k1; j1<=lenA; j1++) {
		for (size_t j2=k2; j2<=lenB; j2++) {
		    
		    // std::cout << i1 << " " << i2 << " "
		    // 	      << k1 << " " << k2 << " "
		    // 	      << j1 << " " << j2 << " " 
		    // 	      << exp_loopE << " " << Q(k1,k2)(j1,j2)
		    // 	      << std::endl;
		    
		    Q(i1,i2)(j1,j2) += 
			exp_loopE *
			Q(k1,k2)(j1,j2);
		    
		}
	    }
	}
    }
}

void
HybridPF::compute_hybrid_pf() { 
    initialize_hybrid_pf();
    
    for (size_t i1=lenA; i1>=1; i1--) {
	for (size_t i2=lenB; i2>=1; i2--) {
	    if (pair_type(i1,i2)>0) {
		compute_hybrid_pf_common_start(i1,i2);
	    }
	}
    }
}

HybridPF::pf_t
HybridPF::partition_function(size_t i1, size_t j1,size_t i2, size_t j2) const {
    
    // std::cout << "HybridPF::partition_function"
    // 	      << i1 << "-" << j1 << "," << i2 << "-" << j2 << std::endl;
    
    assert(i1<=j1);
    assert(i2<=j2);
    return Q(i1,i2)(j1,j2);
}
