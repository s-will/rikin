#include <string.h> // strchr needed in pair_mat.h ! dependency BUG
#include <assert.h>

#include "hybrid_pf.hh"

// #include "LocARNA/stopwatch.hh"

extern "C" {
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loop_energies.h"
}

HybridPF::HybridPF(const std::string &seqA_,
		   const std::string &seqB_,
		   size_t maxsitesize,
		   size_t maxsitesize_diff,
		   size_t region_startA,
		   size_t region_endA,
		   size_t region_startB,
		   size_t region_endB
		   
		   ):
    seqA(seqA_),
    seqB(seqB_),
    maxsitesize_(maxsitesize),
    maxsitesize_diff_(maxsitesize_diff),
    region_startA_(region_startA),
    region_endA_(region_endA),
    region_startB_(region_startB),
    region_endB_(region_endB),
    lenA_(seqA_.length()),
    lenB_(seqB_.length()),
    RT_( (temperature+K0)*GASCONST/1000.0 )
{
    // std::cout << "Create HybridPF from sequences "<<seqA<<" and "<<seqB<<std::endl; 

    // LocARNA::StopWatch stopwatch(true);
    // stopwatch.start("hybridpf_init");

    create_temporary();

    compute_hybrid_pf();

    // stopwatch.stop("hybridpf_init");
    // stopwatch.print_info(std::cerr);
}

HybridPF::~HybridPF() {
    free_temporary();
}

void
HybridPF::create_temporary() {
    SA_ = encode_sequence(seqA.c_str(),0);
    SA1_ = encode_sequence(seqA.c_str(),1);
    
    SB_ = encode_sequence(seqB.c_str(),0);
    SB1_ = encode_sequence(seqB.c_str(),1);
    
    params_ =  scale_parameters();
    pf_params_ = get_scaled_pf_parameters();

    make_pair_matrix();
}

void
HybridPF::free_temporary() {
    free(pf_params_);
    free(params_);
    free(SA_);
    free(SA1_);
    free(SB_);
    free(SB1_);
}

int
HybridPF::pair_type(size_t i1, size_t i2) const {
    assert(1 <= i1);
    assert(i1 <= lenA_);
    assert(1 <= i2);
    
    if (! (i2 <= lenB_) ) {
	std::cerr <<"i2: "<<i2<<" lenB: "<<lenB_<<std::endl;
	assert(i2 <= lenB_);
    }

    return pair[SA_[i1]][SB_[i2]];
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
		      SA_[i1+1],
		      SB_[i2+1],
		      SA_[k1-1],
		      SB_[k2-1],
		      pf_params_);
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
		  SA_[i1+1],
		  SB_[i2+1],
		  SA_[k1-1],
		  SB_[k2-1],
		  params_)/100.0;
}

void
HybridPF::initialize_hybrid_pf() {
    // resize Q
    Q_.resize(region_endA_-region_startA_+1,
	     region_endB_-region_startB_+1,
	     region_startA_,
	     region_startB_);

    for (size_t i1=region_startA_; i1<=region_endA_; i1++) {
	for (size_t i2=region_startB_; i2<=region_endB_; i2++) {
	    Q_(i1,i2).resize(region_endA_-region_startA_+1,
			    region_endB_-region_startB_+1,
			    region_startA_,
			    region_startB_);
	}
    }
 
    // initialisation
    
    // fill with 0
    for (size_t i1=region_startA_; i1<=region_endA_; i1++) {
	for (size_t i2=region_startB_; i2<=region_endB_; i2++) {
	    Q_(i1,i2).fill(0);
	} 
    }

    
    // set Q[i][i2][i][i2] to
    // 1 if (i,i2) is a possible interaction base pair
    // and 0 otherwise
    for (size_t i1=region_startA_; i1<=region_endA_; i1++) {
	for (size_t i2=region_startB_; i2<=region_endB_; i2++) {
	    int ptype = pair_type(i1,i2);
	    
	    // Do we want the Duplex Init Energy added to each hybridisation????
	    // This makes energies more similar to cofold results.
	    // Currently, we add no duplex energy (see below)
	    
	    Q_(i1,i2)(i1,i2) = 
		(ptype > 0) // is i1.i2 pairing canonical ?
		// ?pf_params_->expDuplexInit // instead of 1.0 ?
		?1.0
		:0.0;
	}
    }
}


void
HybridPF::compute_hybrid_pf_common_start(size_t i1, size_t i2) {
    
    // iterate over interaction loops closed by i1.i2 with innner
    // base pair k1.k2 such that k1-i1 + k2-i2 <= MAXLOOP
    
    size_t max_k1 = std::min(region_endA_, i1+MAXLOOP+1);
    for (size_t k1=i1+1; k1<=max_k1; k1++) {
	
	// here compute min_l1 such that maxlooplength constraint holds
	size_t max_k2 = std::min(region_endB_, i2+MAXLOOP+1-(k1-i1-1));
	
	for (size_t k2=i2+1; k2<=max_k2; k2++) {
	    
	    // factor for energy contribution of loop
	    // closed by i2.i2 with inner bp k1.k2
	    
	    if (pair_type(i1,i2)<=0) continue;
	    
	    pf_t exp_loopE = 
		exp_ILoopE(i1,i2,k1,k2);
	    
	    for (size_t j1=k1; j1<=region_endA_; j1++) {
		
		// site_length_diff = |(j1-i1)-(j2-i2)| <= maxsitesize_diff_
		// ==                 (j1-i1+i2)-maxsitesize_diff_ <= j2 <= (j1-i1+i2)+maxsitesize_diff_
		size_t min_j2=std::max(k2+maxsitesize_diff_,(j1-i1+i2))-maxsitesize_diff_;
		size_t max_j2=std::min(region_endB_,(j1-i1+i2)+maxsitesize_diff_);
		
		for (size_t j2=min_j2; j2<=max_j2; j2++) {
		    
		    // std::cout << i1 << " " << i2 << " "
		    // 	      << k1 << " " << k2 << " "
		    // 	      << j1 << " " << j2 << " " 
		    // 	      << exp_loopE << " " << Q_(k1,k2)(j1,j2)
		    // 	      << std::endl;
		    
		    Q_(i1,i2)(j1,j2) += 
			exp_loopE *
			Q_(k1,k2)(j1,j2);
		    
		}
	    }
	}
    }
}

void
HybridPF::compute_hybrid_pf() { 
    initialize_hybrid_pf();
    
    for (size_t i1=region_endA_; i1>=region_startA_; i1--) {
	for (size_t i2=region_endB_; i2>=region_startB_; i2--) {
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
    return Q_(i1,i2)(j1,j2);
}

HybridPF::pf_t
HybridPF::partition_function_cond(size_t i1, 
				  size_t j1,
				  size_t i2,
				  size_t j2,
				  size_t u1,
				  size_t u2,
				  bool left
				  ) const {
    
    // std::cout << "HybridPF::partition_function_cond"
    // 	      << i1 << "-" << j1 << "," << i2 << "-" << j2 << " | "<< u1 <<"/" << u2 << " unpaired at " << (pos==left?"left":"right") " ends"<< std::endl;
    
    assert(i1+u1<j1);
    assert(i2+u2<j2);
    
    pf_t q=0; // resulting partition function
    
    if (left) {
       	size_t max_k1 = std::min(j1, i1+MAXLOOP+1);
	for (size_t k1=i1+u1+1; k1<=max_k1; k1++) {
	    size_t max_k2 = std::min(j2, i2+MAXLOOP+1-(k1-i1-1));
	    for (size_t k2=i2+u2+1; k2<=max_k2; k2++) {
		q += exp_ILoopE(i1,i2,k1,k2) * Q_(k1,k2)(j1,j2);
	    }
	}
    } else {
	size_t min_k1 = std::max(i1+MAXLOOP+1, j1) - MAXLOOP-1;
	for (size_t k1=min_k1; k1<=j1-u1-1; k1++) {
	    size_t min_k2 = std::max(j2+(MAXLOOP+1+(k1-i1-1)), i2) - (MAXLOOP+1+(k1-i1-1));
	    for (size_t k2=min_k2; k2<=j2-u2-1; k2++) {
		q += Q_(i1,i2)(k1,k2) * exp_ILoopE(k1,k2,j1,j2);
	    }
	}
    }

    return q;
}


