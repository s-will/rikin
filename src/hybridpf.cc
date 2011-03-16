#include  <stdlib.h>
#include  <string.h>
#include  <stdio.h>
#include  <math.h>
#include <assert.h>

extern "C" {
#include "ViennaRNA/utils.h"
#include "ViennaRNA/fold_vars.h"

#include "ViennaRNA/data_structures.h"
#include "ViennaRNA/energy_par.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params.h"
#include "ViennaRNA/loop_energies.h"

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/LPfold.h>


}

/*
  Open questions: 
  was wollen wir als energie eines Zustandes definieren?
  
  Ensemble Energie Hyb i1..j1,i2..j2 + Delta E unpaired i1..j1 + Delta E unpaired i2..j2
  
 */


// #ifdef _OPENMP
// #include <omp.h>
// #endif


#include <LocARNA/matrices.hh>


//! class for partition functions of single RNA with unpaired
//! ranges i..j
class UnpairedPF {
    
public:
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    
    //! construct with sequence
    //!
    //! construct with given sequence and compute partition
    //! functions where a range i..j is unpaired for all ranges
    //! i..j, 1<=i<=j<=len
    //! @param seq The RNA sequence
    UnpairedPF(const std::string seq);

    //! get partition function with unpaired range divided by Z
    //! @param i left end of range, position in sequence
    //! @param j right end of range, position in sequence
    //! @returns partition function of object's RNA sequence,
    //! where range i..j is unpaired divided by total partition function Z
    //! @note sequence positions in 1..len 
    pf_t 
    get_pf(size_t i, size_t j) const;
    
private:
    
    const std::string seq; //!< the RNA sequence
    const size_t length; //!< length of the RNA sequence

    LocARNA::Matrix<pf_t> Q; //!< matrix to hold partition functions
    
    //! calculate all partition functions calling plfold subrouting of
    //! libRNA
    //! \note Sequence seq has to be upper case and must not contain Ts
    void
    computePFs();

    //! total partition function of the sequence seq
    //! \note UNUSED
    pf_t
    total_pf() const;    
};

UnpairedPF::UnpairedPF(const std::string seq_) 
    : seq(seq_),
      length(seq.size())
{
    computePFs();
}

UnpairedPF::pf_t 
UnpairedPF::get_pf(size_t i, size_t j) const {
    return Q(i,j);
}


UnpairedPF::pf_t
UnpairedPF::total_pf() const {
    const char *sequence=seq.c_str();
    //char structure[length+2];
    
    /* for longer sequences one should also set a scaling factor for
       partition function folding */
    // double e=fold(sequence,structure);
    // free_arrays();
    
    // double kT = (temperature+K0)*GASCONST/1000.;  /* kT in kcal/mol */
    // pf_scale = exp(-e/kT/length);

    pf_scale = -1;
    
    //update_pf_params(seq.length());
    
    // calculate Gibb's free energy in kcal/mol
    double G = pf_fold(sequence,NULL);

    free_pf_arrays();
    
    return exp( - G / (GASCONST/1000.0*(temperature+K0)) );
}


void
UnpairedPF::computePFs() {
    
    // OBSOLETE:
    // we need the total partition function since we want
    // to convert unpaired probabilities
    // to partition functions.
    // Therefore, calculate partition function first.
    // NOTE: it appears there is no way to get this information form the
    // pfl_fold() function of the RNAlib
    // pf_t Z = total_pf();
    

    int plfW=length; // parameter -W of plfold
    int plfL=length; // parameter -L of plfold
    
    // maximal size of unpaired region for which pf is computed
    int maxUnpairedRegionSize=length;
    
    const char *sequence=seq.c_str();
    
    float cutoff = 0.0;    /* bpcutoff for plfold, does not influence
			      our results */
    
    // declare and alloc pup array
    double ** pup = new double *[length+1];
    
    pup[0] = new double[1]; /* we need only entry 0 */
    pup[0][0] = (double) maxUnpairedRegionSize;
    
    // don't compute conditional probabilities
    plist *dpp=NULL;
    
    pf_scale = -1;
    
    //! call libRNA
    //update_pf_paramsLP(length); //? update_fold_params() was working here
    pfl_fold(const_cast<char *>(sequence), plfW, plfL, cutoff, pup, &dpp, NULL, NULL);
    
    // copy result to matrix Q
    Q.resize(length+1,length+1);
    Q.fill(0.0);
    // i+u: end position of unpaired region counting from 0, u: length of unpaired region
    for(int u=1; u<=maxUnpairedRegionSize; u++) {
	for(size_t i=0; i<length-u+1; i++) {
	    Q(i,i+u-1)=pup[i+u][u]; // this is actually PF/Z !!!
	}
    }
    
    // free pup array
    for (size_t i=0;i<=length;i++) {delete pup[i];}
    delete pup;
}



//! class for hybrid partition functions of all subsequences
//! of two RNA sequences seqA and seqB
class HybridPF {
    
public:
    
    //! type of partition functions
    typedef FLT_OR_DBL pf_t;
    
    //! construct with sequences A and B
    //! @param seqA sequence A (3' - 5')
    //! @param seqB sequence B (3' - 5')
    HybridPF(const std::string &seqA,const std::string &seqB);
    
    //! get partition funtion of hybridization of subsequences
    //! @params i1 start position in sequence A
    //! @params j1 end position in sequence A
    //! @params i2 start position in sequence B
    //! @params j2 end position in sequence B
    //! @returns partition function of hybridization of seqA[i1..j1] and seqB[i2..j2],
    //!          where i1.i2 and j1.j2 interact
    pf_t
    get_pf(size_t i1, size_t j1,size_t i2, size_t j2) const;
    


private:

    //! 4D matrix for holding partition functions
    //! replace by triangle matrix for space/performance fetishism
    typedef LocARNA::Matrix<LocARNA::Matrix<FLT_OR_DBL> > PFMatrix4D;


    //! structure containing Boltzmann weights for energy parameter
    //! scaled to current temperature
    pf_paramT *pf_params; 

    const std::string &seqA; //!< sequence A
    const std::string &seqB; //!< sequence B
    
    const size_t lenA; //!< length of seqA
    const size_t lenB; //!< length of seqB

    short *SA; //!< encoded sequence A (type 0)
    short *SA1; //!< encoded sequence A (type 1)
    short *SB; //!< encoded reverse sequence B (type 0)
    short *SB1; //!< encoded reverse sequence B (type 1)
    
    //! 4-dim matrix that holds all partition functions 
    //! Q[i1][i2][j1][j2] is the hybrid partition function 
    //! for subsequences seqA[i1..j1] and reverse(seqB)[j1..j2].
    PFMatrix4D Q;
    
    //! compute table Q with hybrid partition functions for
    //! all subsequences of seqA and seqB
    void 
    compute_hybrid_pf();
    
    //! intitalize matrix Q
    void
    initialize_hybrid_pf();
    
    //! compute all hybrid partition functions for common left ends
    //! @param i1 left end position in seqA
    //! @param j1 left end position in seqB
    void
    compute_hybrid_pf_common_start(size_t i1, size_t i2);
        
    //! create temporary data structures
    //! that are required in the computation of Q
    void 
    create_temporary();
    
    //! free space of temporary data structures
    void
    free_temporary();
    
    //! Boltzmann weight of interaction loop
    //! @params i1 3' position in sequence A
    //! @params j1 5' position in sequence A
    //! @params i2 5' position in reverse sequence B
    //! @params j2 3' position in reverse sequence B
    //! @returns Boltzmann weight of interaction loop closed by i1.i2 
    //!          base pair with interior base pair j1.j2
    pf_t
    exp_ILoopE(size_t i1, size_t j1, size_t i2, size_t j2) const;
    

    //! pair type for pair of interacting positions
    //! @param i1 position in seqA
    //! @param i1 position in seqB
    //! @returns pair type of i1.i2
    int
    pair_type(size_t i1, size_t i2) const;

};

HybridPF::HybridPF(const std::string &seqA_,const std::string &seqB_):
    seqA(seqA_),
    seqB(seqB_),
    lenA(seqA.length()),
    lenB(seqB.length())
{
    
    create_temporary();
    
    compute_hybrid_pf();
    
    free_temporary();
}


void
HybridPF::create_temporary() {   

    SA = encode_sequence(seqA.c_str(),0);
    SA1 = encode_sequence(seqA.c_str(),1);
    
    std::string revseqB = seqB;
    reverse(revseqB.begin(),revseqB.end());
    
    SB = encode_sequence(revseqB.c_str(),0);
    SB1 = encode_sequence(revseqB.c_str(),1);

    pf_params = get_scaled_pf_parameters();

    make_pair_matrix();

}

void
HybridPF::free_temporary() {
    free(pf_params);
    free(SA);
    free(SA1);
    free(SB);
    free(SB1);
}

int
HybridPF::pair_type(size_t i1, size_t i2) const {
    assert(1 <= i1 && i1 <= lenA);
    assert(1 <= i2 && i2 <= lenB);
    
    return pair[SA[i1]][SB[i2]];
}

HybridPF::pf_t
HybridPF::exp_ILoopE(size_t i1, size_t i2, size_t k1,  size_t k2) const {
        
    int ptype_closing = pair_type(i1,i2);
    int ptype_enclosed = rtype[pair_type(k1,k2)]; // note: enclosed bp type
					   // 'turned around' for lib
					   // call
    
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
HybridPF::get_pf(size_t i1, size_t j1,size_t i2, size_t j2) const {
    assert(i1<=j1);
    assert(i2<=j2);
    return Q(i1,lenB-j2+1)(j1,lenB-i2+1);
}



int
main()
{
    // set some global variables for Vienna libRNA
    dangles=2;
    temperature = 37.;


    //                      0        1         2         3         4         5         6         7         8         9         0         1         2
    //                      12345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345
    const std::string seqA="GUUGGGAACUAGACCGAUCGCCAAUCCGUUUAUCUUUCAUAGAAGCCGGGAUUUAUCAGCUAUGUCGAAGAAUUUUAACUUGCUAUUGGGCACCCUGGUGGGGGUUAGUUUAGUUUUUCCCCAGG";

    //                      123456789012345678901234567
    const std::string seqB="ACAUAGCUGAUAAAUCCCGGCUUCUAU";


    UnpairedPF unpaired_pf_A(seqA);
    UnpairedPF unpaired_pf_B(seqB);
    HybridPF hybrid_pf(seqA,seqB);
    
    // print all subsequence hybrid partition functions (for test purposes)
    double kT = GASCONST/1000.0 * (K0+temperature);

    
    for (size_t i1=1; i1<=seqA.length(); i1++) {
	for (size_t j1=i1; j1<=seqA.length(); j1++) {
	    for (size_t i2=1; i2<=seqB.length(); i2++) {
		for (size_t j2=i2; j2<=seqB.length(); j2++) {
		    HybridPF::pf_t pfhyb = hybrid_pf.get_pf(i1,j1,i2,j2);
		    double ensemble_energy_hyb = - kT * log(pfhyb);
		    
		    UnpairedPF::pf_t upfA=unpaired_pf_A.get_pf(i1,j1);
		    UnpairedPF::pf_t upfB=unpaired_pf_B.get_pf(i2,j2);
		    double ensemble_energy_A = - kT * log(upfA);
		    double ensemble_energy_B = - kT * log(upfB);
		    
		    if (pfhyb >0) {
			std::cout << i1 << " " << j1 << " "
				  << i2 << " " << j2 <<" "
				  << pfhyb << " "
				  << upfA << " " 
			    	  << upfB << " " 
				  << ensemble_energy_hyb << " "
				  << ensemble_energy_A << " "
				  << ensemble_energy_B << " "
				  << std::endl;
		    }
		}
	    }
	}
    }
    
}
