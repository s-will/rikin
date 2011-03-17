#include  <math.h>
// #include <assert.h>

#include "unpaired_pf.hh"

extern "C" {
#include "ViennaRNA/fold_vars.h"    
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/LPfold.h>
}


UnpairedPF::UnpairedPF(const std::string seq_) 
    : seq(seq_),
      length(seq.size())
{
    computePFs();
}

UnpairedPF::pf_t 
UnpairedPF::get_unpaired_prob_single(size_t i, size_t j) const {
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
