#include  <math.h>
// #include <assert.h>

#include "unpaired_pf.hh"

extern "C" {
#include "ViennaRNA/fold_vars.h"    
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/utils.h>
#include <LPfold.h>
}


UnpairedPF::UnpairedPF(const std::string &seq,
		       int orientation,
		       size_t maxsitesize,
		       size_t span,
		       size_t window,
		       bool cond) 
    : seq_(seq),
      //orientation_(orientation),
      RT_( (temperature+K0)*GASCONST/1000.0 ),
      maxsitesize_(maxsitesize),
      span_(span),
      window_(window),
      cond_(cond)
{
    assert(orientation==-1 || orientation==1);
    
    if (orientation==-1) { std::reverse(seq_.begin(),seq_.end()); }
    
    // std::cout << "Create UnpairedPF from sequence "<<seq<<" ("<<seq_.length()<<")"<<std::endl;
    computeSingleProbs();
    if (cond_) { computeCondProbs(); }
    
    if  (orientation==-1) {
	revert_single_probs();
	if (cond_) { revert_cond_probs(); }
    }
}

template<class M>
void
UnpairedPF::revert_upper_triangle_matrix(size_t size, M &m) {
    for (size_t i=1; i<=size/2; i++) {
	for (size_t j=i; j<=(size-i); j++) {
	    size_t i1=size-j+1;
	    size_t j1=size-i+1;
	    std::swap(m(i,j),m(i1,j1));
	}
    }
}

void
UnpairedPF::revert_single_probs() {
    revert_upper_triangle_matrix(seq_.length(),Psingle);
}

void
UnpairedPF::revert_cond_probs() {
    for (size_t i=1; i<=seq_.length(); i++) {
	for (size_t j=i+1; (j-i+1<=maxsitesize_) && (j<=seq_.length()); j++) {
	    revert_upper_triangle_matrix(seq_.length(),Pcond(i,j));
	}
    }
    revert_upper_triangle_matrix(seq_.length(),Pcond);
}

UnpairedPF::prob_t 
UnpairedPF::unpaired_prob_single(size_t i, size_t j) const {
    assert(1<=i);
    assert(i<=j);
    assert(j<=seq_.length());
    
    assert(j-i+1 <= maxsitesize_);
    
    return Psingle(i,j);
}

UnpairedPF::prob_t
UnpairedPF::unpaired_prob_conditional(size_t i, size_t j,size_t k, size_t l) const {
    assert(1<=i);
    assert(i<=j);
    assert(j<=seq_.length());
    
    assert(1<=k);
    assert(k<=l);
    assert(l<=seq_.length());
 
    assert(j<k);
    
    assert(j-i+1 <= maxsitesize_);
    assert(l-k+1 <= maxsitesize_);
    
    if (cond_) {
	return Pcond(k,l)(i,j);
    } else {
	return Psingle(i,j);
    }
}



UnpairedPF::pf_t
UnpairedPF::total_pf() const {
    const char *sequence=seq_.c_str();
    //char structure[seq_.length()+2];
    
    /* for longer sequences one should also set a scaling factor for
       partition function folding */
    // double e=fold(sequence,structure);
    // free_arrays();
    
    // double kT = (temperature+K0)*GASCONST/1000.;  /* kT in kcal/mol */
    // pf_scale = exp(-e/kT/seq_.length());

    pf_scale = -1;
    
    //update_pf_params(seq_.length());
    
    // calculate Gibb's free energy in kcal/mol
    double G = pf_fold(sequence,NULL);

    // free_pf_arrays();
    
    return exp( - G / (GASCONST/1000.0*(temperature+K0)) );
}


void
UnpairedPF::computeSingleProbs() {
    computeProbsGeneric(Psingle,NULL);    
}

void
UnpairedPF::computeProbsGeneric(LocARNA::Matrix<prob_t> &P, const char *structure) {
    
    int fold_constrained_before=fold_constrained;
    if (structure!=NULL) {fold_constrained=-1;}
    
    int plfW=std::min(window_,seq_.length()); // parameter -W of plfold
    int plfL=std::min(span_,seq_.length()); // parameter -L of plfold
    
    // maximal size of unpaired region for which pf is computed
    int maxUnpairedRegionSize=std::min(seq_.length(),maxsitesize_);

    const char *sequence=seq_.c_str();
    
    float cutoff = 1.0;    /* bpcutoff for plfold, don't keep base pairs in pl */
    
    // declare and alloc pup array
    double ** pup = (double **)space( sizeof(double *) * (seq_.length()+1) );
    
    pup[0] = (double *)space( sizeof(double) ); /* we need only entry 0 */
    pup[0][0] = (double) maxUnpairedRegionSize;
    
    // don't compute conditional probabilities
    plist *dpp=NULL;
    
    pf_scale = 1;
    
    //! call libRNA
    plist *pl = pfl_foldC(const_cast<char *>(sequence), structure, plfW, plfL, cutoff, pup, &dpp, NULL, NULL);
    
    // copy result to matrix Psingle
    P.resize(seq_.length()+1,seq_.length()+1);
    P.fill(0.0);
    // i+u: end position of unpaired region counting from 0, u: length of unpaired region
    for(int u=1; u<=maxUnpairedRegionSize; u++) {
	for(size_t i=1; i<=seq_.length()-u+1; i++) {
	    P(i,i+u-1)=pup[i+u-1][u];
	}
    }
    
    // free pl
    if (pl) free(pl);
    
    // free pup array
    for (size_t i=0;i<=seq_.length();i++) {free(pup[i]);}
    //free(pup[0]);
    free(pup);

    // restore state of global variable
    fold_constrained=fold_constrained_before;
}


void
UnpairedPF::computeCondProbs(size_t i, size_t j) {
    assert(i<=j);
    assert(1<=i);
    assert(j<=seq_.length());
    assert(j-i+1 <= maxsitesize_);
    
    // generate constraint structure string
    // of the form .....xxxxx......, where the string
    // contains x at positions i-1..j-1
    std::string structure="";
    for (size_t x=1; x<i; x++) structure+='.';
    for (size_t x=i; x<=j; x++) structure+='x';
    for (size_t x=j+1; x<=seq_.length(); x++) structure+='.';
        
    computeProbsGeneric(Pcond(i,j), structure.c_str());
    
}


void
UnpairedPF::computeCondProbs() {
    
    Pcond.resize(seq_.length()+1,seq_.length()+1);
    
    //std::cout << "Compute conditional unpaired probabilities ..." <<std::endl;
    for (size_t i=1; i<=seq_.length(); i++) {
	//std::cout <<"  "<<i<<" "<<i<<".."<<seq_.length()<<std::endl;
	for (size_t j=i+1; (j-i+1 <= maxsitesize_) && j<=seq_.length(); j++) {
	    
	    computeCondProbs(i,j);
	    
	    /*
	    // report the strongest dependencies for debugging
	    
	    for (size_t k=1; k<=seq_.length(); k++) {
		for (size_t l=k+1; l<=seq_.length(); l++) {
		    
		    bool overlap = (i<=k && k<=j) || (i<=l && l<=j); 
		    
	   	    if (!overlap && Pcond(i,j)(k,l)>0.7) {
			std::cout << i << " "<< j << " "<< k << " "<< l << " " << Pcond(i,j)(k,l) << std::endl;
		    }
		}
	    }
	    */
	}
    }
}
