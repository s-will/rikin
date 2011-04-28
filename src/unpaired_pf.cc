#include  <math.h>
// #include <assert.h>

#include "unpaired_pf.hh"

extern "C" {
#include "ViennaRNA/fold_vars.h"    
#include <ViennaRNA/part_func.h>
#include <LPfold.h>
}


UnpairedPF::UnpairedPF(const std::string &seq_,int orientation_) 
    : seq(seq_),
      orientation(orientation_),
      RT_( (temperature+K0)*GASCONST/1000.0 )
{
    assert(orientation==-1 || orientation==1);
    
    if (orientation==-1) { std::reverse(seq.begin(),seq.end()); }
    
    // std::cout << "Create UnpairedPF from sequence "<<seq<<" ("<<seq.size()<<")"<<std::endl;
    computeSingleProbs();
    computeCondProbs();
}

UnpairedPF::prob_t 
UnpairedPF::unpaired_prob_single(size_t i, size_t j) const {
    assert(1<=i);
    assert(i<=j);
    assert(j<=seq.size());
    
    if (orientation==-1) {
	i=seq.length()-i+1;
	j=seq.length()-j+1;
	std::swap(i,j);
    }
    
    return Psingle(i,j);
}

UnpairedPF::prob_t
UnpairedPF::unpaired_prob_conditional(size_t i, size_t j,size_t k, size_t l) const {
    if (orientation==-1) {
	i=seq.length()-i+1;
	j=seq.length()-j+1;
	std::swap(i,j);
	k=seq.length()-k+1;
	l=seq.length()-l+1;
	std::swap(k,l);
    }
    
    return Pcond(k,l)(i,j);
}



UnpairedPF::pf_t
UnpairedPF::total_pf() const {
    const char *sequence=seq.c_str();
    //char structure[seq.size()+2];
    
    /* for longer sequences one should also set a scaling factor for
       partition function folding */
    // double e=fold(sequence,structure);
    // free_arrays();
    
    // double kT = (temperature+K0)*GASCONST/1000.;  /* kT in kcal/mol */
    // pf_scale = exp(-e/kT/seq.size());

    pf_scale = -1;
    
    //update_pf_params(seq.size());
    
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
    
    int plfW=seq.size(); // parameter -W of plfold
    int plfL=seq.size(); // parameter -L of plfold
    
    // maximal size of unpaired region for which pf is computed
    int maxUnpairedRegionSize=seq.size();
    
    const char *sequence=seq.c_str();
    
    float cutoff = 0.0;    /* bpcutoff for plfold, does not influence
			      our results */
    
    // declare and alloc pup array
    double ** pup = new double *[seq.size()+1];
    
    pup[0] = new double[1]; /* we need only entry 0 */
    pup[0][0] = (double) maxUnpairedRegionSize;
    
    // don't compute conditional probabilities
    plist *dpp=NULL;
    
    pf_scale = -1;
    
    //! call libRNA
    pfl_fold(const_cast<char *>(sequence), structure, plfW, plfL, cutoff, pup, &dpp, NULL, NULL);
    
    // copy result to matrix Psingle
    P.resize(seq.size()+1,seq.size()+1);
    P.fill(0.0);
    // i+u: end position of unpaired region counting from 0, u: length of unpaired region
    for(int u=1; u<=maxUnpairedRegionSize; u++) {
	for(size_t i=1; i<=seq.size()-u+1; i++) {
	    P(i,i+u-1)=pup[i+u-1][u];
	}
    }
    
    // free pup array
    for (size_t i=0;i<=seq.size();i++) {delete pup[i];}
    delete pup;
    
    // restore state of global variable
    fold_constrained=fold_constrained_before;
}


void
UnpairedPF::computeCondProbs(size_t i, size_t j) {
    assert(i<=j);
    assert(1<=i);
    assert(j<=seq.size());
    
    // generate constraint structure string
    // of the form .....xxxxx......, where the string
    // contains x at positions i-1..j-1
    std::string structure="";
    for (size_t x=1; x<i; x++) structure+='.';
    for (size_t x=i; x<=j; x++) structure+='x';
    for (size_t x=j+1; x<=seq.size(); x++) structure+='.';
        
    computeProbsGeneric(Pcond(i,j), structure.c_str());
    
}


void
UnpairedPF::computeCondProbs() {
    
    Pcond.resize(seq.size()+1,seq.size()+1);
    
    //std::cout << "Compute conditional unpaired probabilities ..." <<std::endl;
    for (size_t i=1; i<=seq.size(); i++) {
	//std::cout <<"  "<<i<<" "<<i<<".."<<seq.size()<<std::endl;
	for (size_t j=i+1; j<=seq.size(); j++) {
	    
	    computeCondProbs(i,j);
	    
	    /*
	    // report the strongest dependencies for debugging
	    
	    for (size_t k=1; k<=seq.size(); k++) {
		for (size_t l=k+1; l<=seq.size(); l++) {
		    
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
