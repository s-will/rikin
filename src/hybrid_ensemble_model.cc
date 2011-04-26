#include <stdlib.h>
#include <math.h>


#include "hybrid_ensemble_model.hh"

#include "aux.hh"

// ----------------------------------------
// HybEnsModel::StateDescription
//

HybEnsModel::StateDescription::StateDescription() 
    :isites(0)
{
}


HybEnsModel::StateDescription::StateDescription(size_t i1, size_t j1, size_t i2, size_t j2)
    : isites(1,ISite(i1,i2,j1,j2))
{
}


HybEnsModel::StateDescription::StateDescription(size_t i1, size_t j1, size_t i2, size_t j2,
						size_t k1, size_t l1, size_t k2, size_t l2)
    : isites()
{
    isites.reserve(2);
    isites.push_back(ISite(i1,i2,j1,j2));
    isites.push_back(ISite(k1,k2,l1,l2));
}


size_t
HybEnsModel::StateDescription::num_sites() const {
    return isites.size();
}


// --------------------
// encoding and decoding compressed representation
//
// use very simple code:
// always use 1 byte per position and 8 bytes in total (two sites),
// encode empty sites by 4 0-entries.
// fail if positions are too large.
//

// implement the test for too large positions
unsigned char
fail_on_toolarge(size_t x) {
    if (x>255) {
	std::cerr << "HybEnsModel::StateDescription::encode: "
		  << "Encoding does not support positions >=256."<<std::endl
		  << "  Use smaller sequence or implement more flexible coding." << std::endl;
	exit(-1);
    }
    return x;
}

// idea for more flexible code:
// the encoding is a sequence of the positions in the vector isites in
// the order specified by vector and ISite
// the first unsigned character encodes the number of sites and
// the number of bits per sequence position
//
//
void
HybEnsModel::StateDescription::encode(code_t &code) const {

    code.resize(num_sites()*2,0);

    for (size_t i=0; i<num_sites(); i++) {
	code[i*4+0]=fail_on_toolarge(isites[i].i1);
	code[i*4+1]=fail_on_toolarge(isites[i].i2);
	code[i*4+2]=fail_on_toolarge(isites[i].j1);
	code[i*4+3]=fail_on_toolarge(isites[i].j2);
    }
    
    // sketch for more flexible encoding/decoding
    // code.bitresize(0);
    // code.push_back(num_sites(),2); // number of sites
    // // determine maximal position in the vector
    // size_t maxpos=...;
    // size_t bitsperpos=(int)ceiling(log(maxpos)/log(2));
    // code.push_back(bitsperpos,6); // bits per element
    
}

void
HybEnsModel::StateDescription::decode(const code_t &code) {
    assert(code.size()%4==0);
    
    isites.reserve(code.size()/4);
    
    for (size_t i=0; i<code.size()/4; i++) {
	isites.push_back(ISite(code[i*4+0],
			       code[i*4+1],
			       code[i*4+2],
			       code[i*4+3]));
    }
}


// ------------------------------------------------------------
// Implementation of HybEnsModel

HybEnsModel::HybEnsModel(std::string seqA, std::string seqB)
    : seqA_(seqA),
      seqB_(seqB),
      uppfA_(seqA),
      uppfB_(seqB),
      hybridpf_(seqA,seqB),
      maxunpinloop_(5),
      minsitesize_(3),
      minsitedist_(6)
{
}


HybEnsModel::energy_t
HybEnsModel::energy_hybrid(size_t i1,size_t i2,size_t j1,size_t j2) const {
    return - hybridpf_.RT() * log( hybridpf_.partition_function(i1,i2,j1,j2) );
}

HybEnsModel::energy_t
HybEnsModel::energy_unpair(const StateDescription::ISite &is) const {
    return
	- uppfA_.RT() * log( uppfA_.unpaired_prob_single(is.i1,is.j1) )
	- uppfB_.RT() * log( uppfB_.unpaired_prob_single(is.i2,is.j2) );
}


HybEnsModel::energy_t
HybEnsModel::energy_unpair(const StateDescription::ISite &is1,const StateDescription::ISite &is2) const {
    return
	- uppfA_.RT() * log( uppfA_.unpaired_prob_joint(is1.i1,is1.j1,is2.i1,is2.j1) ) //seq 1
	- uppfB_.RT() * log( uppfB_.unpaired_prob_joint(is1.i2,is1.j2,is2.i2,is2.j2) ) //seq 2
	;
}

