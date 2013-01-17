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


HybEnsModel::StateDescription::StateDescription(size_t i1, size_t i2, size_t j1, size_t j2)
    : isites(1,ISite(i1,i2,j1,j2))
{
}


HybEnsModel::StateDescription::StateDescription(size_t i1, size_t i2, size_t j1, size_t j2,
						size_t k1, size_t k2, size_t l1, size_t l2)
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

std::string
HybEnsModel::StateDescription::toString() const {
    std::ostringstream out;
    out << (*this);
    return out.str();
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
HybEnsModel::StateDescription::code_t &
HybEnsModel::StateDescription::encode(code_t &code) const {

    code.resize(num_sites()*4);
    
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
    
    return code;
}

HybEnsModel::StateDescription::code_t
HybEnsModel::StateDescription::encode() const {

    std::string code;
    
    encode(code);
    
    return code;
}



HybEnsModel::StateDescription &
HybEnsModel::StateDescription::decode(const code_t &code) {
    assert(code.size()%4==0);
    
    isites.reserve(code.size()/4);
    isites.resize(0);

    for (size_t i=0; i<code.size()/4; i++) {
	isites.push_back(ISite((unsigned char)code[i*4+0],
			       (unsigned char)code[i*4+1],
			       (unsigned char)code[i*4+2],
			       (unsigned char)code[i*4+3]));
    }
    return *this;
}

bool
HybEnsModel::StateDescription::is_valid(const HybEnsModel &model) const {
    bool valid=true;
    
    // at most 2 sites
    valid = valid && num_sites()<=2;
	    
    // check site size and boundaries
    for (size_t i=0; valid && i<num_sites(); ++i) {
	valid = valid
	    && 1 <= isites[i].i1 
	    && isites[i].i1 <= isites[i].j1
	    && isites[i].j1 <= model.seqA().length();
	valid = valid 
	    && 1 <= isites[i].i2
	    && isites[i].i2 <= isites[i].j2 
	    && isites[i].j2 <= model.seqB().length();
	valid = valid 
	    && isites[i].j1  >= isites[i].i1 + model.minsitesize() - 1;
	valid = valid
	    && isites[i].j2  >= isites[i].i2 + model.minsitesize() - 1;
    }
	    
    // check site distance
    if (num_sites()==2) {
	valid = valid
	    && isites[0].j1 + model.minsitedist() + 1 <= isites[1].i1;
	valid = valid
	    && isites[0].j2 + model.minsitedist() + 1 <= isites[1].i2;
    }

    return valid;
}


// ------------------------------------------------------------
// Implementation of HybEnsModel

HybEnsModel::HybEnsModel(std::string seqA, std::string seqB)
    : seqA_(seqA),
      seqB_(seqB),
      uppfA_(seqA,+1),
      uppfB_(seqB,-1),
      hybridpf_(seqA,seqB),
      maxunpinloop_(6),
      minsitesize_(3),
      minsitedist_(6)
{
}

HybEnsModel::energy_t
HybEnsModel::energy_hybrid(const StateDescription::ISite &is) const {
    return energy_hybrid(is.i1,is.i2,is.j1,is.j2);
}


HybEnsModel::energy_t
HybEnsModel::energy_hybrid(size_t i1,size_t i2,size_t j1,size_t j2) const {    
    return - hybridpf_.RT() * log( hybridpf_.partition_function(i1,j1,i2,j2) );
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


HybEnsModel::energy_t
HybEnsModel::energy(const StateDescription &sd) const {
    switch(sd.num_sites()) {
    case 0:
	return - energy_duplex_init()/100.0; // energy penalty for first interaction is added to empty state!!!
    case 1:
	// std::cout << "E = "<<energy_unpair(sd[0])<<"(unp) + "<<energy_hybrid(sd[0])<<"(hyb)"<<std::endl;
	return energy_unpair(sd[0])+energy_hybrid(sd[0]);
    case 2:
	return energy_unpair(sd[0],sd[1]) + energy_hybrid(sd[0]) + energy_hybrid(sd[1]);
    default:
	assert(false);
    }
}


std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription::ISite &isite) {
    out << "("
	<< isite.i1 << "," 
	<< isite.i2 << ")-(" 
	<< isite.j1 << "," 
	<< isite.j2 << ")";
    return out;
}

std::ostream &
operator << (std::ostream &out, const HybEnsModel::StateDescription &sd) {
    out << "{";
    for (size_t i=0; i<sd.num_sites(); i++) {
	out << sd[i];
	if (i<sd.num_sites()-1) {out << ",";}
    }
    out << "}";
    return out;
}

