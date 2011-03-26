#include <stdlib.h>

#include "hybrid_ensemble.hh"

#include "aux.hh"

// ----------------------------------------
// HybridEnsembleModel::StateDescription
//

HybridEnsembleModel::StateDescription::StateDescription() {
}


HybridEnsembleModel::StateDescription::StateDescription(size_t i1, size_t j1, size_t i2, size_t j2)
    : isites(1,ISite(i1,j1,i2,j2))
{
}


HybridEnsembleModel::StateDescription::StateDescription(size_t i1, size_t j1, size_t i2, size_t j2,
							size_t k1, size_t l1, size_t k2, size_t l2)
    : isites()
{
    isites.reserve(2);
    isites.push_back(ISite(i1,j1,i2,j2));
    isites.push_back(ISite(k1,l1,k2,l2));
}


size_t
HybridEnsembleModel::StateDescription::num_sites() const {
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
	std::cerr << "HybridEnsembleModel::StateDescription::encode: "
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
HybridEnsembleModel::StateDescription::encode(code_t &code) const {

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
HybridEnsembleModel::StateDescription::decode(const code_t &code) {
    assert(code.size()%4==0);
    
    isites.reserve(code.size()/4);
    
    for (size_t i=0; i<code.size()/4; i++) {
	isites.push_back(ISite(code[i*4+0],
			       code[i*4+1],
			       code[i*4+2],
			       code[i*4+3]));
    }    
}




// ----------------------------------------
// HybridEnsembleModel::StateDescription::NeighborIterator
//


HybridEnsembleModel::StateDescription::NeighborIterator::NeighborIterator(const StateDescription &origin_)
    :origin(origin_)
{    
}

const HybridEnsembleModel::StateDescription::NeighborIterator::ItState
HybridEnsembleModel::StateDescription::NeighborIterator::firstItState() const {
    ...
}


const HybridEnsembleModel::StateDescription::NeighborIterator::ItState
HybridEnsembleModel::StateDescription::NeighborIterator::nextItState(const ItState &itstate) const {
    
    if ( itstate.isites.size() == origin.num_sites() ) {
	// resize phase
    } else if ( itstate.isites.size() > origin.num_sites() ) {
	// insert new site phase
    } else if ( itstate.isites.size() < origin.num_sites() ) {
	// remove site phase
    }
}

void
HybridEnsembleModel::StateDescription::NeighborIterator::applyItState(const ItState &itstate, StateDescription &sd) const {
    ...
}

StateDescription
HybridEnsembleModel::StateDescription::NeighborIterator::applyItState(const ItState &itstate) const {
    ...
}
