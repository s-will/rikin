#include <stdlib.h>

#include "hybrid_ensemble_model.hh"

#include "aux.hh"

// ----------------------------------------
// HybridEnsembleModel::StateDescription
//

HybridEnsembleModel::StateDescription::StateDescription() 
    :isites(0)
{
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
// HybridEnsembleModel::ILMove/ILMoveIterator
//

HybridEnsembleModel::ILMove::mk_split=4;
HybridEnsembleModel::ILMove::mk_newsite=5;
HybridEnsembleModel::ILMove::mk_empty=6;



HybridEnsembleModel::ILMoveIterator::
ILMoveIterator(const StateDescription &origin_, size_t maxunpinloop_)
    :origin(origin_),
     maxunpinloop(maxunpinloop_)
{
}

HybridEnsembleModel::ILMove
HybridEnsembleModel::ILMoveIterator::firstMove() const {
    HybridEnsembleModel::ILMove m;
    
    m=firstGrowShrinkMove();
    
    if (isEndMove(m)) {
	m=firstSplitMove();
    }
    
    if (isEndMove(m)) {
	m=firstNewSiteMove();
    }
    
    return m;
}


const HybridEnsembleModel::ILMove &
HybridEnsembleModel::ILMoveIterator::nextMove(ILMove &m) const {

    assert(!isEndMove(m));
    
    size_t mk=m.move_kind;
    
    if (mk < mk_split) {
	nextGrowShrinkMove(m);
	
	if (m.isEndMove()) {
	    mk=mk_split;
	} else {
	    return m;
	}
    }
    
    if (mk == mk_split) {
	if (m.isEndMove()) {
	    m=firstSplitMove();
	} else {
	    nextSplitMove(m);
	}
	
	if (m.isEndMove()) {
	    mk=mk_newsite;
	} else {
	    return m;
	}
    }

    if (mk==mk_newsite) {
	if (m.isEndMove()) {
	    m=firstNewSiteMove();
	} else {
	    nextNewSiteMove(m);
	}
    }
    
    return m;
}

HybridEnsembleModel::ILMove
firstGrowShrinkMove() const {
    ILMove m;
    
    //unless empty (=non-interacting) state
    if (origin.num_sites()>0) {
	m.set_empty();
    }

    // set to growing left end of first site
    move_kind=0;
    
    xs.resize(4);
    
    // specify maximal (due to maxunpinloop) allowed growth
    xs[0]=max(origin.sites[0].i1,maxunpinloop)-maxunpinloop;
    xs[1]=max(origin.sites[0].i2,maxunpinloop2)-maxunpinloop2;
    
    xs[2]=origin.sites[0].j1;
    xs[3]=origin.sites[0].j2;
    
    
    // we assume that there is always a grow shrink move at the first site.
    // This fails in the case of very small sequences and sites that span the whole sequences.

    return m;
}


void
HybridEnsembleModel::ILMove::applyMove(StateDescription &sd) const {
    assert(false);
}

bool
HybridEnsembleModel::ILMoveIterator::isEndMove(const ILMove &move) const {
    return move.is_empty();
}

double
HybridEnsembleModel::Move::
energyOfTransitionState(StateDescription &sd,
			const HybridEnsembleModel &model) const {
    assert(false);
}
