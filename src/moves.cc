#include "hybrid_ensemble_model.hh"

/**
 * @file moves.cc
 *
 * @brief Moves and Move Iteration for HybridEnsembleModel
 */

// ------------------------------------------------------------
// Grow-Shrink Move

GrowShrinkMove()
{    
}

bool
HybridEnsembleModel::GrowShrinkMove::firstLeft() {        
    // try add/remove largest loop to the left
    k1 = min1;
    k2 = min2;
    
    return (k1<i1 && k2<i2);
}
    
bool
HybridEnsembleModel::GrowShrinkMove::firstRight() {        
    // try add/remove smallest loop to the right
    k1 = std::min(max1,i1+1);
    k2 = std::min(max2,i2+1);
    
    // succeed, if there is a move to the right, otherwise fail
    return (k1>i1 && k2>i2);
}

bool
HybridEnsembleModel::GrowShrinkMove::first(size_t i1_,
					   size_t i2_,
					   size_t min1_,
					   size_t min2_,
					   size_t max1_,
					   size_t max2_) {
    i1=i1_;
    i2=i2_;
    min1=min1_;
    min2=min2_;
    min1=min1_;
    max2=max2_;
    
    // constrain max and min by max loop size
    min1 = std::max(min1+maxunpinloop+1,i1)-maxunpinloop-1;
    min2 = std::max(min2+maxunpinloop+1,i2)-maxunpinloop-1;
    max1 = std::min(max1,i1+maxunpinloop+1);
    max2 = std::min(max2,i2+maxunpinloop+1);
    
    return firstLeft() || firstRight();
}

bool
HybridEnsembleModel::GrowShrinkMove::nextLeft() {
    assert(k1<i1);
       
    // mode "grow/shrink to left"
    
    k2++;
    if (k2<i2) {
	return true;
    } else {
	k1++;
	k2=min2;
	return (k1<i1);
    }
}

bool
HybridEnsembleModel::GrowShrinkMove::nextRight() {
    assert(k1>i1);

    // mode "grow/shrink to right"
    
    assert(max2 > i2); // this code is accessed only after success of first or
    // mode switch in next, in both cases, the condition is implied
    
    k2++;
    if (k2<max2) {
	return true;
    } else {
	k1++;
	k2 = i2+1;
	
	return (k1<=max1);
    }
}

bool
HybridEnsembleModel::GrowShrinkMove::next() {        
    if (k1<i1) {
	return nextLeft() || firstRight();
    } else {
	return nextRight();
    }
}


// ------------------------------------------------------------
// Grow-Shrink Move First site Left end

HybridEnsembleModel::GrowShrinkMoveFL::~GrowShrinkMoveFL() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveFL::nextMoveType() const {
    return new GrowShrinkMoveFR(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveFL::first() {
    
    // set to first move of this type if there is one
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }

    return 
	HybridEnsembleModel::GrowShrinkMove::first(mi.origin().isites[0].i1,
						   mi.origin().isites[0].i2, 
						   1, 
						   1,
						   o.isites[0].j1-minsize+1,
						   o.isites[0].j2-minsize+1
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveFL::transitionEnergy(const StateDescription &sd) const {
}

void
HybridEnsembleModel::GrowShrinkMoveFL::apply(StateDescription &sd) const {
    sd.isites[0].i1=k1;
    sd.isites[0].i2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move First site Right end

HybridEnsembleModel::GrowShrinkMoveFR::~GrowShrinkMoveFR() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveFR::nextMoveType() const {
    return new GrowShrinkMoveSL(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveFR::first() {
    
    // set to first move of this type if there is one
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }
    
    size_t len1 = mi.model().seqA().length();
    size_t len2 = mi.model().seqB().length();

    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   mi.origin().isites[0].j1,
						   mi.origin().isites[0].j2, 
						   mi.origin().isites[0].i1,
						   mi.origin().isites[0].i2,
						   (o.num_sites()==2)
						   ? o.isites[1].i1-minsize+1
						   : len1,
						   (o.num_sites()==2)
						   ? o.isites[1].i2-minsize+1
						   : len2
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveFR::transitionEnergy(const StateDescription &sd) const {
}

void
HybridEnsembleModel::GrowShrinkMoveFR::apply(StateDescription &sd) const {
    sd.isites[0].j1=k1;
    sd.isites[0].j2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Left end

HybridEnsembleModel::GrowShrinkMoveSL::~GrowShrinkMoveSL() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveSL::nextMoveType() const {
    return new GrowShrinkMoveSR(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveSL::first() {
    
    // set to first move of this type if there is one
    
    // fail if there are less than two sites
    if (o.num_sites() < 2) {
	return false;
    }
    
    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   mi.origin().isites[1].i1,
						   mi.origin().isites[1].i2, 
						   mi.origin().isites[0].j1+mindistance,
						   mi.origin().isites[0].j2+mindistance,
						   o.isites[1].j1-minsize+1,
						   o.isites[1].j2-minsize+1
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveSL::transitionEnergy(const StateDescription &sd) const {
}

void
HybridEnsembleModel::GrowShrinkMoveSL::apply(StateDescription &sd) const {
    sd.isites[1].i1 = k1;
    sd.isites[1].i2 = k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Right end

HybridEnsembleModel::GrowShrinkMoveSR::~GrowShrinkMoveSR() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveSR::nextMoveType() const {
    return new SplitMove(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveSR::first() {
    
    // set to first move of this type if there is one

    size_t len1 = mi.model().seqA().length();
    size_t len2 = mi.model().seqB().length();
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }
    
    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   mi.origin().isites[1].j1,
						   mi.origin().isites[1].j2, 
						   mi.origin().isites[1].i1+minsize+1,
						   mi.origin().isites[1].i2+minsize+1,
						   len1,
						   len2
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveSR::transitionEnergy(const StateDescription &sd) const {
}

void
HybridEnsembleModel::GrowShrinkMoveSR::apply(StateDescription &sd) const {
    sd.isites[1].j1=k1;
    sd.isites[1].j2=k2;
}


.
.
.

// ------------------------------------------------------------
// Stop Move

HybridEnsembleModel::StopMove::~StopMove() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::StopMove::nextMoveType() const {
    return NULL;
}

bool
HybridEnsembleModel::StopMove::first() {
    return false;
}

bool
HybridEnsembleModel::StopMove::next() {
    // this should never be called
    assert(false);
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::StopMove::transitionEnergy(const StateDescription &sd) const {
    // this should never be called
    assert(false);
}

void
HybridEnsembleModel::StopMove::apply(StateDescription &sd) const {
    // this should never be called
    assert(false);
}


// ------------------------------------------------------------
// Move Iterator

HybridEnsembleModel::Move *
HybridEnsembleModel::MoveIterator::firstMove() const {
    Move *m = new GrowShrinkMoveFL(*this);
    
    while( ! m->first() ) {
	Move *m2 = m->nextMoveType();
	delete m;
	m=m2;
    }
    return m;
}


HybridEnsembleModel::Move *
HybridEnsembleModel::MoveIterator::nextMove(Move *m) const {
    assert( m != NULL );

    if ( m->next() ) {
	return m;
    };
	    
    do {
	Move *m2 = m->nextMoveType(); // get new move object of next move type
	delete m; // and delete old object
	if (m2==NULL) {return NULL;} // return NULL, if there is no next move type
	m=m2;
    } while( ! m->first() ); // there is no first move of the next move type 
    
    return m;
}

void
HybridEnsembleModel::MoveIterator::disposeMove(Move *move) const {
    if (move!=NULL) delete move;
}
