#include "hybrid_ensemble_model.hh"

/**
 * @file moves.cc
 *
 * @brief Moves and Move Iteration for HybridEnsembleModel
 */

// ------------------------------------------------------------
// Move
HybridEnsembleModel::Move::~Move() {}

// ------------------------------------------------------------
// Grow-Shrink Move

HybridEnsembleModel::GrowShrinkMove::GrowShrinkMove(const MoveIterator &mi): Move(mi) {}

HybridEnsembleModel::GrowShrinkMove::~GrowShrinkMove() {}

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

    const HybridEnsembleModel &m=mi.model();
        
    // constrain max and min by max loop size
    min1 = std::max(min1+m.maxunpinloop()+1,i1)-m.maxunpinloop()-1;
    min2 = std::max(min2+m.maxunpinloop()+1,i2)-m.maxunpinloop()-1;
    max1 = std::min(max1,i1+m.maxunpinloop()+1);
    max2 = std::min(max2,i2+m.maxunpinloop()+1);
    
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

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMove::transitionEnergy(size_t site, size_t left) const {
    
    const StateDescription &sd=mi.origin();
    StateDescription sd2(sd);
    apply(sd2);
    
    bool grow=( (left && k1<i1) || (!left && k1>i1) );

    const StateDescription &sd_small=grow ? sd  : sd2;
    const StateDescription &sd_large=grow ? sd2 : sd ;
    
    return
	left
	?
	transitionEnergy(sd_small,sd_large,k1,k2,i1,i2)
	:
	transitionEnergy(sd_small,sd_large,i1,i2,k1,k2);
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMove::transitionEnergy(const StateDescription &sd_small,
						      const StateDescription &sd_large,
						      size_t loop_i1, size_t loop_i2, size_t loop_j1, size_t loop_j2
						      ) const {
    // assert some of the preconditions
    assert(sd_small.num_sites()!=0);
    assert(sd_large.num_sites()==sd_small.num_sites());
        
    const HybridEnsembleModel &model = mi.model();
        
    energy_t e_hyb;
    energy_t e_loop;
    energy_t e_unp;
    
    	
    // pf of transition state all sites without added/removed loop
    e_hyb=
	model.energy_hybrid(sd_small.isites[0])
	+
	((sd_small.num_sites()==2)?model.energy_hybrid(sd_small.isites[1]):0);
	
    // Energy of added hybridization loop
    e_loop=model.energy_hybrid_loop(loop_i1,loop_i2,loop_j1,loop_j2);
    
    // energy difference for unpairing
    e_unp=
	(sd_small.num_sites()==1)
	?
	model.energy_unpair(sd_large.isites[0])
	:
	model.energy_unpair(sd_large.isites[0],
			    sd_large.isites[1]);
    
    return
	e_hyb
	+
	e_loop
	+
	e_unp;
}


// ------------------------------------------------------------
// Grow-Shrink Move First site Left end

HybridEnsembleModel::GrowShrinkMoveFL::GrowShrinkMoveFL(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybridEnsembleModel::GrowShrinkMoveFL::~GrowShrinkMoveFL() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveFL::nextMoveType() const {
    return new GrowShrinkMoveFR(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveFL::first() {
    
    // set to first move of this type if there is one
    
    const StateDescription &o=mi.origin();
    const HybridEnsembleModel &m=mi.model();

    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }

    return 
	HybridEnsembleModel::GrowShrinkMove::first(o.isites[0].i1,
						   o.isites[0].i2, 
						   1, 
						   1,
						   o.isites[0].j1-m.minsitesize()+1,
						   o.isites[0].j2-m.minsitesize()+1
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveFL::transitionEnergy() const {
    return
	HybridEnsembleModel::GrowShrinkMove::transitionEnergy(0,true);
}
    

void
HybridEnsembleModel::GrowShrinkMoveFL::apply(StateDescription &sd) const {
    sd.isites[0].i1=k1;
    sd.isites[0].i2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move First site Right end


HybridEnsembleModel::GrowShrinkMoveFR::GrowShrinkMoveFR(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybridEnsembleModel::GrowShrinkMoveFR::~GrowShrinkMoveFR() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveFR::nextMoveType() const {
    return new GrowShrinkMoveSL(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveFR::first() {
    
    // set to first move of this type if there is one
    const StateDescription &o=mi.origin();
    const HybridEnsembleModel &m=mi.model();
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }
    
    size_t len1 = m.seqA().length();
    size_t len2 = m.seqB().length();

    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   o.isites[0].j1,
						   o.isites[0].j2, 
						   o.isites[0].i1,
						   o.isites[0].i2,
						   (o.num_sites()==2)
						   ? o.isites[1].i1-m.minsitesize()+1
						   : len1,
						   (o.num_sites()==2)
						   ? o.isites[1].i2-m.minsitesize()+1
						   : len2
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveFR::transitionEnergy() const {
    return
	HybridEnsembleModel::GrowShrinkMove::transitionEnergy(0,false);
}

void
HybridEnsembleModel::GrowShrinkMoveFR::apply(StateDescription &sd) const {
    sd.isites[0].j1=k1;
    sd.isites[0].j2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Left end

HybridEnsembleModel::GrowShrinkMoveSL::GrowShrinkMoveSL(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybridEnsembleModel::GrowShrinkMoveSL::~GrowShrinkMoveSL() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveSL::nextMoveType() const {
    return new GrowShrinkMoveSR(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveSL::first() {
    
    // set to first move of this type if there is one
    
    const StateDescription &o=mi.origin();
    const HybridEnsembleModel &m=mi.model();

    // fail if there are less than two sites
    if (o.num_sites() < 2) {
	return false;
    }
    
    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   o.isites[1].i1,
						   o.isites[1].i2, 
						   o.isites[0].j1+m.minsitedist(),
						   o.isites[0].j2+m.minsitedist(),
						   o.isites[1].j1-m.minsitesize()+1,
						   o.isites[1].j2-m.minsitesize()+1
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveSL::transitionEnergy() const {
    return
	HybridEnsembleModel::GrowShrinkMove::transitionEnergy(1,true);
}

void
HybridEnsembleModel::GrowShrinkMoveSL::apply(StateDescription &sd) const {
    sd.isites[1].i1 = k1;
    sd.isites[1].i2 = k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Right end

HybridEnsembleModel::GrowShrinkMoveSR::GrowShrinkMoveSR(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybridEnsembleModel::GrowShrinkMoveSR::~GrowShrinkMoveSR() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMoveSR::nextMoveType() const {
    return new StopMove(mi);
}

bool
HybridEnsembleModel::GrowShrinkMoveSR::first() {
    
    // set to first move of this type if there is one

    const StateDescription &o=mi.origin();
    const HybridEnsembleModel &m=mi.model();

    size_t len1 = m.seqA().length();
    size_t len2 = m.seqB().length();
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }
    
    return 
	HybridEnsembleModel::GrowShrinkMove::first(
						   o.isites[1].j1,
						   o.isites[1].j2, 
						   o.isites[1].i1+m.minsitesize()+1,
						   o.isites[1].i2+m.minsitesize()+1,
						   len1,
						   len2
						   );
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMoveSR::transitionEnergy() const {    
    return
	HybridEnsembleModel::GrowShrinkMove::transitionEnergy(1,false);
}

void
HybridEnsembleModel::GrowShrinkMoveSR::apply(StateDescription &sd) const {
    sd.isites[1].j1=k1;
    sd.isites[1].j2=k2;
}


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
HybridEnsembleModel::StopMove::transitionEnergy() const {
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


HybridEnsembleModel::MoveIterator::
MoveIterator(const StateDescription &origin,
	     const HybridEnsembleModel &model)
    :origin_(origin),
     model_(model)
{
}


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
