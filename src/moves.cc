#include "hybrid_ensemble_model.hh"

/**
 * @file moves.cc
 *
 * @brief Moves and Move Iteration for HybridEnsembleModel
 */

// ------------------------------------------------------------
// Grow-Shrink Move 0

HybridEnsembleModel::GrowShrinkMove0::~GrowShrinkMove0() {}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMove0::nextMoveType() const {
    return new GrowShrinkMove1(mi);
}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMove0::first() {
    // set to first move of this type if there is one
    
    // on succes 
    return this; 
}

HybridEnsembleModel::Move *
HybridEnsembleModel::GrowShrinkMove0::next() {
    // generate next move of this type if there is one
    
    // on succes 
    return this; 
}

HybridEnsembleModel::energy_t
HybridEnsembleModel::GrowShrinkMove0::transitionEnergy(const StateDescription &sd) const {
}

void
HybridEnsembleModel::GrowShrinkMove0::apply(const StateDescription &sd) const {
}

// ------------------------------------------------------------
// Grow-Shrink Move 1

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

HybridEnsembleModel::Move *
HybridEnsembleModel::StopMove::first() {
    delete this;
    return NULL;
}

HybridEnsembleModel::Move *
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
HybridEnsembleModel::StopMove::apply(const StateDescription &sd) const {
    // this should never be called
    assert(false);
}


// ------------------------------------------------------------
// Move Iterator

HybridEnsembleModel::Move *
HybridEnsembleModel::MoveIterator::firstMove() const {
    Move *m = new GrowShrinkMove0(*this);
    Move *m2;
    
    while( ( m2 = m->first() ) == NULL ) {
	m2 = m->nextMoveType();
	delete m;
	m=m2;
    }
    return m2;
}


HybridEnsembleModel::Move *
HybridEnsembleModel::MoveIterator::nextMove(Move *m) const {
    assert( m != NULL );
	    
    Move *m2 = m->next();
	    
    if (m2!=NULL) {
	return m2;
    };
	    
    do {
	m2 = m->nextMoveType(); // get new move object of next move type
	delete m; // and delete old object
	if (m2==NULL) {return NULL;} // return NULL, if there is no next move type
	m=m2;
    } while( ( m2 = m->first() ) == NULL ); // there is no first move of the next move type 
    
    return m2;
}

void
HybridEnsembleModel::MoveIterator::disposeMove(Move *move) const {
    if (move!=NULL) delete move;
}
