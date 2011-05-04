#include "hybrid_ensemble_state.hh"


HybridEnsembleState* 
HybridEnsembleState::NeighborList::first(State::NeighborList::ItState** itstate) const {
    const HybridEnsembleModel::NeighborIterator::ItState &first=nit.firstItState();
    if (nit.isEndItState(first)) return NULL;
    
    *itstate = new ItState(first);
    
    HybridEnsembleState *elem=new HybridEnsembleState();
    
    nit.applyItState(first,elem->sd);
    
    return elem;
}

HybridEnsembleState*
HybridEnsembleState::NeighborList::next(State::NeighborList::ItState* itstate_,
			  State* elem_) const {
    
    HybridEnsembleState::NeighborList::ItState *itstate=
	static_cast<HybridEnsembleState::NeighborList::ItState *>(itstate_);
    
    HybridEnsembleState *elem = static_cast<HybridEnsembleState *>(elem_);
    
    nit.nextItState(*itstate);
    if (nit.isEndItState(*itstate)) return NULL;
    
    // if (elem!=NULL) delete elem;
    // HybridEnsembleState *
    // new_elem=new HybridEnsembleState(nit.applyItState(*itstate));
    // return new_elem;
    
    assert(elem);
    // Change the StateDescriptor in place in the state elem
    // This is potentially more dangerous than the above code,
    // since it makes assumptions about the internals of class
    // HybridEnsembleState
    // Note that this is /legal/ in ISO C++, since access to private members
    // is granted to nested classes.
    nit.applyItState(*itstate, elem->sd);
    
    return elem;
    
}

const std::string& 
HybridEnsembleState::getID() const {
    assert(false);
    std::string ID="";
    return ID;
}

// virtual desctructor
HybridEnsembleState::~HybridEnsembleState()
{}
