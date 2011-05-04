#include "hybrid_ensemble_model.hh"

/**
 * @file moves.cc
 *
 * @brief Moves and Move Iteration for HybEnsModel
 */

// ------------------------------------------------------------
// Move
HybEnsModel::Move::~Move() {}


std::ostream &
HybEnsModel::Move::print(std::ostream &out) const {
    out << "[Generic Move]";
    return out;
}

// ------------------------------------------------------------
// Grow-Shrink Move

HybEnsModel::GrowShrinkMove::GrowShrinkMove(const MoveIterator &mi): Move(mi) {}

HybEnsModel::GrowShrinkMove::~GrowShrinkMove() {}

std::ostream &
HybEnsModel::GrowShrinkMove::print(std::ostream &out) const {
    out << "[Grow Shrink Move "
	<<k1<<" "<<k2 //<<" "
	// <<min1<<" "<<min2<<" "
	// <<i1<<" "<<i2 <<" "
	// <<max1<<" "<<max2
	<<"]";
    return out;
}


bool
HybEnsModel::GrowShrinkMove::firstLeft() {        
    // try add/remove largest loop to the left
    k1 = min1;
    k2 = min2;
    
    return (k1<i1 && k2<i2);
}
    
bool
HybEnsModel::GrowShrinkMove::firstRight() {        
    // try add/remove smallest loop to the right
    k1 = std::min(max1,i1+1);
    k2 = std::min(max2,i2+1);
    
    // succeed, if there is a move to the right, otherwise fail
    return (k1>i1 && k2>i2);
}

bool
HybEnsModel::GrowShrinkMove::first(size_t i1_,
				   size_t i2_,
				   size_t min1_,
				   size_t min2_,
				   size_t max1_,
				   size_t max2_) {
    i1=i1_;
    i2=i2_;
    min1=min1_;
    min2=min2_;
    max1=max1_;
    max2=max2_;

    const HybEnsModel &m=mi.model();
        
    // constrain max and min by max loop size
    min1 = std::max(min1+m.maxunpinloop()+1,i1)-m.maxunpinloop()-1;
    min2 = std::max(min2+m.maxunpinloop()+1,i2)-m.maxunpinloop()-1;
    max1 = std::min(max1,i1+m.maxunpinloop()+1);
    max2 = std::min(max2,i2+m.maxunpinloop()+1);
    
    return firstLeft() || firstRight();
}

bool
HybEnsModel::GrowShrinkMove::nextLeft() {
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
HybEnsModel::GrowShrinkMove::nextRight() {
    assert(k1>i1);

    // mode "grow/shrink to right"
    
    assert(max2 > i2); // this code is accessed only after success of first or
    // mode switch in next, in both cases, the condition is implied
    
    k2++;
    if (k2<=max2) {
	return true;
    } else {
	k1++;
	k2 = i2+1;
	
	return (k1<=max1);
    }
}

bool
HybEnsModel::GrowShrinkMove::next() {        
    if (k1<i1) {
	return nextLeft() || firstRight();
    } else {
	return nextRight();
    }
}

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMove::transitionEnergy(size_t site, size_t left) const {
    
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

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMove::transitionEnergy(const StateDescription &sd_small,
					      const StateDescription &sd_large,
					      size_t loop_i1, size_t loop_i2, size_t loop_j1, size_t loop_j2
					      ) const {
    // assert some of the preconditions
    assert(sd_small.num_sites()!=0);
    assert(sd_large.num_sites()==sd_small.num_sites());
    
    /*
    std::cout <<std::endl << "GrowShrinkMove::transitionEnergy" << sd_small << " " << sd_large<< " " 
     	      << loop_i1 << " "<< loop_i2 << " " << loop_j1 << " " << loop_j2
     	      <<std::endl;
    */
    
    const HybEnsModel &model = mi.model();
    
    
    // pf of transition state all sites without added/removed loop
    energy_t e_hyb =
	model.energy_hybrid(sd_small[0])
	+
	((sd_small.num_sites()==2)?model.energy_hybrid(sd_small[1]):0);
    
    // Energy of added hybridization loop
    energy_t e_loop = 
	model.energy_hybrid_loop(loop_i1,loop_i2,loop_j1,loop_j2);
    
    // energy difference for unpairing
    energy_t e_unp =
	(sd_small.num_sites()==1)
	?
	model.energy_unpair(sd_large[0])
	:
	model.energy_unpair(sd_large[0],
			    sd_large[1]);

    //std::cout << "GrowShrinkMove::transitionEnergy returns " << e_hyb <<" + " << e_loop <<" + "<< e_unp << std::endl;

    return
	e_hyb
	+
	e_loop
	+
	e_unp;
}


// ------------------------------------------------------------
// Grow-Shrink Move First site Left end

HybEnsModel::GrowShrinkMoveFL::GrowShrinkMoveFL(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybEnsModel::GrowShrinkMoveFL::~GrowShrinkMoveFL() {}

HybEnsModel::Move *
HybEnsModel::GrowShrinkMoveFL::nextMoveType() const {
    return new GrowShrinkMoveFR(mi);
}

bool
HybEnsModel::GrowShrinkMoveFL::first() {
    
    // set to first move of this type if there is one
    
    const StateDescription &o=mi.origin();
    const HybEnsModel &m=mi.model();

    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }

    return 
	HybEnsModel::GrowShrinkMove::first(o[0].i1,
					   o[0].i2, 
					   1, 
					   1,
					   o[0].j1-m.minsitesize()+1,
					   o[0].j2-m.minsitesize()+1
					   );
}

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMoveFL::transitionEnergy() const {
    return
	HybEnsModel::GrowShrinkMove::transitionEnergy(0,true);
}
    

void
HybEnsModel::GrowShrinkMoveFL::apply(StateDescription &sd) const {
    sd[0].i1=k1;
    sd[0].i2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move First site Right end


HybEnsModel::GrowShrinkMoveFR::GrowShrinkMoveFR(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybEnsModel::GrowShrinkMoveFR::~GrowShrinkMoveFR() {}

HybEnsModel::Move *
HybEnsModel::GrowShrinkMoveFR::nextMoveType() const {
    return new GrowShrinkMoveSL(mi);
}

bool
HybEnsModel::GrowShrinkMoveFR::first() {
    
    // set to first move of this type if there is one
    const StateDescription &o=mi.origin();
    const HybEnsModel &m=mi.model();
    
    // fail if there is no site
    if (o.num_sites() == 0) {
	return false;
    }
    
    size_t len1 = m.seqA().length();
    size_t len2 = m.seqB().length();

    return 
	HybEnsModel::GrowShrinkMove::first(
					   o[0].j1,
					   o[0].j2, 
					   o[0].i1,
					   o[0].i2,
					   (o.num_sites()==2)
					   ? o[1].i1-m.minsitesize()+1
					   : len1,
					   (o.num_sites()==2)
					   ? o[1].i2-m.minsitesize()+1
					   : len2
					   );
}

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMoveFR::transitionEnergy() const {
    return
	HybEnsModel::GrowShrinkMove::transitionEnergy(0,false);
}

void
HybEnsModel::GrowShrinkMoveFR::apply(StateDescription &sd) const {
    sd[0].j1=k1;
    sd[0].j2=k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Left end

HybEnsModel::GrowShrinkMoveSL::GrowShrinkMoveSL(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybEnsModel::GrowShrinkMoveSL::~GrowShrinkMoveSL() {}

HybEnsModel::Move *
HybEnsModel::GrowShrinkMoveSL::nextMoveType() const {
    return new GrowShrinkMoveSR(mi);
}

bool
HybEnsModel::GrowShrinkMoveSL::first() {
    
    // set to first move of this type if there is one
    
    const StateDescription &o=mi.origin();
    const HybEnsModel &m=mi.model();

    // fail if there are less than two sites
    if (o.num_sites() < 2) {
	return false;
    }
    
    return 
	HybEnsModel::GrowShrinkMove::first(
					   o[1].i1,
					   o[1].i2, 
					   o[0].j1+m.minsitedist(),
					   o[0].j2+m.minsitedist(),
					   o[1].j1-m.minsitesize()+1,
					   o[1].j2-m.minsitesize()+1
					   );
}

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMoveSL::transitionEnergy() const {
    return
	HybEnsModel::GrowShrinkMove::transitionEnergy(1,true);
}

void
HybEnsModel::GrowShrinkMoveSL::apply(StateDescription &sd) const {
    sd[1].i1 = k1;
    sd[1].i2 = k2;
}

// ------------------------------------------------------------
// Grow-Shrink Move Second site Right end

HybEnsModel::GrowShrinkMoveSR::GrowShrinkMoveSR(const MoveIterator &mi):GrowShrinkMove(mi) {}

HybEnsModel::GrowShrinkMoveSR::~GrowShrinkMoveSR() {}

HybEnsModel::Move *
HybEnsModel::GrowShrinkMoveSR::nextMoveType() const {
    return new RemoveSiteMove(mi);
}


bool
HybEnsModel::GrowShrinkMoveSR::first() {
    
    // set to first move of this type if there is one

    const StateDescription &o=mi.origin();
    const HybEnsModel &m=mi.model();

    size_t len1 = m.seqA().length();
    size_t len2 = m.seqB().length();
    
    // fail if there are less then two sites
    if (o.num_sites() < 2)  {
	return false;
    }
    
    return 
	HybEnsModel::GrowShrinkMove::first(
					   o[1].j1,
					   o[1].j2, 
					   o[1].i1+m.minsitesize()+1,
					   o[1].i2+m.minsitesize()+1,
					   len1,
					   len2
					   );
}

HybEnsModel::energy_t
HybEnsModel::GrowShrinkMoveSR::transitionEnergy() const {    
    return
	HybEnsModel::GrowShrinkMove::transitionEnergy(1,false);
}

void
HybEnsModel::GrowShrinkMoveSR::apply(StateDescription &sd) const {
    sd[1].j1=k1;
    sd[1].j2=k2;
}

// ------------------------------------------------------------
// Move that removes a site
//

HybEnsModel::RemoveSiteMove::RemoveSiteMove(const MoveIterator &mi) : Move(mi) {}

HybEnsModel::RemoveSiteMove::~RemoveSiteMove() {}

std::ostream &
HybEnsModel::RemoveSiteMove::print(std::ostream &out) const {
    out << "[Remove Site Move "<<site<<"]";
    return out;
}


HybEnsModel::Move *
HybEnsModel::RemoveSiteMove::nextMoveType() const {
    return new NewSiteMoveF(mi);
}


bool
HybEnsModel::RemoveSiteMove::thisornext() {
    const HybEnsModel &m=mi.model();
       
    for(;site<mi.origin().num_sites();site++) {
	// test whether it is allowed to remove the site,
	// and return true on success
	
	const StateDescription::ISite &is=mi.origin()[site];
	
	if ((is.j1-is.i1-1 <= m.maxunpinloop())
	    &&
	    (is.j2-is.i2-1 <= m.maxunpinloop()))
	    {
		return true;
	    }
    }
    return false;
}

bool
HybEnsModel::RemoveSiteMove::first() {
    site=0;
    
    return thisornext();
}

bool
HybEnsModel::RemoveSiteMove::next() {
    site++;
    return thisornext();
}

HybEnsModel::energy_t
HybEnsModel::RemoveSiteMove::transitionEnergy() const {
    // make removed site interaction loop and score
    const HybEnsModel &model = mi.model();

    assert(site<2);

    const StateDescription::ISite &is=mi.origin()[site];
    const StateDescription::ISite &is2=mi.origin()[1-site];
    
    energy_t E_unp =
	    (mi.origin().num_sites()==1)
	    ?
	    model.energy_unpair(is)
	    :
	    model.energy_unpair(is, is2);
        	
    energy_t E_loop =
	model.energy_hybrid_loop(is.i1,is.i2,is.j1,is.j2);
    
    energy_t E_hyb = 0;
    if (mi.origin().num_sites()==2) {
	E_hyb = model.energy_hybrid(is2);	    
    }
    
    return E_unp + E_loop + E_hyb;
}

void
HybEnsModel::RemoveSiteMove::apply(StateDescription &sd) const {
    assert(site<2);
    
    if (site==0) {
	if (sd.num_sites()==2) {
	    sd[0]=sd[1];
	    
	}
	sd.resize(sd.num_sites()-1);
    } else {
	sd.resize(1);
    }
}


// ------------------------------------------------------------
// Move that creates a new site
//

HybEnsModel::NewSiteMove::NewSiteMove(const MoveIterator &mi) : Move(mi) {
}

HybEnsModel::NewSiteMove::~NewSiteMove() {}


std::ostream &
HybEnsModel::NewSiteMove::print(std::ostream &out) const {
    out << "[New Site Move "<<i1<<" "<<i2<<"]";
    return out;
}

// ------------------------------------------------------------
// Move that creates a new site as the first interaction site
//

HybEnsModel::NewSiteMoveF::NewSiteMoveF(const MoveIterator &mi) : NewSiteMove(mi) {
}

HybEnsModel::NewSiteMoveF::~NewSiteMoveF() {}



HybEnsModel::Move *
HybEnsModel::NewSiteMoveF::nextMoveType() const {
    return new NewSiteMoveL(mi);
}

bool
HybEnsModel::NewSiteMoveF::first() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();

    size_t len1 = model.seqA().length();
    
    if (o.num_sites()!=0) {
	return false;
    }
    
    i1=1;
    i2=0;
    
    return (i1 + model.minsitesize() < len1) 
	&& 
	next();
}

bool
HybEnsModel::NewSiteMoveF::next() {
    const HybEnsModel &model=mi.model();
    
    size_t len1 = model.seqA().length();
    size_t len2 = model.seqB().length();
    
    i2++;
    if (i2 + model.minsitesize() < len2) {
	return true;
    } else {
	i1++;
	i2=1;
	return (i1 + model.minsitesize()  < len1);
    }
}

HybEnsModel::energy_t
HybEnsModel::NewSiteMoveF::transitionEnergy() const {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    energy_t E_loop =
	model.energy_hybrid_loop(is_new.i1,is_new.i2,is_new.j1,is_new.j2);
    
    assert (o.num_sites()==0);
    
    energy_t E_unp =
	model.energy_unpair(is_new);
    
    return E_unp + E_loop;
    
}

void
HybEnsModel::NewSiteMoveF::apply(StateDescription &sd) const {
    const HybEnsModel &model=mi.model();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    assert (sd.num_sites()==0);
    
    sd.resize(1);
    sd[0]=is_new;
    
}


// ------------------------------------------------------------
// Move that creates a new site on the left
//

HybEnsModel::NewSiteMoveL::NewSiteMoveL(const MoveIterator &mi) : NewSiteMove(mi) {}

HybEnsModel::NewSiteMoveL::~NewSiteMoveL() {}

HybEnsModel::Move *
HybEnsModel::NewSiteMoveL::nextMoveType() const {
    return new NewSiteMoveR(mi);
}

bool
HybEnsModel::NewSiteMoveL::first() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    if (o.num_sites()!=1) {
	return false;
    }
    
    i1=1; 
    i2=0;
    
    return (i1 + model.minsitesize() + model.minsitedist() < o[0].i1) 
	&&
	( i2 + 1 + model.minsitesize() + model.minsitedist() < o[0].i2 )
	&& 
	next();
}

bool
HybEnsModel::NewSiteMoveL::next() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
        
    i2++;
    if (i2 + model.minsitesize() + model.minsitedist() < o[0].i2) {
	return true;
    } else {
	i1++;
	i2=1;
	return (i1 + model.minsitesize() + model.minsitedist() < o[0].i1);
    }
}

HybEnsModel::energy_t
HybEnsModel::NewSiteMoveL::transitionEnergy() const {
    const StateDescription &o=mi.origin();
    assert(o.num_sites()==1);
    
    const HybEnsModel &model=mi.model();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    const StateDescription::ISite &is=mi.origin()[0];
    
    energy_t E_loop =
	model.energy_hybrid_loop(is_new.i1,is_new.i2,is_new.j1,is_new.j2);
    
    energy_t E_unp =
	model.energy_unpair(is_new,is);
    
    energy_t E_hyb = model.energy_hybrid(is);	    
    
    return E_unp + E_loop + E_hyb;
}

void
HybEnsModel::NewSiteMoveL::apply(StateDescription &sd) const {
    const HybEnsModel &model=mi.model();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    const StateDescription::ISite is=sd[0];
    sd.resize(2);
    sd[0] = is_new;
    sd[1]=is;       
}

// ------------------------------------------------------------
// Move that creates a new site on the right
//

HybEnsModel::NewSiteMoveR::NewSiteMoveR(const MoveIterator &mi) : NewSiteMove(mi) {}

HybEnsModel::NewSiteMoveR::~NewSiteMoveR() {}

HybEnsModel::Move *
HybEnsModel::NewSiteMoveR::nextMoveType() const {
    return new MergeMove(mi);
}

bool
HybEnsModel::NewSiteMoveR::first() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();

    size_t len1 = model.seqA().length();
    size_t len2 = model.seqB().length();
    
    const StateDescription::ISite &is=o[0];
    
    if (o.num_sites()!=1) {
	return false;
    }
    
    i1=is.j1+model.minsitedist()+2; 
    i2=is.j2+model.minsitedist()+1;
    
    return (i1 + model.minsitesize() < len1) 
	&&
	( i2 + 1 + model.minsitesize() < len2 )
	&& 
	next();
}

bool
HybEnsModel::NewSiteMoveR::next() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    size_t len1 = model.seqA().length();
    size_t len2 = model.seqB().length();

    const StateDescription::ISite &is=o[0];
    
    i2++;
    if (i2 + model.minsitesize() < len2) {
	return true;
    } else {
	i1++;
	i2=is.j2+model.minsitedist()+2;
	
	return ( i1 + model.minsitesize() < len1 );
    }
}

HybEnsModel::energy_t
HybEnsModel::NewSiteMoveR::transitionEnergy() const {
    const StateDescription &o=mi.origin();
    assert(o.num_sites()==1);
    
    const HybEnsModel &model=mi.model();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    const StateDescription::ISite &is=mi.origin()[0];
    
    energy_t E_loop =
	model.energy_hybrid_loop(is_new.i1,is_new.i2,is_new.j1,is_new.j2);
    
    energy_t E_unp =
	model.energy_unpair(is,is_new);
    
    energy_t E_hyb = model.energy_hybrid(is);	    
    
    return E_unp + E_loop + E_hyb;
}

void
HybEnsModel::NewSiteMoveR::apply(StateDescription &sd) const {
    const HybEnsModel &model=mi.model();
    
    const StateDescription::ISite is_new=
	StateDescription::ISite(i1,i2,i1+model.minsitesize()-1,i2+model.minsitesize()-1);
    
    sd.resize(2);
    sd[1]=is_new;
}



// ------------------------------------------------------------
// Move that merges two sites
//

HybEnsModel::MergeMove::MergeMove(const MoveIterator &mi) : Move(mi) {}

HybEnsModel::MergeMove::~MergeMove() {}

HybEnsModel::Move *
HybEnsModel::MergeMove::nextMoveType() const {
    return new SplitMove(mi);
}

bool
HybEnsModel::MergeMove::first() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    if (o.num_sites()!=2) {
	return false;
    }
    
    return 
	( o[0].j1 + model.minsitedist() + 1 == o[1].i1 )
	&&
	( o[0].j2 + model.minsitedist() + 1 == o[1].i2 );
}

bool
HybEnsModel::MergeMove::next() {
    return false;
}

HybEnsModel::energy_t
HybEnsModel::MergeMove::transitionEnergy() const {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    
    const StateDescription::ISite is_between=
	StateDescription::ISite(o[0].j1,o[0].j2,o[1].i1,o[1].i2);

    const StateDescription::ISite is_merged=
	StateDescription::ISite(o[0].i1,o[0].i2,o[1].j1,o[1].j2);
    
    energy_t E_unp =
	model.energy_unpair(is_merged);
    
    energy_t E_hyb =
	model.energy_hybrid(is_between) 
	+ model.energy_hybrid(mi.origin()[0]) 
	+ model.energy_hybrid(mi.origin()[1]);
    
    return E_unp + E_hyb;
}

void
HybEnsModel::MergeMove::apply(StateDescription &sd) const {   
    const StateDescription::ISite is_merged=
	StateDescription::ISite(sd[0].i1,sd[0].i2,sd[1].j1,sd[1].j2);
    
    sd[0] = is_merged;
    sd.resize(1);
}

std::ostream &
HybEnsModel::MergeMove::print(std::ostream &out) const {
    out << "[Merge Move]";
    return out;
}

// ------------------------------------------------------------
// Move that splits one site into two
//

HybEnsModel::SplitMove::SplitMove(const MoveIterator &mi) : Move(mi) {}

HybEnsModel::SplitMove::~SplitMove() {}

std::ostream &
HybEnsModel::SplitMove::print(std::ostream &out) const {
    out << "[Split Move "<<i1<<" "<<i2<<"]";
    return out;
}


HybEnsModel::Move *
HybEnsModel::SplitMove::nextMoveType() const {
    return NULL;
}

bool
HybEnsModel::SplitMove::first() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    
    if (o.num_sites()!=1) {
	return false;
    }
    
    i1 = o[0].i1 + model.minsitesize() - 1;
    i2 = o[0].i2 + model.minsitesize() - 1;
    
    return 
	( i1 + model.minsitedist() + model.minsitesize() <= o[0].j1 )
	&&
	( i2 + model.minsitedist() + model.minsitesize() <= o[0].j2 );
}

bool
HybEnsModel::SplitMove::next() {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
        
    i2++;
    if (i2 + model.minsitedist() + model.minsitesize() <= o[0].j2) {
	return true;
    } else {
	i1++;
	i2=o[0].i2 + model.minsitesize() - 1;
	return (i1 + model.minsitedist() + model.minsitesize() <= o[0].j1);
    }
}

HybEnsModel::energy_t
HybEnsModel::SplitMove::transitionEnergy() const {
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
        
    const StateDescription::ISite &is_orig = o[0];

    const StateDescription::ISite is_between=
	StateDescription::ISite(i1,i2,
				i1+model.minsitedist()+1,
				i2+model.minsitedist()+1);

    const StateDescription::ISite is1=
	StateDescription::ISite(is_orig.i1,is_orig.i2,is_between.i1,is_between.i2);

    const StateDescription::ISite is2=
	StateDescription::ISite(is_between.j1,is_between.j2,is_orig.j1,is_orig.j2);
    
    energy_t E_unp =
	model.energy_unpair(is_orig);
    
    energy_t E_hyb =
	model.energy_hybrid(is_between) 
	+ model.energy_hybrid(is1) 
	+ model.energy_hybrid(is2);
    
    return E_unp + E_hyb;
}

void
HybEnsModel::SplitMove::apply(StateDescription &sd) const {
    
    const HybEnsModel &model=mi.model();
    const StateDescription &o=mi.origin();
    const StateDescription::ISite &is_orig = o[0];

    const StateDescription::ISite is_between=
	StateDescription::ISite(i1,i2,
				i1+model.minsitedist()+1,
				i2+model.minsitedist()+1);

    const StateDescription::ISite is1=
	StateDescription::ISite(is_orig.i1,is_orig.i2,is_between.i1,is_between.i2);

    const StateDescription::ISite is2=
	StateDescription::ISite(is_between.j1,is_between.j2,is_orig.j1,is_orig.j2);
    
    sd.resize(2);
    sd[0] = is1;
    sd[1] = is2;
}



// ------------------------------------------------------------
// Move Iterator


HybEnsModel::MoveIterator::
MoveIterator(const StateDescription &origin,
	     const HybEnsModel &model)
    :origin_(origin),
     model_(model)
{
}


HybEnsModel::Move *
HybEnsModel::MoveIterator::firstMove() const {
    Move *m = new GrowShrinkMoveFL(*this);
    
    while( (m!=NULL) && (! m->first()) ) {
	Move *m2 = m->nextMoveType();
	delete m;
	m=m2;
    }
    return m;
}


HybEnsModel::Move *
HybEnsModel::MoveIterator::nextMove(Move *m) const {
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
HybEnsModel::MoveIterator::disposeMove(Move *move) const {
    if (move!=NULL) delete move;
}
