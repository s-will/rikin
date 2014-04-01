#include <stdlib.h>
#include <math.h>

#include  <stdlib.h>
#include  <stdio.h>


#include "hybrid_ensemble_model.hh"

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
HybEnsModel::StateDescription::size() const {
    return isites.size();
}

std::string
HybEnsModel::StateDescription::to_string() const {
    std::ostringstream out;
    out << (*this);
    return out.str();
}

size_t
HybEnsModel::StateDescription::max_site_size() const {
    size_t mss=0;
    for(std::vector<ISite>::const_iterator it=isites.begin(); it!=isites.end(); ++it) {
	mss = std::max(mss,it->j1-it->i1+1);
	mss = std::max(mss,it->j2-it->i2+1);
    }
    return mss;
}


// --------------------
// encoding and decoding compressed representation
//
// use very simple code:
// always use 2 bytes per position and 8 * 2 bytes in total (two sites)
//
// encode empty sites by invalid entry where i1 > j1 to avoid 0 in
// encoding; as long as positions are >0 this allows to write encoding
// to file and still use gnu's "sort -z"! ---  CURRENTLY BROKEN
//

HybEnsModel::StateDescription::code_t &
HybEnsModel::StateDescription::encode(code_t &the_code) const {
    
    for (size_t i=0; i<size(); i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((size()==0)?(&the_code.first):(&the_code.second));
    
	*(code++)=isites[i].i1;
	*(code++)=isites[i].i2;
	*(code++)=isites[i].j1;
	*(code++)=isites[i].j2;
    }
    
    // if site()!=2, mark the second site or both sites as empty
    for (size_t i=0; i<(2-size());i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((size()==1)?(&the_code.first):(&the_code.second));

	*(code++)=2;
	*(code++)=2;
	*(code++)=1;
	*(code++)=1;
    }
    
    return the_code;
}

HybEnsModel::StateDescription::code_t
HybEnsModel::StateDescription::encode() const {

    code_t code;
    
    encode(code);
    
    return code;
}



HybEnsModel::StateDescription &
HybEnsModel::StateDescription::decode(const code_t &the_code) {
    size_t num_sites=0;
    for (size_t i=0; i<2; i++) {
	const unsigned short *code = 
	    reinterpret_cast<const unsigned short *>((size()==0)?(&the_code.first):(&the_code.second));
	if (code[0]<=code[2]) num_sites++;
    }
    
    isites.reserve(num_sites);
    isites.resize(0);
    
    for (size_t i=0; i<num_sites; i++) {
	const unsigned short *code = 
	    reinterpret_cast<const unsigned short *>((size()==0)?(&the_code.first):(&the_code.second));
	isites.push_back(ISite(*(code),*(code+1),*(code+2),*(code+3)));
    }
    return *this;
}

void
HybEnsModel::StateDescription::write_binary() const {
    code_t code;
    encode(code);
    fwrite(reinterpret_cast<char *>(&code.first),sizeof(char),8,stdout);
    fwrite(reinterpret_cast<char *>(&code.second),sizeof(char),8,stdout);
    fputc(0,stdout);
}

bool
HybEnsModel::StateDescription::read_binary(std::istream &in) {
    HybEnsModel::StateDescription::code_t code;
    
    char *codebuf;
    codebuf = reinterpret_cast<char *>(&code.first);
    in.read(codebuf,8);
    codebuf = reinterpret_cast<char *>(&code.second);
    in.read(codebuf,8);
        
    if(in.get()!=0) {
	return false;
    }
    
    decode(code);
    return true;
}

bool
HybEnsModel::StateDescription::is_valid(const HybEnsModel &model) const {
    bool valid=true;
    
    // at most 2 sites
    valid = valid && size()<=2;
	    
    // check site size and boundaries
    for (size_t i=0; valid && i<size(); ++i) {
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
    if (size()==2) {
	valid = valid
	    && isites[0].j1 + model.minsitedist() + 1 <= isites[1].i1;
	valid = valid
	    && isites[0].j2 + model.minsitedist() + 1 <= isites[1].i2;
    }

    return valid;
}


// ------------------------------------------------------------
// Implementation of HybEnsModel

HybEnsModel::HybEnsModel(std::string seqA,
			 std::string seqB,
			 size_t maxsitesize,
			 size_t maxsitesize_diff,
			 bool cond
			 )
    : seqA_(seqA),
      seqB_(seqB),
      uppfA_(seqA,+1,maxsitesize,cond),
      uppfB_(seqB,-1,maxsitesize,cond),
      hybridpf_(seqA,seqB,maxsitesize,maxsitesize_diff),
      maxunpinloop_(6),
      minsitesize_(3),
      minsitedist_(6),
      maxsitesize_( maxsitesize ),
      homodimer_(false)
{
    std::string seqA1=seqA;
    reverse(seqA1);
    if (seqA1==seqB) {
	const_cast<bool&>(homodimer_)=true;
    }
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
    switch(sd.size()) {
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
    for (size_t i=0; i<sd.size(); i++) {
	out << sd[i];
	if (i<sd.size()-1) {out << ",";}
    }
    out << "}";
    return out;
}

