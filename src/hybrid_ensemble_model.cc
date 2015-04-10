#include <cmath>

//#include  <cstdlib>
//#include  <cstdio>
#include <iostream>
#include <iomanip>

#include "hybrid_ensemble_model.hh"

// ----------------------------------------
// HybEnsModel::StateDescription
//

HybEnsModel::StateDescription::StateDescription() 
    :isites(0)
{
}

/**
 * @brief Construct state from binary code
 */
HybEnsModel::StateDescription::StateDescription(const code_t &code) 
    :isites(0)
{
    decode(code);
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

std::string
HybEnsModel::to_dotbracket(const StateDescription &sd) const {
    size_t n = seqA().length();
    size_t m = seqB().length();
    
    std::string s(n+m+1,'.');
    
    s[n] = '&';
    
    for (size_t k=0; k<sd.size(); k++) {
	s[sd[k].i1-1]='(';
	s[sd[k].j1-1]=')';
	s[n+1+sd[k].i2-1]='(';
	s[n+1+sd[k].j2-1]=')';
    }
    
    return s;
}


// --------------------
// encoding and decoding compressed representation
//
// use simple code:
// always use 2 bytes per position and 8 * 2 bytes in total (two sites)
//
// encode empty sites by invalid entry where i1 > j1
//
HybEnsModel::StateDescription::code_t &
HybEnsModel::StateDescription::encode(code_t &the_code) const {
    
    for (size_t i=0; i<size(); i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((i==0)?(&the_code.first):(&the_code.second));
	
	*(code++)=isites[i].i1;
	*(code++)=isites[i].i2;
	*(code++)=isites[i].j1;
	*(code++)=isites[i].j2;
    }
    
    // if site()!=2, mark the second site or both sites as empty
    for (size_t i=size(); i<2; i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((i==0)?(&the_code.first):(&the_code.second));

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

    // decode and determine number of encoded sites
    size_t num_sites=0;
    for (size_t i=0; i<2; i++) {
        const unsigned short *code = 
	    reinterpret_cast<const unsigned short *>((i==0)?(&the_code.first):(&the_code.second));
		
	if (code[0]<=code[2]) { num_sites++; } 
	else { break; } // break at first invalid site
    }
    
    isites.reserve(num_sites);
    isites.resize(0);
    
    for (size_t i=0; i<num_sites; i++) {
	const unsigned short *code = 
	    reinterpret_cast<const unsigned short *>((i==0)?(&the_code.first):(&the_code.second));
	isites.push_back(ISite(*(code),*(code+1),*(code+2),*(code+3)));
    }
    return *this;
}


// to avoid 0 bytes in binary file output: encode unsigned short int
// x=00aaaaaaabbbbbbb as bit sequence aaaaaaa1bbbbbbb1
void
encode_ushort(unsigned short &x) {
    x = 1 | ( (x & 0x7F) << 1 ) | 0x100 | ( ((x>>7) & 0x7F) << 9 );
}
void
decode_ushort(unsigned short &x) {
    x = ((x>>1) & 0x7F) | (((x>>9) & 0x7F)<<7);
}

std::ostream &
HybEnsModel::StateDescription::write_binary(std::ostream &out) const {
    code_t the_code;
    encode(the_code);

    for (size_t i=0; i<2; i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((i==0)?(&the_code.first):(&the_code.second));
		
	encode_ushort(code[0]);
	encode_ushort(code[1]);
	encode_ushort(code[2]);
	encode_ushort(code[3]);
	
	out.write(reinterpret_cast<char *>(&code),4*sizeof(unsigned short));
    }
    
    return out;
}

bool
HybEnsModel::StateDescription::read_binary(std::istream &in) {
    HybEnsModel::StateDescription::code_t the_code;
    
    for (size_t i=0; i<2; i++) {
	unsigned short *code = 
	    reinterpret_cast<unsigned short *>((i==0)?(&the_code.first):(&the_code.second));
	
	in.read(reinterpret_cast<char *>(&code),4*sizeof(unsigned short));
	
	decode_ushort(code[0]);
	decode_ushort(code[1]);
	decode_ushort(code[2]);
	decode_ushort(code[3]);

    }

    decode(the_code);

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
			 size_t region_startA,
			 size_t region_endA,
			 size_t region_startB,
			 size_t region_endB,
			 size_t span,
			 size_t window,
			 bool cond
			 )
    : 
    seqA_(seqA),
    seqB_(seqB),
    uppfA_(seqA,+1,maxsitesize,span,window,cond),
    uppfB_(seqB,-1,maxsitesize,span,window,cond),
    hybrid_pf_(seqA,seqB,
	      maxsitesize,
	      maxsitesize_diff,
	      region_startA,
	      region_endA,
	      region_startB,
	      region_endB),
    maxunpinloop_(6),
    minsitesize_(3),
    minsitedist_(6),
    maxsitesize_( maxsitesize ),
    maxsitesize_diff_( maxsitesize_diff ),
    region_startA_(region_startA),
    region_endA_(region_endA),
    region_startB_(region_startB),
    region_endB_(region_endB),
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
    return - hybrid_pf_.RT() * log( hybrid_pf_.partition_function(i1,j1,i2,j2) );
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

double
HybEnsModel::interaction_probability(size_t k1, 
				     size_t k2,
				     const StateDescription &sd) const {
    assert(pair_type(k1,k2)>0);
    
    double prob=0;
    for ( auto &isite : sd ) {
	if (isite.i1<=k1 && k1<=isite.j1
	    &&
	    isite.i2<=k2 && k2<=isite.j2) {
	    
	    prob =
		hybrid_pf_.partition_function(isite.i1,k1,isite.i2,k2) *
		hybrid_pf_.partition_function(k1,isite.j1,k2,isite.j2)
		/ hybrid_pf_.partition_function(isite.i1,isite.j1,isite.i2,isite.j2);
	}
    }
    assert(0.0<=prob);
    if (prob>1.0) {
	std::cerr <<sd<<" "<<k1<<" "<<k2<<" ";
	for ( auto &isite : sd ) {
	    if (isite.i1<=k1 && k1<=isite.j1
		&&
		isite.i2<=k2 && k2<=isite.j2) {
		std::cerr << hybrid_pf_.partition_function(isite.i1,k1,isite.i2,k2) << "*" << hybrid_pf_.partition_function(k1,isite.j1,k2,isite.j2) << "/" << hybrid_pf_.partition_function(isite.i1,isite.j1,isite.i2,isite.j2) <<" = "<< prob << " ";
	    }
	}
	std::cout << std::endl;
    }
    assert(prob<=1.0);
    return prob;
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

