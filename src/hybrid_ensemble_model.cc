#include <cmath>
#include <cstdint>

//#include  <cstdlib>
//#include  <cstdio>
#include <iostream>
#include <iomanip>

#include "hybrid_ensemble_model.hh"

// ----------------------------------------
// HybEnsModel::StateDescription
//

HybEnsModel::StateDescription::StateDescription() 
{
}

/**
 * @brief Construct state from binary code
 */
HybEnsModel::StateDescription::StateDescription(const code_t &code) 
{
    decode(code);
}


HybEnsModel::StateDescription::StateDescription(size_t i1, size_t i2, size_t j1, size_t j2)
{
    isites.push_back(ISite(i1,i2,j1,j2));
}



HybEnsModel::StateDescription::StateDescription(size_t i1, size_t i2, size_t j1, size_t j2,
						size_t k1, size_t k2, size_t l1, size_t l2)
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
namespace {

    // Pack four 16-bit fields (i1,i2,j1,j2) into one 64-bit word.
    inline uint64_t
    pack_fields(uint16_t i1, uint16_t i2, uint16_t j1, uint16_t j2) {
	return (uint64_t(i1) << 48)
	     | (uint64_t(i2) << 32)
	     | (uint64_t(j1) << 16)
	     |  uint64_t(j2);
    }

    // Unpack a 64-bit word into four 16-bit fields (i1,i2,j1,j2).
    inline void
    unpack_fields(uint64_t word, uint16_t &i1, uint16_t &i2, uint16_t &j1, uint16_t &j2) {
	i1 = static_cast<uint16_t>((word >> 48) & 0xFFFFu);
	i2 = static_cast<uint16_t>((word >> 32) & 0xFFFFu);
	j1 = static_cast<uint16_t>((word >> 16) & 0xFFFFu);
	j2 = static_cast<uint16_t>( word         & 0xFFFFu);
    }

    // Convert a sequence position to the 16-bit field used by code_t,
    // asserting that no information is lost (instead of silently
    // truncating, as a raw static_cast to unsigned short would).
    inline uint16_t
    to_field(size_t x) {
	assert(x <= 0xFFFFu && "position exceeds 16-bit code_t encoding range");
	return static_cast<uint16_t>(x);
    }

    // To avoid 0x00 bytes in the binary file output, each 16-bit field
    // is spread over two bytes whose most significant bit is always 1:
    //   x = 00aaaaaa abbbbbbb  (only bits 0..13 of x are used)
    //   -> stored as bit pattern  aaaaaaa1 bbbbbbb1
    // IMPORTANT: this scheme only preserves 14 bits. Values >= 0x4000
    // (16384) would previously be truncated *silently*; here that is
    // turned into an assertion so it fails loudly during development
    // instead of corrupting positions for large sequences.
    inline void
    encode_ushort(uint16_t &x) {
	assert(x < 0x4000u && "encode_ushort: value exceeds 14-bit range, would be truncated");
	x = 1 | ( (x & 0x7F) << 1 ) | 0x100 | ( ((x>>7) & 0x7F) << 9 );
    }

    inline void
    decode_ushort(uint16_t &x) {
	x = ((x>>1) & 0x7F) | (((x>>9) & 0x7F)<<7);
    }

} // end anonymous namespace

HybEnsModel::StateDescription::code_t &
HybEnsModel::StateDescription::encode(code_t &the_code) const {

    assert(size()<=2);

    uint64_t *words[2] = { &the_code.first, &the_code.second };

    for (size_t i=0; i<2; i++) {
	if (i<size()) {
	    const ISite &is = isites[i];
	    *words[i] = pack_fields(to_field(is.i1), to_field(is.i2),
				     to_field(is.j1), to_field(is.j2));
	} else {
	    // mark unused site slot as empty: i1=2,j1=1 so that i1>j1,
	    // which decode() recognizes as an invalid/empty site.
	    *words[i] = pack_fields(2,2,1,1);
	}
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

    assert(the_code.first!=0 and the_code.second!=0);

    const uint64_t words[2] = { the_code.first, the_code.second };

    // decode and determine number of encoded sites
    size_t num_sites=0;
    for (size_t i=0; i<2; i++) {
	uint16_t i1,i2,j1,j2;
	unpack_fields(words[i], i1, i2, j1, j2);

	if (i1<=j1) { num_sites++; }
	else { break; } // break at first invalid site
    }

    isites.clear();
    isites.reserve(num_sites);

    for (size_t i=0; i<num_sites; i++) {
	uint16_t i1,i2,j1,j2;
	unpack_fields(words[i], i1, i2, j1, j2);
	isites.push_back(ISite(i1,i2,j1,j2));
    }
    return *this;
}

std::ostream &
HybEnsModel::StateDescription::write_binary(std::ostream &out) const {
    code_t the_code;
    encode(the_code);

    const uint64_t words[2] = { the_code.first, the_code.second };

    for (uint64_t word : words) {
	uint16_t fields[4];
	unpack_fields(word, fields[0], fields[1], fields[2], fields[3]);

	for (uint16_t &f : fields) {
	    encode_ushort(f);
	}

	// Serialize explicitly as little-endian bytes. This avoids
	// aliasing a uint16_t through the raw storage of an unrelated
	// object and makes the on-disk format independent of the host's
	// native endianness/short size.
	char buf[8];
	for (size_t k=0; k<4; k++) {
	    buf[2*k]   = static_cast<char>( fields[k]       & 0xFF);
	    buf[2*k+1] = static_cast<char>((fields[k] >> 8) & 0xFF);
	}
	out.write(buf, sizeof(buf));
    }

    return out;
}

bool
HybEnsModel::StateDescription::read_binary(std::istream &in) {
    uint64_t words[2];

    for (uint64_t &word : words) {
	char buf[8];
	in.read(buf, sizeof(buf));
	if (!in) { return false; }

	uint16_t fields[4];
	for (size_t k=0; k<4; k++) {
	    fields[k] = static_cast<uint16_t>(static_cast<unsigned char>(buf[2*k]))
	              | static_cast<uint16_t>(static_cast<unsigned char>(buf[2*k+1]) << 8);
	    decode_ushort(fields[k]);
	}

	word = pack_fields(fields[0], fields[1], fields[2], fields[3]);
    }

    HybEnsModel::StateDescription::code_t the_code(words[0], words[1]);
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
			 bool model_double_sites
			 )
    : 
    seqA_(seqA),
    seqB_(seqB),
    uppfA_(seqA,+1,maxsitesize,span,window,model_double_sites),
    uppfB_(seqB,-1,maxsitesize,span,window,model_double_sites),
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
    assert(region_startA_<=i1);
    assert(region_startB_<=i2);
    assert(i1<=j1);
    assert(i2<=j2);
    assert(j1<=region_endA_);
    assert(j2<=region_endB_);
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
        return 0;
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

