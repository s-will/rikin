#include "pair_pfs.hh"

#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include <zlib.h>


PairPfs::PairPfs(const std::string &filename) {
    // note that gzopen recognizes uncompressed files; if the input is
    // uncompressed, gzread will directly return the uncompressed input.
    gzFile fh=gzopen(filename.c_str(),"r");
    if (fh==NULL) {
	std::cerr << "Cannot read basin ipps from (optionally compressed) file "
		  << filename<<"."<<std::endl;
    } else {
	gzread_pps(fh);
    }
    gzclose(fh);
}

gzFile
PairPfs::gzread_pps(gzFile fh) {
    std::string line;
    size_t bufsiz=512; // any buffer size should work, since we read longer lines in chunks
    char buf[bufsiz];
    
    // read line by line
    do {
	line="";
	// read one line (in chunks of at most bufsiz)
	do {
	    // TO CHECK: what happens if last line in file is not terminated by new line?
	    if (gzgets(fh,buf,bufsiz) == NULL) { return NULL; }
	    line+=buf;
	    if (gzeof(fh)) break;
	} while(line.length()==0 || line[line.length()-1]!='\n');
	
	// remove trailing newline
	if (line.length()>0 && line[line.length()-1]=='\n') {
	    line.pop_back();
	}
	
	// parse non-empty lines
	if (line.length()>0) {
	    size_t state_idx;
	    double state_pf;
	    pp_matrix_t pps;
	    std::istringstream linein(line);
	    read_state(linein,state_idx,state_pf,pps);
	    if (pps_.size()<=state_idx) {
		pps_.resize(state_idx+1);
	    }
	    if (Z_.size()<=state_idx) {
		Z_.resize(state_idx+1);
	    }
	    pps_[state_idx] = pps;
	    Z_[state_idx]   = state_pf;
	}
	
    } while (!gzeof(fh));

    return fh;
}

std::ostream &
PairPfs::write_pps(std::ostream &out, double min_prob) const {

    for ( size_t idx=0; idx<pps_.size(); idx++ ) {
	if (!pps_[idx].empty()) {
	    PairPfs::write_state(out,idx,Z_[idx],pps_[idx],min_prob,false);
	}
    }

    return out;
}

gzFile
PairPfs::gzwrite_pps(gzFile fh, double min_prob) const {

    for ( size_t idx=0; idx<pps_.size(); idx++ ) {
	if (!pps_[idx].empty()) {
	    std::ostringstream out;
	    PairPfs::write_state(out,idx,Z_[idx],pps_[idx],min_prob,false);
	    gzputs(fh,out.str().c_str());
	}
    }

    return fh;
}

std::istream &
PairPfs::read_state(std::istream &in,
		    size_t &state_idx,
		    double &state_pf,
		    pp_matrix_t &pps) {
    
    pps.clear();
    in >> state_idx;
    in >> state_pf;

    std::string line;
    if (!getline(in,line)) return in;
    
    std::istringstream linein(line);
    
    size_t i;
    size_t j;
    double pij;
    while(linein >> i >> j >> pij) {
	pps(i,j) = pij;
    }
    
    return in;
}
