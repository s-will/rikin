#include <iostream>
#include <sstream>
#include "util.hh"

int
str_to_int(const std::string &s) {
    std::istringstream ss(s);
    int x;
    ss>>x;
    return x;
}

void
parse_region(const std::string &s, size_t &start, size_t &end) {
    size_t dash = s.find('-');
    
    if (dash>0 && dash<s.length()) {
	start= str_to_int(s.substr(0,dash));
	end  = str_to_int(s.substr(dash+1));
	return;
    }
    
    int x=str_to_int(s);
    if (x>0) {
	end = std::min(end,(size_t)x);
    } else if (x<0) {
	if ((size_t)(-x)<end) {
	    start = end+x+1;
	}
    }
}
