#include <basin.hh>
#include <iostream>
#include <iomanip>
#include <cmath>

/* Methods of basin */

void Basin::print_header(std::ostream &out) const {
    out <<std::setw(5)<<"idx"<<" "
	<<std::setw(10)<<"n_s"<<" "
	<<std::setw(6)<<"Z";
}

void
Basin::print(std::ostream &out) const {
    out <<std::setw(5)<<index_<<" "
	<<std::setw(10)<<states_<<" "
	<<std::setw(6)<<std::setprecision(2)<<log(Z_);
}
