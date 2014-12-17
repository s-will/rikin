#include <basin.hh>
#include <cstdio>
#include <cmath>

/* Methods of basin */

void Basin::print_header(std::ostream &out) const {
    printf("%5s %10s %6s",
	   "idx",
	   "n_s",
	   "Z"
	   );
}

void
Basin::print(std::ostream &out) const {
    printf("%5lu %10.2f %6.2f",
	   index_,
	   states_,
	   log(Z_)
	   );
}
