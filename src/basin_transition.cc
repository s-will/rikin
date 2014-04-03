#include <basin_transition.hh>
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
	   states,
	   log(Z)
	   );
}

// BasinTransition::BasinTransition(std::istream &in) {
//     in.read(reinterpret_cast<char *>(&Z_),sizeof(double));
// }

// std::ostream &
// BasinTransition::write_binary(std::ostream &out) {
//     out.write(reinterpret_cast<char *>(&Z_),sizeof(double));
//     return out;
// }
