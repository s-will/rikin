#ifndef PFHYB_AUX_HH
#define PFHYB_AUX_HH

#include <vector>



/**
 * Defines an interface for a unsigned char vector that supports
 * encoding and decoding of sequences of integer numbers with given
 * number of bits per element
 *
 * \note This is currently not used, but an idea to support more flexible encoding/decoding of the states
 */
class CodeVector : public std::vector<unsigned char> {
    size_t bitsize;
    
public:
    CodeVector();
    
    void
    push_back(unsigned int x, size_t bits);
    
    void
    set(unsigned int x, size_t bitindex, size_t bits);
    
    unsigned int
    get(size_t bitindex, size_t bits) const;
    
};


#endif
