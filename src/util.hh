#ifndef UTIL_HH
#define UTIL_HH

#include <string>

/** 
 * @brief Convert string to int
 * @param s string
 * @return integer
 */
int
str_to_int(const std::string &s);


/** 
 * @brief Parse region description
 * @param s string
 * @param[out] start start of region
 * @param[out] end end of region
 */
void
parse_region(const std::string &s, size_t &start, size_t &end);

#endif // UTIL_HH
