// code guard
#if !defined(prdf_utils_h)
#define prdf_utils_h

#include <iostream>
#include <string.h>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <dirent.h>
#include <math.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <stdlib.h>

#include "common.h"


// place everything in the local namespace
namespace prdf {
    namespace utils {
				
				// error handling
				void Error(const std::string &);
				/*
				// substing replacement
				std::string replaceSubString(std::string s, const std::string & olds, const std::string & news);
				// line parser
				class LineParser;
				// convert string to a template type
				template <typename T> T convertStringTo(const std::string &);
				*/

	} // namespace utils 
} // namespace prdf

// get the inline definitions
#define prdf_utils_icc
#include "util.icc"
#undef prdf_utils_icc

#endif
