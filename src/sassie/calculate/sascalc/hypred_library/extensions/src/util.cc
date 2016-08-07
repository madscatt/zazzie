#include "../include/util.h"

// the file parser class that serves the constructor of class with multiple constant variables
void
prdf::utils::
Error(const std::string & s)
{
	std::cerr<<s<<std::endl;
	exit(1);
};

/*
// substing replacement
std::string
prdf::utils::
replaceSubString(std::string s, const std::string & olds, const std::string & news)
{
	size_t idx=0;
	while (true)
	{
		idx = s.find(olds, idx);
		if (idx==std::string::npos) break;
		s.replace(idx, olds.length(), news);
		idx += news.length();
	}
	return s;
}

*/
