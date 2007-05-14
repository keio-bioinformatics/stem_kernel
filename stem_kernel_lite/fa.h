// $Id:$

#ifndef __INC_FA_H__
#define __INC_FA_H__

#include <list>
#include <boost/spirit/iterator/file_iterator.hpp>
#include "../common/rna.h"

template < class Seq >
bool
load_fa(Seq& s, boost::spirit::file_iterator<>& fi);

template < class Seq >
bool
load_fa(std::list<Seq>& ma, const char* filename);

#endif //  __INC_FA_H__

// Local Variables:
// mode: C++
// End:
