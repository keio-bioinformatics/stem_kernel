// $Id:$

#ifndef __INC_FA_H__
#define __INC_FA_H__

#include <boost/spirit/iterator/file_iterator.hpp>
#include <list>
#include "../common/rna.h"

template < class Seq >
bool
load_fa(Seq& s, boost::spirit::file_iterator<>& fi);

template < class Seq >
bool
load_fa(std::list<Seq>& ma, const char* filename);

template < class Seq >
bool
load_fa(std::list< MASequence<Seq> >& ma, const char* filename);

#endif //  __INC_FA_H__

// Local Variables:
// mode: C++
// End:
