// $Id$

#ifndef __INC_ALN_H__
#define __INC_ALN_H__

#include <boost/spirit/iterator/file_iterator.hpp>
#include <list>
#include "../common/rna.h"

template < class Seq >
bool
load_aln(MASequence<Seq>& ma, boost::spirit::file_iterator<>& fi);

template < class Seq >
bool
load_aln(std::list< MASequence<Seq> >& ma, const char* filename);

#endif	// __INC_ALN_H__

// Local Variables:
// mode: C++
// End:
