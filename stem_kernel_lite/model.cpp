// $Id$

#include "model.h"
#include <list>
#include <string>
#include <sstream>
#include <boost/spirit.hpp>
#include <boost/lambda/lambda.hpp>

using namespace boost::spirit;
using namespace boost::lambda;

struct model_parser : public grammar< model_parser >
{
  model_parser(std::list<uint>& l) : sv_list(l) { }

  std::list<uint>& sv_list;

  template < class ScannerT >
  struct definition
  {
    typedef rule<ScannerT> rule_t;
    rule_t model;
    rule_t line;
    rule_t sv_start;
    rule_t sv;

    definition(const model_parser& self)
    {
      model = *(line-sv_start) >> sv_start >> *sv >> end_p;
      sv_start = str_p("SV") >> eol_p;
      line = *print_p >> eol_p;
      sv = real_p >> +blank_p >> str_p("0:")
		  >> uint_p[push_back_a(self.sv_list)] 
		  >> +blank_p >> eol_p;
    }

    const rule_t& start() const { return model; }
  };
};

static
bool
load_sv_index(std::vector<uint>& sv_index, const char* filename)
{
  file_iterator<> fi(filename);
  if (!fi) {
    std::ostringstream os;
    os << filename << ": no such file";
    throw os.str().c_str();
    //return false;
  }
  std::list<uint> sv_list;
  model_parser parser(sv_list);
  parse_info<file_iterator<> > info = parse(fi, fi.make_end(), parser);
  if (!info.hit) {
    std::ostringstream os;
    os << filename << ": bad format";
    throw os.str().c_str();
    //return false;
  }

  sv_index.resize(sv_list.size());
  std::transform(sv_list.begin(), sv_list.end(), sv_index.begin(), _1-1);

  return true;
}

bool
load_sv_index(std::vector<uint>& sv_index,
	      const std::vector<std::string>& models)
{
  std::vector<std::string>::const_iterator x;
  bool res=true;
  for (x=models.begin(); x!=models.end(); ++x) {
    std::vector<uint> s;
    res=load_sv_index(s, x->c_str());
    if (!res) return res;

    if (sv_index.empty()) {
      sv_index.swap(s);
    } else {
      sv_index.insert(sv_index.end(), s.begin(), s.end());
    }
  }

  std::sort(sv_index.begin(), sv_index.end());
  std::vector<uint>::iterator e;
  e = std::unique(sv_index.begin(), sv_index.end());
  sv_index.erase(e, sv_index.end());
  
  return res;
}

#if 0
#include <iterator>

int
main(int argc, char* argv[])
{
  std::vector<std::string> models;
  for (uint i=1; i!=argc; ++i) {
    models.push_back(std::string(argv[i]));
  }
  
  std::vector<uint> sv_index;
  std::cout << load_sv_index(sv_index, models) << std::endl;
  std::copy(sv_index.begin(), sv_index.end(),
	    std::ostream_iterator<uint>(std::cout, " "));
  std::cout << std::endl;
  return 0;
}
#endif
