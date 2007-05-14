// $Id$

#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "example.h"
#include "fasta.h"

uint
load_examples(std::ifstream& in, ExampleSet& ex)
{
  std::string buf;
  while (std::getline(in, buf)) {
    std::istringstream s(buf);
    std::string c;
    std::string seq;
    s >> c >> seq;
    if (!seq.empty()) {
      boost::algorithm::to_lower(seq);
      ex.push_back(Example(c,seq));
    }
  }
  return ex.size();
}

uint
load_examples(const std::string& label, Fasta& f, ExampleSet& ex)
{
  while (!f.IsEOF() && f.GetNextSeq()) {
    f.ToLower();
    ex.push_back(Example(label,f.Seq()));
  }
  return ex.size();
}

