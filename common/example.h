// $Id$

#ifndef __INC_EXAMPLE_H__
#define __INC_EXAMPLE_H__

#include <fstream>
#include <vector>
#include <string>
#include <utility>
#include "fasta.h"

typedef std::pair<std::string,std::string> Example;
typedef std::vector< Example > ExampleSet;

uint
load_examples(std::ifstream& in, ExampleSet& ex);

uint
load_examples(const std::string& label, Fasta& f, ExampleSet& ex);

#endif
