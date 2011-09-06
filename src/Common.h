#ifndef transabs_Common_h
#define transabs_Common_h

#include <boost/program_options.hpp>

#include <armadillo>

using namespace arma;
namespace po = boost::program_options;

extern po::variables_map vm;

extern double soft_core;

#endif