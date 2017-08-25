#include <Rcpp.h>


#include <neurostr/core/log.h>
#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/io/JSONWriter.h>

#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <iostream>

using namespace Rcpp;

// [[Rcpp::export]]
std::string c_neuro_converter2(std::string ifile, bool correct, float eps) {
  std::ostringstream oss;

  auto r = neurostr::io::read_file_by_ext(ifile);

  // Simpify / correct
  for(auto it = r->begin(); it != r->end(); ++it){
    if(correct) it->correct();
    if(eps != 0.0 ){
      it->simplify(eps);
    }
  }

  neurostr::io::JSONWriter writer(oss);
  writer.write(*r);

  std::string features = oss.str();

  return features;
}
