#include <Rcpp.h>

#include <neurostr/core/neuron.h>
#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/methods/boxCutter.h>
#include <neurostr/io/JSONWriter.h>
#include <neurostr/io/JSONParser.h>
#include <rapidjson/document.h>


using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]

//This function cuts a neuron reconstruction according to the boundaries stipulated by the user
// [[Rcpp::export]]
std::string c_box_cutter(std::string json_info, NumericVector mins, NumericVector maxs) {

  std::ifstream aux;//Auxiliar variable needed to use the constructor of the neurostr library. Does not have meaning for the context of this package
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());

  //Read the JSON file and return a reconstruction
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);
  auto r = p->read_string("reconstruction",doc);

  delete p;
  std::ostringstream oss;

  //Get the boundaries defined by the user
  std::vector<float> min_corner = Rcpp::as<std::vector<float> >(mins);
  std::vector<float> max_corner = Rcpp::as<std::vector<float>>(maxs);

  // For each neuron
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){

    // Cut every neuron in the reconstruction
    neurostr::Neuron& n = *n_it;

    // Get the values of the bounding box for the complete neuron
    auto box = n.boundingBox();
    neurostr::point_type original_min_corner = box.min_corner();
    neurostr::point_type original_max_corner = box.max_corner();

    // If no value was given to some dimensions replace by the boundary box value
    if(!R_FINITE(mins[0])) min_corner[0] = neurostr::geometry::getx(original_min_corner)-1;
    if(!R_FINITE(mins[1])) min_corner[1] = neurostr::geometry::gety(original_min_corner)-1;
    if(!R_FINITE(mins[2])) min_corner[2] = neurostr::geometry::getz(original_min_corner)-1;

    if(!R_FINITE(maxs[0])) max_corner[0] = neurostr::geometry::getx(original_max_corner)+1;
    if(!R_FINITE(maxs[1])) max_corner[1] = neurostr::geometry::gety(original_max_corner)+1;
    if(!R_FINITE(maxs[2])) max_corner[2] = neurostr::geometry::getz(original_max_corner)+1;

    // Cut neuron
    neurostr::methods::neuronBoxCutter(n,
                                       neurostr::point_type(min_corner[0],min_corner[1],min_corner[2]),
                                       neurostr::point_type(max_corner[0],max_corner[1],max_corner[2]));

  } // End neuron for
  neurostr::io::JSONWriter writer(oss);
  writer.write(*r);

  std::string cut_neuron = oss.str();

  return(cut_neuron);
}
