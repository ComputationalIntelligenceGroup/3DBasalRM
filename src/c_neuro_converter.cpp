#include <Rcpp.h>

#include <neurostr/core/log.h>
#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/io/JSONWriter.h>

#include <rapidjson/ostreamwrapper.h>
#include <rapidjson/writer.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <iostream>

#include <neurostr/selector/neuron_selector.h>

using namespace Rcpp;

namespace ns = neurostr::selector;
using namespace Rcpp;
namespace bg =        boost::geometry;
using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;
template <typename T>
using const_selector_reference = std::reference_wrapper<const T>;
using const_node_reference  = const_selector_reference<neurostr::Node>;

// [[Rcpp::export]]
std::string c_neuro_converter(std::string ifile, bool correct, float eps) {
  std::ostringstream oss;

  //Read the file DAT or JSON file
  auto r = neurostr::io::read_file_by_ext(ifile);

  //If there is at least one contour and there is not a soma it is asummed that the soma was saved as a contour
  if((r->n_contours() > 0) && (!r->begin()->has_soma()))
  {
    //Get the id of the last node
    std::vector<neurostr::Node> soma_nodes;
    int init_id = 0;

    //Get all the nodes in the neuron and search for the max id. It is needed to give an id to the soma
    std::vector<const_node_reference> neuron_nodes=ns::neuron_node_selector(*(r->begin()));
    neurostr::Node::id_type max_id_node = 0;
    for(auto it = neuron_nodes.begin(); it != neuron_nodes.end(); ++it)
    {
      if(max_id_node < it->get().id())
      {
        max_id_node = it->get().id();
      }
    }
    init_id = max_id_node;

    auto std_contour = r->contour_begin();
    for(auto it = std_contour->begin(); it != std_contour->end(); ++it)
    {

      neurostr::Node new_node(init_id, *it, 0.2);
      soma_nodes.push_back(new_node);
      init_id++;
    }

    // Simpify, correct and add the nodes in the contour to the soma
    for(auto it = r->begin(); it != r->end(); ++it){
      if(correct) it->correct();
      if(eps != 0.0 ){
        it->simplify(eps);
      }
      it->erase_apical();
      it->erase_axon();
      it->add_soma(soma_nodes);
      it->center();
    }

  }else{ //If the file does not contain contours or have a defined soma

    for(auto it = r->begin(); it != r->end(); ++it){
      if(correct) it->correct();
      if(eps != 0.0 ){
        it->simplify(eps);
      }
    }
  }

  r->erase_contour();
  neurostr::io::JSONWriter writer(oss);
  writer.write(*r);

  std::string features = oss.str();
  return features;
}
