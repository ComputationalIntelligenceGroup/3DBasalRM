#include <Rcpp.h>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>

#include <neurostr/io/parser_dispatcher.h>
#include <neurostr/core/log.h>
#include <neurostr/core/neuron.h>
#include <neurostr/core/node.h>
#include <neurostr/core/contour.h>

#include <neurostr/selector/selector.h>
#include <neurostr/selector/universal_selector.h>
#include <neurostr/selector/neurite_selector.h>
#include <neurostr/selector/node_selector.h>

#include <neurostr/measure/universal_measure.h>
#include <neurostr/measure/branch_measure.h>
#include <neurostr/measure/neurite_measure.h>
#include <neurostr/measure/node_measure.h>
#include <neurostr/measure/aggregate.h>
#include <neurostr/measure/measure_operations.h>

namespace ns = neurostr::selector;
namespace nm = neurostr::measure;
namespace bg = boost::geometry;

using point_type =    bg::model::point<float, 3, bg::cs::cartesian>;
using namespace Rcpp;

// [[Rcpp::plugins(cpp14)]]
template <typename T>
using const_selector_reference = std::reference_wrapper<const T>;
using const_node_reference  = const_selector_reference<neurostr::Node>;

//Get the coordinates of the nodes
NumericMatrix get_root_node_coords(std::vector<neurostr::Node> root_nodes)
{

  auto n = std::distance(root_nodes.begin(), root_nodes.end());

  // Create nx3 matrix
  NumericMatrix m (n,3);

  // Fill it
  int i = 0;
  for (auto it = root_nodes.begin(); it != root_nodes.end(); ++it, ++i) {
    m(i, 0) = it->x();
    m(i, 1) = it->y();
    m(i, 2) = it->z();
  }

  return(m);
}

//Get the coordinates of the nodes
NumericMatrix get_contour_coords(std::vector<point_type> contour_points)
{

  auto n = std::distance(contour_points.begin(), contour_points.end());
  // Create nx3 matrix
  NumericMatrix m (n,3);

  // Fill it
  int i = 0;
  for (auto it = contour_points.begin(); it != contour_points.end(); ++it, ++i) {
    m(i, 0) = it->get<0>();
    m(i, 1) = it->get<1>();
    m(i, 2) = it->get<2>();
  }

  return(m);
}


// [[Rcpp::export]]
List c_get_root_neurite_nodes(std::string json_info) {
  std::ifstream aux;
  rapidjson::Document doc;
  doc.Parse(json_info.c_str());
  neurostr::io::JSONParser* p = new neurostr::io::JSONParser(aux);

  std::vector<neurostr::Node> root_nodes;//Terminal nodes of neurite before cutting
  std::vector<neurostr::Node> soma_nodes;
  std::vector<point_type> contour_points;

  auto r = p->read_string("reconstruction",doc);

  delete p;
  bool first = true;
  Rcout<<r->n_contours();
  for(auto contour = r->contour_begin(); contour!= r->contour_end(); ++contour)
  {
    neurostr::Contour c = *contour;
    for(auto it = c.begin(); it!= c.end(); ++it)
    {
      contour_points.emplace_back(*it);
    }
  }
  for(auto n_it = r->begin(); n_it != r->end(); ++n_it){
    neurostr::Neuron& n = *n_it;


    if(n.has_soma()){
      for(auto it = n.begin_soma(); it != n.end_soma(); ++it){
        soma_nodes.emplace_back(*it);
      }
    }

    for(auto it = n.begin_neurite(); it != n.end_neurite(); ++it){
      if(it->has_root())
        root_nodes.emplace_back(it->root());
      else
        root_nodes.emplace_back(* (it->begin_node()));
    }
  } // End neuron for
  NumericMatrix root_node_coords=get_root_node_coords(root_nodes);
  NumericMatrix soma_nodes_coords=get_root_node_coords(soma_nodes);
  NumericMatrix contour_nodes_coords=get_contour_coords(contour_points);


  return(List::create(Named("root_node")=root_node_coords,Named("soma_nodes")=soma_nodes_coords,Named("contour_nodes")=contour_nodes_coords));
}
