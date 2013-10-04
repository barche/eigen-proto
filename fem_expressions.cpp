#include <iostream>
#include <boost/proto/proto.hpp>

using namespace boost;

// Different types to distinguish terminals at compile time
struct element_coords_tag {};
struct centroid_tag {};
struct shape_func_tag {};

// Some statically created terminals
proto::terminal< centroid_tag >::type const centroid = {};
proto::terminal< element_coords_tag >::type const element_coords = {};
proto::terminal< shape_func_tag >::type const N = {};
proto::terminal< std::ostream & >::type cout_ = { std::cout };

int main(void)
{
  // The terminals can be combined into any valid C++ expression
  (centroid + element_coords * cout_) / N;
  cout_ << N(centroid)*element_coords << "\n";
  proto::display_expr(cout_ << N(centroid)*element_coords << "\n");
}
