#include <boost/proto/proto.hpp>

using namespace boost;

// Different types to distinguish terminals at compile time
struct shape_func_tag {};
struct element_quadrature_tag {};
struct transpose_tag {};

// Some statically created terminals
proto::terminal< shape_func_tag >::type const N = {};
proto::terminal< element_quadrature_tag >::type const element_quadrature = {};
proto::terminal< transpose_tag >::type const transpose = {};

int main(void)
{
  double* M;
  // This is a proto expression:
  element_quadrature(M += transpose(N)*N);
  proto::display_expr(element_quadrature(M += transpose(N)*N));
  return 0;
}
