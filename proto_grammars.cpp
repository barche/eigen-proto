#include <boost/proto/proto.hpp>

using namespace boost;

namespace eigen_proto
{

/// Simple grammar that only allows addition of integers
struct add_ints_grammar :
  proto::plus< proto::terminal<int>, proto::terminal<int> >
{
};

/// Helper function to check if an expression matches the grammar
template<typename ExprT>
void check_expr(const ExprT&)
{
  BOOST_MPL_ASSERT((
    proto::matches<ExprT, add_ints_grammar>
  ));
}

}

/// Demonstrate some basic proto expression properties
int main(void)
{
  // We will demonstrate the expressions using a simple integer value
  int i = 1;
  
  // Adding two ints is OK
  eigen_proto::check_expr(proto::lit(i) + 2);

  // Anything else including an unsigned int isn't:
  //eigen_proto::check_expr(proto::lit(i) + 2u);

  return 0;
}
