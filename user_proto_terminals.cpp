#include <boost/proto/proto.hpp>

using namespace boost;

namespace eigen_proto
{

template<typename T>
struct user_op {};

struct my_callable {};

template<typename GrammarT, typename ExprT>
void assert_match(const ExprT&)
{
  static_assert(proto::matches<ExprT, GrammarT>::value, "Expression doesn't match grammar");
}

}

int main(void)
{
  using namespace boost;
  using namespace eigen_proto;

  // Create a terminal
  boost::proto::terminal< user_op<my_callable> >::type const my_op = {};

  // This is a valid expression
  proto::display_expr(my_op(1,2,3));

  // This grammar matches it:
  assert_match< proto::function< proto::terminal< user_op<proto::_> >, proto::vararg<proto::_> > >(my_op(1,2,3));

  return 0;
}
