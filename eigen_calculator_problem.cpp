#include <boost/proto/proto.hpp>
#include <Eigen/Dense>

using namespace boost;

namespace eigen_proto
{

/// Simple calculator grammar and transform for any type supporting arithmetic operators
struct calculator_transform :
  proto::or_
  <
    proto::when< proto::terminal<proto::_>, proto::_value >, // Replace terminals with their value
    proto::when<proto::or_ // When we have +, - * or / ...
    <
      proto::plus<calculator_transform, calculator_transform>,
      proto::minus<calculator_transform, calculator_transform>,
      proto::multiplies<calculator_transform, calculator_transform>,
      proto::divides<calculator_transform, calculator_transform>
    >,
    proto::_default<calculator_transform> > // ... Do what C++ would do
  >
{
};

// Print the result type of an expression evaluated with calculator_transform
template<typename ExprT>
void print_result_type(const ExprT& expr)
{
  boost::result_of<calculator_transform(ExprT)>::type::print_error();
}

}

/// Demonstrate some basic proto expression properties
int main(void)
{
  using namespace eigen_proto;
  using proto::lit;

  calculator_transform eval;

  // Typedefs for matrix types
  typedef Eigen::Matrix<double, 1, 2> AT;
  typedef Eigen::Matrix<double, 2, 1> BT;
  typedef Eigen::Matrix<double, 2, 2> CT;

  // Construct the matrices
  AT a_mat; a_mat.setConstant(2);
  BT b_mat; b_mat.setConstant(2);
  CT c_mat; c_mat.setConstant(1);

  // Build proto terminals
  proto::literal<AT&> a(a_mat);
  proto::literal<BT&> b(b_mat);
  proto::literal<CT&> c(c_mat);

  const CT expected1 = b_mat*a_mat;
  const CT result1 = eval(b*a);
  std::cout << "First test expected:\n" << expected1 << "\nobtained:\n" << result1 << "\n" << std::endl;
  BOOST_ASSERT(result1 == expected1);

  const CT expected2 = (b_mat*a_mat)*c_mat;

  proto::display_expr((b*a)*c);

  std::cout << "Second test expected:\n" << expected2 << std::endl;
  const CT result2 = eval((b*a)*c);
  std::cout << "obtained:\n" << result2 << "\n" << std::endl;
  BOOST_ASSERT(result2 == expected2);

  // Print types for analysis
  //print_result_type(b*a);
  //print_result_type((b*a)*c);

  return 0;
}
