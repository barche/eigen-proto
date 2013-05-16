#include <boost/assert.hpp>
#include <boost/proto/proto.hpp>
#include <Eigen/Dense>

using namespace boost;

namespace eigen_proto
{

/// Helper to get the actual matrix type that comes out of a product
template<typename LeftT, typename RightT>
struct product_value_type
{
  typedef typename remove_reference<LeftT>::type left_unref;
  typedef typename remove_reference<RightT>::type right_unref;
  typedef typename Eigen::ProductReturnType<left_unref, right_unref>::Type product_type;
  typedef typename Eigen::MatrixBase<product_type>::PlainObject type;
};

struct do_eigen_multiply : proto::callable
{
  template<typename Signature>
  struct result;

  template<class ThisT, class ExprT, class LeftT, class RightT>
  struct result<ThisT(ExprT, LeftT, RightT)>
  {
    typedef typename product_value_type<LeftT, RightT>::type& type;
  };

  template<typename ExprT, typename LeftT, typename RightT>
  typename product_value_type<LeftT, RightT>::type& operator()(ExprT& expr, const LeftT& left, const RightT& right) const
  {
    expr.value = left*right;
    return expr.value;
  }
};

/// Simple calculator grammar and transform for any type supporting arithmetic operators
struct calculator_transform :
  proto::or_
  <
    proto::when< proto::terminal<proto::_>, proto::_value >, // Replace terminals with their value
    proto::when
    <
       proto::multiplies< calculator_transform, calculator_transform >,
       do_eigen_multiply(proto::_, calculator_transform(proto::_left), calculator_transform(proto::_right))
    >,
    proto::when<proto::or_ // When we have +, - or / ...
    <
      proto::plus<calculator_transform, calculator_transform>,
      proto::minus<calculator_transform, calculator_transform>,
      proto::divides<calculator_transform, calculator_transform>
    >,
    proto::_default<calculator_transform> > // ... Do what C++ would do
  >
{
};

/// Wraps a given expression, so the value that it represents can be stored inside the expression itself
template<typename ExprT, typename ValueT>
struct stored_result_expression :
  proto::extends< ExprT, stored_result_expression<ExprT, ValueT> >
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  typedef proto::extends< ExprT, stored_result_expression<ExprT, ValueT> > base_type;

  explicit stored_result_expression(ExprT const &expr = ExprT())
    : base_type(expr)
  {
  }

  /// Temporary storage for the result of the expression
  mutable ValueT value;
};

struct do_wrap_expression : proto::transform< do_wrap_expression >
{
  template<typename ExprT, typename StateT, typename DataT>
  struct impl : proto::transform_impl<ExprT, StateT, DataT>
  {
    typedef typename result_of<calculator_transform(ExprT, StateT, DataT)>::type result_ref_type;
    typedef typename remove_reference<result_ref_type>::type value_type;
    typedef typename remove_const<typename remove_reference<ExprT>::type>::type expr_val_type;
    typedef stored_result_expression<expr_val_type, value_type> result_type;

    result_type operator()(typename impl::expr_param expr, typename impl::state_param state, typename impl::data_param data)
    {
      return result_type(expr);
    }
  };
};

/// Wrap multiplies expressions so they can store a temporary result
struct wrap_expression :
  proto::or_
  <
    proto::terminal<proto::_>,
    proto::when
    <
      proto::multiplies<proto::_, proto::_>,
      do_wrap_expression(proto::functional::make_multiplies
      (
        wrap_expression(proto::_left), wrap_expression(proto::_right)
      ))
    >,
    proto::nary_expr< proto::_, proto::vararg<wrap_expression> >
  >
{
};

// Print the result type of an expression evaluated with calculator_transform
template<typename ExprT>
void print_result_type(const ExprT& expr)
{
  result_of<calculator_transform(ExprT)>::type::print_error();
}

}

/// Demonstrate some basic proto expression properties
int main(void)
{
  using namespace eigen_proto;
  using proto::lit;

  wrap_expression wrap;
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
  const CT result1 = eval(wrap(b*a));
  std::cout << "First test expected:\n" << expected1 << "\nobtained:\n" << result1 << "\n" << std::endl;
  BOOST_ASSERT(result1 == expected1);

  const CT expected2 = (b_mat*a_mat)*c_mat;

  proto::display_expr((b*a)*c);

  std::cout << "Second test expected:\n" << expected2 << std::endl;
  const CT result2 = eval(wrap((b*a)*c));
  std::cout << "obtained:\n" << result2 << "\n" << std::endl;
  BOOST_ASSERT(result2 == expected2);

  const CT expected3 = (b_mat*a_mat)*c_mat*(b_mat*a_mat)+c_mat*b_mat*a_mat-(c_mat+c_mat)*b_mat*a_mat;
  std::cout << "Third test expected:\n" << expected3 << std::endl;
  const CT result3 = eval(wrap((b*a)*c*(b*a)+c*b*a-(c+c)*b*a));
  std::cout << "obtained:\n" << result3 << "\n" << std::endl;
  BOOST_ASSERT(result3 == expected3);

  return 0;
}
