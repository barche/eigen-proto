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

/// Functor to evaluate a product wrapped to hold its result
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

/// Placeholder for user-defined operations
template<typename CallableT>
struct user_op
{
  typedef CallableT callable_t;
};

/// Primitive transform to evaluate user operations
struct evaluate_user_op : proto::transform< evaluate_user_op >
{
  template<typename ExprT, typename StateT, typename DataT>
  struct impl : proto::transform_impl<ExprT, StateT, DataT>
  {
    // Calculate the type of the functor that the user supplied
    typedef typename result_of<proto::_child0(ExprT)>::type callable_term_t;
    typedef typename result_of<proto::_value(callable_term_t)>::type callable_ref_t;
    typedef typename remove_reference<callable_ref_t>::type::callable_t callable_t;

    // Helper struct to calculate result types
    template<int NbArgs, int dummy = 0>
    struct compute_result_type;

    // Specialization for 1 argument
    template<int dummy>
    struct compute_result_type<1, dummy>
    {
      typedef typename result_of<proto::_child1(ExprT)>::type child1_t;
      typedef typename result_of<callable_t(child1_t)>::type type;
    };

    // Specialization for 2 arguments
    template<int dummy>
    struct compute_result_type<2, dummy>
    {
      typedef typename result_of<proto::_child1(ExprT)>::type child1_t;
      typedef typename result_of<proto::_child2(ExprT)>::type child2_t;
      typedef typename result_of<callable_t(child1_t, child2_t)>::type type;
    };

    // Number of function arguments
    static const int nb_args = proto::arity_of<ExprT>::value-1;

    // Result type of our user-supplied callable
    typedef typename compute_result_type<nb_args>::type result_type;

    result_type operator()(typename impl::expr_param expr, typename impl::state_param state, typename impl::data_param data)
    {
      return dispatch(mpl::int_<nb_args>(), expr);
    }

    // Dispatch on the number of arguments
    result_type dispatch(mpl::int_<1>, typename impl::expr_param expr)
    {
      return callable_t()(proto::child_c<1>(expr));
    }

    result_type dispatch(mpl::int_<2>, typename impl::expr_param expr)
    {
      return callable_t()(proto::child_c<1>(expr), proto::child_c<1>(expr));
    }
  };
};

/// Simple calculator grammar and transform for any type supporting arithmetic operators
struct calculator_transform :
  proto::or_
  <
    proto::when
    <
      proto::function< proto::terminal< user_op<proto::_> >, proto::vararg<proto::_> >,
      evaluate_user_op(proto::function< proto::_, proto::vararg<calculator_transform> >)
    >,
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

/// Primitive transform wrapping expressions to hold their value
struct do_wrap_expression : proto::transform< do_wrap_expression >
{
  template<typename ExprT, typename StateT, typename DataT>
  struct impl : proto::transform_impl<ExprT, StateT, DataT>
  {
    typedef typename result_of<calculator_transform(ExprT, StateT, DataT)>::type result_ref_type;
    typedef typename remove_reference<result_ref_type>::type value_type;
    typedef typename remove_const<typename remove_reference<ExprT>::type>::type expr_val_type;
    typedef const stored_result_expression<expr_val_type, value_type> result_type;

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

struct do_transpose
{
  template<typename Signature>
  struct result;

  template<class ThisT, typename MatrixT>
  struct result<ThisT(MatrixT)>
  {
    typedef Eigen::Transpose<MatrixT> type;
  };

  template<typename MatrixT>
  Eigen::Transpose<MatrixT> operator()(MatrixT mat)
  {
    return mat.transpose();
  }

  template<int Rows, int Cols>
  Eigen::Transpose< Eigen::Matrix<double, Rows, Cols> > operator()(Eigen::Matrix<double, Rows, Cols>& mat)
  {
    return mat.transpose();
  }
};

// Create the transpose op
proto::terminal< user_op<do_transpose> >::type const transpose = {};

}

int main(void)
{
  using namespace boost;
  using namespace eigen_proto;

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

  // Transforms to evaluate our expressions
  wrap_expression wrap;
  calculator_transform eval;

  // This is a valid expression
  c_mat.col(0).setConstant(5.);
  BOOST_ASSERT(c_mat.transpose() == eval(wrap(transpose(c))));
  BOOST_ASSERT(c_mat == eval(wrap(transpose(transpose(c)))));
  BOOST_ASSERT(c_mat.transpose() == eval(wrap(transpose(transpose(transpose(c))))));

  std::cout << eval(wrap(transpose(transpose(transpose(c))))) << std::endl;

  const CT expected2 = (b_mat*a_mat)*c_mat;
  std::cout << "expected:\n" << expected2 << std::endl;
  const CT result2 = eval(wrap((transpose(a)*transpose(b))*c));
  std::cout << "obtained:\n" << result2 << "\n" << std::endl;
  BOOST_ASSERT(result2 == expected2);

  return 0;
}
