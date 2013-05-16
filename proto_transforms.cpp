#include <boost/assert.hpp>
#include <boost/proto/proto.hpp>

using namespace boost;

namespace eigen_proto
{

/// Simple grammar that only allows addition of integers
struct add_ints_grammar :
  proto::plus< proto::terminal<int>, proto::terminal<int> >
{
};

/// We need a functor that we can call to do the actual evaluation
struct evaluate_plus : proto::callable
{
  typedef int result_type;

  int operator()(const int a, const int b) const
  {
    return a + b;
  }
};

/// Now we can create a transform that can evaluate our expressions:
struct add_ints_transform :
  proto::when<
    add_ints_grammar,
    evaluate_plus(proto::_value(proto::_left), proto::_value(proto::_right))
  >
{
};

/// We can also do this simpler, combining different types of arithmetic
/// and using the proto default transform
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

}

/// Demonstrate some basic proto expression properties
int main(void)
{
  using namespace eigen_proto;
  using proto::lit;

  // Create an integer
  int i = 1;
  
  // Evaluate the expression using our transform
  const int result = add_ints_transform()(lit(i) + 2);

  BOOST_ASSERT(result == 3);
  std::cout << "i + 2 = " << result << std::endl;

  // The calculator is much more powerful:
  int a = 2; double b = 0.5; unsigned int c = 3;

  const double result2 = calculator_transform()((lit(a) + lit(c)) * lit(b));
  BOOST_ASSERT(result2 == 2.5);
  std::cout << "lit(a) + lit(c)) * lit(b) = " << result2 << std::endl;

  return 0;
}
