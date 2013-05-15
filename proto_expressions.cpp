#include <boost/proto/proto.hpp>

using namespace boost;

/// Demonstrate some basic proto expression properties
int main(void)
{
  // We will demonstrate the expressions using a simple integer value
  int i = 1;
  
  // We can easily make a terminal out of it, held by reference here
  proto::literal<int&> int_term(i);
  
  // or do this in-place
  proto::display_expr(proto::lit(i));
  
  // We can get the value and change it using proto::value
  proto::value(int_term) = 2;
  BOOST_ASSERT(i == 2);
  std::cout << "Value after proto::value(int_term) = 2: " << i << std::endl;
  
  // The terminal can also be constructed in-place using proto::lit()
  proto::value(proto::lit(i)) = 3;
  BOOST_ASSERT(i == 3);
  std::cout << "Value after proto::value(proto::lit(a)) = 3: " << i << std::endl;
  
  // The following compiles, just building a "multiplies" expression:
  int_term * 2;
  // Same for the function call:
  int_term(2);
  // (Obviously, this won't work on the raw int)
  //i(lit(2));
  
  // Let's print some expressions:
  std::cout << "A simple terminal:" << std::endl;
  proto::display_expr(int_term);
  std::cout << "A multiplication and sum..." << std::endl;
  proto::display_expr(int_term * 2 + 3);
  std::cout << "Assignment..." << std::endl;
  proto::display_expr(int_term = 2);
  std::cout << "Function calls, with the first child being our terminal..." << std::endl;
  proto::display_expr(int_term(1,2));
  // This is one too many, redefine BOOST_PROTO_MAX_ARITY to something higher than 10 for it to work
  //int_term(1,2,3,4,5,6,7,8,9,10);

  proto::display_expr(int_term << "shift" << "things");
  
  return 0;
}
