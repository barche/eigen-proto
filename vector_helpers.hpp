#include <vector>
#include <iostream>


namespace std
{

}
namespace eigen_proto
{

typedef std::vector<double> double_vec_t;
  
/// Initialize a vector with a sequence
inline void make_sequence(double_vec_t& v, int size)
{
  v.reserve(size);
  for(int i = 0; i != size; ++i)
    v.push_back(static_cast<double>(i));
}
  
}
