#include <iostream>
#include <boost/multi_array.hpp>

#include <boost/proto/proto.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/for_each.hpp>

#include <Eigen/Dense>


using namespace boost;

namespace fem
{

/// Unstructured mesh data
struct mesh_data
{
  mesh_data(const int nb_nodes, const int nb_elems, const int nb_nodes_per_elem, const int dimension) :
    coordinates(boost::extents[nb_nodes][dimension]),
    connectivity(boost::extents[nb_elems][nb_nodes_per_elem])
  {
  }

  typedef boost::multi_array<double, 2> coordinates_t;
  typedef boost::multi_array<int, 2> connectivity_t;
  coordinates_t coordinates;
  connectivity_t connectivity;
};

/// 1D Line shape function
struct line1d
{
  static const int nb_nodes = 2;
  static const int dimension = 1;

  typedef Eigen::Matrix<double, 1, dimension> coord_t;
  typedef Eigen::Matrix<double, 1, nb_nodes> shape_func_t;

  static shape_func_t shape_function(const coord_t& c)
  {
    const double xi = c[0];
    shape_func_t result;
    result[0] = 0.5*(1.-xi);
    result[1] = 0.5*(1.+xi);
    return result;
  }

  static const coord_t& centroid()
  {
    static const coord_t c = coord_t::Constant(0.);
    return c;
  }
};

/// 2D Triangle shape function
struct triag2d
{
  static const int nb_nodes = 3;
  static const int dimension = 2;

  typedef Eigen::Matrix<double, 1, dimension> coord_t;
  typedef Eigen::Matrix<double, 1, nb_nodes> shape_func_t;

  static shape_func_t shape_function(const coord_t& c)
  {
    const double xi = c[0]; const double eta = c[1];
    shape_func_t result;
    result[0] = xi;
    result[1] = eta;
    result[2] = 1. - xi - eta;
    return result;
  }

  static const coord_t& centroid()
  {
    static const coord_t c = coord_t::Constant(1./3.);
    return c;
  }
};

}

namespace dsl
{

/// Helper data for proto
template<typename ElementT>
struct dsl_data
{
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  typedef ElementT element_t;
  typedef Eigen::Matrix<double, element_t::nb_nodes, element_t::dimension> coord_mat_t;

  dsl_data(const fem::mesh_data& d) : mesh_data(d)
  {
  }

  void set_element(const int e)
  {
    for(int i = 0; i != element_t::nb_nodes; ++i)
      for(int j = 0; j != element_t::dimension; ++j)
        coord_mat(i,j) = mesh_data.coordinates[mesh_data.connectivity[e][i]][j];
  }

  const fem::mesh_data& mesh_data;
  coord_mat_t coord_mat;
  typename element_t::shape_func_t shape_func;
};

struct element_coords_tag {};

struct eval_element_coord : proto::callable
{
  template<typename Signature>
  struct result;

  template<class ThisT, typename DataT>
  struct result<ThisT(DataT)>
  {
    typedef const typename boost::remove_reference<DataT>::type::coord_mat_t& type;
  };

  template<typename DataT>
  const typename DataT::coord_mat_t& operator()(const DataT& data) const
  {
    return data.coord_mat;
  }
};

/// Return the element coordinates matrix
proto::terminal< element_coords_tag >::type const element_coords = {};

struct centroid_tag {};

struct eval_centroid : proto::callable
{
  template<typename Signature>
  struct result;

  template<class ThisT, typename DataT>
  struct result<ThisT(DataT)>
  {
    typedef const typename boost::remove_reference<DataT>::type::element_t::coord_t& type;
  };

  template<typename DataT>
  const typename DataT::element_t::coord_t& operator()(const DataT&) const
  {
    return DataT::element_t::centroid();
  }
};

/// Return the element coordinates matrix
proto::terminal< centroid_tag >::type const centroid = {};

struct shape_func_tag {};

struct eval_shape_func : proto::callable
{
  template<typename Signature>
  struct result;

  template<class ThisT, typename CoordT, typename DataT>
  struct result<ThisT(CoordT, DataT)>
  {
    typedef const typename boost::remove_reference<DataT>::type::element_t::shape_func_t& type;
  };

  template<typename CoordT, typename DataT>
  const typename DataT::element_t::shape_func_t& operator()(const CoordT& coord, DataT& data) const
  {
    data.shape_func = DataT::element_t::shape_function(coord);
    return data.shape_func;
  }
};

/// Return the element coordinates matrix
proto::terminal< shape_func_tag >::type const N = {};

/// Helper for output
proto::terminal< std::ostream & >::type cout_ = { std::cout };

/// Proto grammar
struct fem_grammar :
  proto::or_
  <
    proto::when<proto::terminal<element_coords_tag>, eval_element_coord(proto::_data)>,
    proto::when<proto::terminal<centroid_tag>, eval_centroid(proto::_data)>,
    proto::when
    <
      proto::function<proto::terminal<shape_func_tag>, proto::_>,
      eval_shape_func(fem_grammar(proto::_child1), proto::_data)
    >,
    proto::_default<fem_grammar>
  >
{
};

/// MPL functor to loop over elements
template<typename ExprT>
struct element_looper
{
  element_looper(const fem::mesh_data& m, const ExprT& e) : mesh(m), expr(e)
  {
  }

  template < typename ElemT >
  void operator()(ElemT) const
  {
    // Bail out if the shape function doesn't match
    if(ElemT::dimension != mesh.coordinates.shape()[1] || ElemT::nb_nodes != mesh.connectivity.shape()[1])
      return;

    // Construct helper data
    dsl_data<ElemT> data(mesh);

    // Evaluate for all elements
    const int nb_elems = mesh.connectivity.size();
    for(int i = 0; i != nb_elems; ++i)
    {
      data.set_element(i);
      fem_grammar()(expr, 0, data);
    }
  }

  const fem::mesh_data& mesh;
  const ExprT& expr;
};

/// Execute the given expression for every element in the mesh
template<typename ExprT>
void for_each_element(const fem::mesh_data& mesh, const ExprT& expr)
{
  // Allowed element types
  typedef mpl::vector2<fem::line1d, fem::triag2d> element_types;
  mpl::for_each<element_types>(element_looper<ExprT>(mesh, expr));
}

}

/// Demonstrate some basic proto expression properties
int main(void)
{
  using namespace dsl;
  using namespace fem;

  mesh_data line_mesh(4, 3, 2, 1);
  line_mesh.coordinates[0][0] = 0.;
  line_mesh.coordinates[1][0] = 1.;
  line_mesh.coordinates[2][0] = 2.;
  line_mesh.coordinates[3][0] = 3.;

  line_mesh.connectivity[0][0] = 0;
  line_mesh.connectivity[0][1] = 1;
  line_mesh.connectivity[1][0] = 1;
  line_mesh.connectivity[1][1] = 2;
  line_mesh.connectivity[2][0] = 2;
  line_mesh.connectivity[2][1] = 3;

  std::cout << "line mesh centroids:" << std::endl;
  for_each_element(line_mesh, cout_ << N(centroid)*element_coords << "\n");

  mesh_data triag_mesh(4, 2, 3, 2);
  triag_mesh.coordinates[0][0] = 0.;
  triag_mesh.coordinates[0][1] = 0.;
  triag_mesh.coordinates[1][0] = 1.;
  triag_mesh.coordinates[1][1] = 0.;
  triag_mesh.coordinates[2][0] = 1.;
  triag_mesh.coordinates[2][1] = 1.;
  triag_mesh.coordinates[3][0] = 0.;
  triag_mesh.coordinates[3][1] = 1.;

  triag_mesh.connectivity[0][0] = 0;
  triag_mesh.connectivity[0][1] = 1;
  triag_mesh.connectivity[0][2] = 2;
  triag_mesh.connectivity[1][0] = 0;
  triag_mesh.connectivity[1][1] = 2;
  triag_mesh.connectivity[1][2] = 3;

  std::cout << "triag mesh centroids:" << std::endl;
  for_each_element(triag_mesh, cout_ << N(centroid)*element_coords << "\n");
}
