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
  // Construct using a number of nodes and element properties
  mesh_data(const int nb_nodes, const int nb_elems, const int nb_nodes_per_elem, const int dimension) :
    coordinates(boost::extents[nb_nodes][dimension]),
    connectivity(boost::extents[nb_elems][nb_nodes_per_elem])
  {
  }

  // Global coordinates array
  boost::multi_array<double, 2> coordinates;

  // Global connectivity array
  boost::multi_array<int, 2> connectivity;
};

/// 1D Line shape function
struct line1d
{
  static const int nb_nodes = 2; // Number of nodes
  static const int dimension = 1;

  // Type of the mapped coordinates
  typedef Eigen::Matrix<double, 1, dimension> coord_t;
  // Type of the shape function vector
  typedef Eigen::Matrix<double, 1, nb_nodes> shape_func_t;
  // Type of the coordinates matrix
  typedef Eigen::Matrix<double, nb_nodes, dimension> coord_mat_t;

  // Compute the shape function vector at mapped coordinate c
  static shape_func_t shape_function(const coord_t& c)
  {
    const double xi = c[0];
    shape_func_t result;
    result[0] = 0.5*(1.-xi);
    result[1] = 0.5*(1.+xi);
    return result;
  }
  
  // Compute the jacobian determinant
  static double jacobian_determinant(const coord_t& mapped_coord, const coord_mat_t& node_coords)
  {
    return 0.5*(node_coords[1] - node_coords[0]);
  }
  
  static const int nb_gauss_points = 2;
  
  // Type of the matrix with the Gauss points
  typedef Eigen::Matrix<double, nb_gauss_points, 1> gauss_points_t;
  
  // The Gauss points for the current shape function
  static const gauss_points_t gauss_points()
  {
    gauss_points_t points;
    points[0] = 0.5773502691896257645091488;
    points[1] = -points[0];
    return points;
  }
  
  // Type for the weights
  typedef Eigen::Matrix<double, nb_gauss_points, 1> gauss_weights_t;
  
  // The Gauss weights
  static const gauss_weights_t gauss_weights()
  {
    gauss_weights_t weights;
    weights.setConstant(1.);
    return weights;
  }
};

/// 2D Triangle shape function
struct triag2d
{
  static const int nb_nodes = 3;
  static const int dimension = 2;

  typedef Eigen::Matrix<double, 1, dimension> coord_t;
  typedef Eigen::Matrix<double, 1, nb_nodes> shape_func_t;
  typedef Eigen::Matrix<double, nb_nodes, dimension> coord_mat_t;

  static shape_func_t shape_function(const coord_t& c)
  {
    const double xi = c[0]; const double eta = c[1];
    shape_func_t result;
    result[0] = xi;
    result[1] = eta;
    result[2] = 1. - xi - eta;
    return result;
  }

  // Compute the jacobian determinant
  static double jacobian_determinant(const coord_t& mapped_coord, const coord_mat_t& node_coords)
  {
    const Eigen::Vector3d x = node_coords.col(0);
    const Eigen::Vector3d y = node_coords.col(1);
    return (x[1] - x[0]) * (y[2] - y[0]) - (x[2] - x[0]) * (y[1] - y[0]);
  }

  static const int nb_gauss_points = 3;

  // Type of the matrix with the Gauss points
  typedef Eigen::Matrix<double, nb_gauss_points, 2> gauss_points_t;

  // The Gauss points for the current shape function
  static const gauss_points_t gauss_points()
  {
    gauss_points_t points;
    points(0,0) = 0.5; points(0,1) = 0.;
    points(1,0) = 0.5; points(1,1) = 0.5;
    points(2,0) = 0.; points(2,1) = 0.5;
    return points;
  }

  // Type for the weights
  typedef Eigen::Matrix<double, nb_gauss_points, 1> gauss_weights_t;

  // The Gauss weights
  static const gauss_weights_t gauss_weights()
  {
    gauss_weights_t weights;
    weights.setConstant(1./6.);
    return weights;
  }
};

}

namespace dsl
{

/// Helper data for proto
template<typename ElementT>
struct dsl_data
{
  // Required by Eigen to store fixed-size data
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // The concrete element type
  typedef ElementT element_t;
  
  // Construct using the mesh
  dsl_data(const fem::mesh_data& d) : mesh_data(d)
  {
  }

  // Set the current element
  void set_element(const int e)
  {
    for(int i = 0; i != element_t::nb_nodes; ++i)
    {
      const int node_idx = mesh_data.connectivity[e][i];
      node_indices[i] = node_idx;
      for(int j = 0; j != element_t::dimension; ++j)
        coord_mat(i,j) = mesh_data.coordinates[node_idx][j];
    }
  }

  // Reference to the mesh
  const fem::mesh_data& mesh_data;
  // Storage for the coordinates of the current element nodes
  typename ElementT::coord_mat_t coord_mat;
  // Global indices of the nodes of the current element
  int node_indices[element_t::nb_nodes];
  // Value of the last shape function computation
  typename element_t::shape_func_t shape_func;
  // Value of the last Jacobian determinant computation
  double det_jacobian;
};

struct shape_func_tag {};

struct eval_shape_func : proto::callable
{
  // C++ result_of declaration
  template<typename Signature>
  struct result;

  // C++ result_of implementation
  template<class ThisT, typename DataT>
  struct result<ThisT(DataT)>
  {
    typedef const typename
      boost::remove_reference<DataT>::type::element_t::shape_func_t& type;
  };

  template<typename DataT>
  const typename DataT::element_t::shape_func_t& operator()(DataT& data) const
  {
    // Return a reference to the result, stored in data
    return data.shape_func;
  }
};

// Evaluation of element quadrature
struct element_quadrature_tag {};

// Forward declaration
struct eval_element_quadrature;

/// Perform integration of an expression of the form A += B
proto::terminal< element_quadrature_tag >::type const element_quadrature = {};

struct transpose_tag {};

// Evaluation of matrix transpose
struct eval_transpose : proto::callable
{
  // C++ result_of declaration
  template<typename Signature>
  struct result;

  // C++ result_of implementation
  template<class ThisT, typename MatT>
  struct result<ThisT(MatT)>
  {
		typedef Eigen::Transpose<const typename boost::remove_reference<MatT>::type> type;
  };

  template<typename MatT>
  Eigen::Transpose<MatT> operator()(MatT& mat) const
  {
    return mat.transpose();
  }
};

proto::terminal<transpose_tag>::type const transpose = {};

/// Return the shape function value at the given mapped coordinate
proto::terminal< shape_func_tag >::type const N = {};

/// Helper for output
proto::terminal< std::ostream & >::type cout_ = { std::cout };

// Grammar to evaluate the expressions
struct fem_grammar :
  // Match the following rules in order:
  proto::or_
	<
    // Evaluate shape functions using the eval_shape_func transform
    proto::when
    <
      proto::terminal<shape_func_tag>,
      eval_shape_func(proto::_data)
    >,
    // Evaluate transpose expressions using eval_transpose
    proto::when
    <
			proto::function<proto::terminal<transpose_tag>, fem_grammar >,
      eval_transpose(fem_grammar(proto::_child1))
    >,
		// Evaluate element_quadrature using eval_element_quadrature
		proto::when
		<
			proto::function< proto::terminal<element_quadrature_tag>, proto::plus_assign<proto::terminal<Eigen::MatrixXd>, fem_grammar> >,
			eval_element_quadrature(proto::_value(proto::_left(proto::_child1)), proto::_right(proto::_child1), proto::_data)
		>,
    // On any other expression: perform the default C++ action
    proto::_default<fem_grammar>
  >
{
};

/// Integral evaluation. Needs to be after fem_grammar because it uses it internally.
struct eval_element_quadrature : proto::callable
{
  typedef void result_type;

  template<typename ExprT, typename DataT>
  void operator()(Eigen::MatrixXd& mat, const ExprT& expr, DataT& data) const
  {
    // Temporary storage for the result from the RHS
    Eigen::Matrix<double, DataT::element_t::nb_nodes, DataT::element_t::nb_nodes> rhs_result;
    rhs_result.setZero();

    // Loop over the gauss points
    const typename DataT::element_t::gauss_points_t gauss_points = DataT::element_t::gauss_points();
    const typename DataT::element_t::gauss_weights_t gauss_weights = DataT::element_t::gauss_weights();
    for(int i = 0; i != DataT::element_t::nb_gauss_points; ++i)
    {
      data.shape_func = DataT::element_t::shape_function(gauss_points.row(i));
      data.det_jacobian = DataT::element_t::jacobian_determinant(gauss_points.row(i), data.coord_mat);
      rhs_result += gauss_weights[i] * data.det_jacobian * fem_grammar()(expr, 0, data);
    }

    // Place the result in the global matrix
    for(int i = 0; i != DataT::element_t::nb_nodes; ++i)
      for(int j = 0; j != DataT::element_t::nb_nodes; ++j)
        mat(data.node_indices[i], data.node_indices[j]) += rhs_result(i,j);
  }
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

  Eigen::MatrixXd M(4,4); M.setZero();
  for_each_element(line_mesh, element_quadrature(M += transpose(N)*N));

  std::cout << "Line mass matrix:\n" << M << std::endl;

   mesh_data triag_mesh(4, 2, 3, 2);
   triag_mesh.coordinates[0][0] = 0.;
   triag_mesh.coordinates[0][1] = 0.;
   triag_mesh.coordinates[1][0] = 1.;
   triag_mesh.coordinates[1][1] = 0.;
   triag_mesh.coordinates[2][0] = 0.;
   triag_mesh.coordinates[2][1] = 1.;
   triag_mesh.coordinates[3][0] = 1.;
   triag_mesh.coordinates[3][1] = 1.;

   triag_mesh.connectivity[0][0] = 0;
   triag_mesh.connectivity[0][1] = 3;
   triag_mesh.connectivity[0][2] = 2;
   triag_mesh.connectivity[1][0] = 0;
   triag_mesh.connectivity[1][1] = 1;
   triag_mesh.connectivity[1][2] = 3;

   Eigen::MatrixXd L(4,4); L.setZero();
   for_each_element(triag_mesh, element_quadrature(L += transpose(transpose(transpose(N)))*transpose(transpose(N))));

   std::cout << "Triangle mass matrix:\n" << L << std::endl;
}
