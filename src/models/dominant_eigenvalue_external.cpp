#include <stan/math.hpp>
#include <ostream>
#include <boost/math/tools/promotion.hpp>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace cambridgeModel_model_namespace {

template <typename T0__>
stan::promote_args_t<stan::value_type_t<T0__>>
dominant_eigenvalue_external(const T0__& A_arg__, std::ostream* pstream__){
  //Return type
  typedef typename stan::promote_args_t<stan::value_type_t<T0__>> scalar_t;
  //inputtype
  typedef T0__ matrix_t;

  //Eigensolver
  Eigen::EigenSolver<matrix_t> es(A_arg__);

  //Eigenvalues
  Eigen::Matrix<std::complex<scalar_t>, Eigen::Dynamic, 1> lambdas = es.eigenvalues();

  //Just take the real parts
  std::vector<scalar_t> lambdas_real(lambdas.rows());
  for (int i = 0; i < lambdas.rows(); i++){
    lambdas_real.at(i) = lambdas(i,0).real();
  }


  scalar_t returnValue;
  //Maximum eigenvalue (basic reproduction number is positive so
  //so we don't have to worry about sign)
  returnValue = stan::math::max(lambdas_real);

  return returnValue;

}

} //cambridgemodel_model_namespace
