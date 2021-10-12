
//
// Created by Alec Glisman on 07/30/21
// @SOURCE: https://stackoverflow.com/a/48963679
//

#ifndef HELPER_EIGEN_TENSOR_CONVERSION_H
#define HELPER_EIGEN_TENSOR_CONVERSION_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_MKL_ALL
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>

template <typename T> using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/**
 * @brief Convert between Eigen::Tensor and Eigen::Matrix classes
 *
 * @note example useage:
            Eigen::Tensor<double, 4> my_rank4(2, 2, 2, 2);
            my_rank4.setRandom();

            Eigen::MatrixXd          mymatrix = MatrixCast(my_rank4, 4, 4);    // Cast Eigen::Tensor --> Eigen::Matrix
            Eigen::Tensor<double, 3> my_rank3 = TensorCast(mymatrix, 2, 2, 4); // Cast Eigen::Matrix --> Eigen::Tensor
 *
 */

/**
 *
 * @brief Convert Eigen::Tensor  --> Eigen::Matrix
 *
 */

// Evaluates tensor expressions if needed
template <typename T, typename Device = Eigen::DefaultDevice>
auto
asEval(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>&
                     expr,             // An Eigen::TensorBase object (Tensor, TensorMap, TensorExpr... )
       const Device& device = Device() // Override to evaluate on another device, e.g. thread pool or gpu.
)
{
    using Evaluator                         = Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const T>, Device>;
    Eigen::TensorForcedEvalOp<const T> eval = expr.eval();
    Evaluator                          tensor(eval, device);
    tensor.evalSubExprsIfNeeded(nullptr);
    return tensor;
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Matrix with shape rows/cols
template <typename T, typename sizeType, typename Device = Eigen::DefaultDevice>
auto
MatrixCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const sizeType rows, const sizeType cols,
           const Device& device = Device())
{
    auto tensor  = asEval(expr, device);
    using Scalar = typename Eigen::internal::remove_const<typename decltype(tensor)::Scalar>::type;
    return static_cast<MatrixType<Scalar>>(Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols));
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Vector with the same size
template <typename T, typename Device = Eigen::DefaultDevice>
auto
VectorCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const Device& device = Device())
{
    auto tensor  = asEval(expr, device);
    auto size    = Eigen::internal::array_prod(tensor.dimensions());
    using Scalar = typename Eigen::internal::remove_const<typename decltype(tensor)::Scalar>::type;
    return static_cast<VectorType<Scalar>>(Eigen::Map<const VectorType<Scalar>>(tensor.data(), size));
}

// View an existing Eigen::Tensor as an Eigen::Map<Eigen::Matrix>
template <typename Scalar, auto rank, typename sizeType>
auto
MatrixMap(const Eigen::Tensor<Scalar, rank>& tensor, const sizeType rows, const sizeType cols)
{
    return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}

// View an existing Eigen::Tensor of rank 2 as an Eigen::Map<Eigen::Matrix>
// Rows/Cols are determined from the matrix
template <typename Scalar>
auto
MatrixMap(const Eigen::Tensor<Scalar, 2>& tensor)
{
    return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), tensor.dimension(0), tensor.dimension(1));
}

// View an existing Eigen::Tensor of rank 1 as an Eigen::Map<Eigen::Vector>
// Rows is the same as the size of the tensor.
template <typename Scalar, auto rank>
auto
VectorMap(const Eigen::Tensor<Scalar, rank>& tensor)
{
    return Eigen::Map<const VectorType<Scalar>>(tensor.data(), tensor.size());
}

/**
 *
 * @brief Convert Eigen::Matrix --> Eigen::Tensor
 *
 */

// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions specified in std::array
template <typename Derived, typename T, auto rank>
Eigen::Tensor<typename Derived::Scalar, rank>
TensorCast(const Eigen::EigenBase<Derived>& matrix, const std::array<T, rank>& dims)
{
    return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(),
                                                                                       dims);
}

// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions specified in Eigen::DSizes
template <typename Derived, typename T, auto rank>
Eigen::Tensor<typename Derived::Scalar, rank>
TensorCast(const Eigen::EigenBase<Derived>& matrix, const Eigen::DSizes<T, rank>& dims)
{
    return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().eval().data(),
                                                                                       dims);
}

// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions as variadic arguments
template <typename Derived, typename... Dims>
auto
TensorCast(const Eigen::EigenBase<Derived>& matrix, const Dims... dims)
{
    static_assert(sizeof...(Dims) > 0, "TensorCast: sizeof... (Dims) must be larger than 0");
    return TensorCast(matrix, std::array<Eigen::Index, sizeof...(Dims)>{dims...});
}

// Converts an Eigen::Matrix (or expression) to Eigen::Tensor
// with dimensions directly as arguments in a variadic template
template <typename Derived>
auto
TensorCast(const Eigen::EigenBase<Derived>& matrix)
{
    if constexpr (Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1)
    {
        return TensorCast(matrix, matrix.size());
    }
    else
    {
        return TensorCast(matrix, matrix.rows(), matrix.cols());
    }
}

// View an existing Eigen::Matrix as Eigen::TensorMap
// with dimensions specified in std::array
template <typename Derived, auto rank>
auto
TensorMap(const Eigen::PlainObjectBase<Derived>& matrix, const std::array<long, rank>& dims)
{
    return Eigen::TensorMap<const Eigen::Tensor<const typename Derived::Scalar, rank>>(matrix.derived().data(), dims);
}

// View an existing Eigen::Matrix as Eigen::TensorMap
// with dimensions as variadic arguments
template <typename Derived, typename... Dims>
auto
TensorMap(const Eigen::PlainObjectBase<Derived>& matrix, const Dims... dims)
{
    return TensorMap(matrix, std::array<long, static_cast<int>(sizeof...(Dims))>{dims...});
}

// View an existing Eigen::Matrix as Eigen::TensorMap
// with dimensions determined automatically from the given matrix
template <typename Derived>
auto
TensorMap(const Eigen::PlainObjectBase<Derived>& matrix)
{
    if constexpr (Derived::ColsAtCompileTime == 1 or Derived::RowsAtCompileTime == 1)
    {
        return TensorMap(matrix, matrix.size());
    }
    else
    {
        return TensorMap(matrix, matrix.rows(), matrix.cols());
    }
}

#endif