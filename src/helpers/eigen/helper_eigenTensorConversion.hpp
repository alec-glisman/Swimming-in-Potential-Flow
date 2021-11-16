
//
// Created by Alec Glisman on 07/30/21
// @SOURCE: https://stackoverflow.com/a/48963679
//

#ifndef HELPER_EIGEN_TENSOR_CONVERSION_H
#define HELPER_EIGEN_TENSOR_CONVERSION_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

// Intel MKL
#if __has_include("mkl.h")
#define EIGEN_USE_MKL_ALL
#else
#pragma message(" !! COMPILING WITHOUT INTEL MKL OPTIMIZATIONS !! ")
#endif
// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
// eigen3 conversion between Eigen::Tensor (unsupported) and Eigen::Matrix
#include <helper_eigenTensorConversion.hpp>

template <typename T> using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template <typename T> using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;

/**
 * @brief Convert between Eigen::Tensor and Eigen::Matrix classes
 *
 * @see https://stackoverflow.com/questions/48795789/eigen-unsupported-tensor-to-eigen-matrix
 *
 * @note example useage:
 *
 *     Eigen::Tensor<double, 4> my_rank4(2, 2, 2, 2);
 *     my_rank4.setRandom();
 *
 *     Eigen::MatrixXd          mymatrix = MatrixCast(my_rank4, 4, 4);    // Cast Eigen::Tensor --> Eigen::Matrix
 *     Eigen::Tensor<double, 3> my_rank3 = TensorCast(mymatrix, 2, 2, 4); // Cast Eigen::Matrix --> Eigen::Tensor
 *
 */

/**
 *
 * @brief Convert Eigen::Tensor  --> Eigen::Matrix
 *
 */

template <typename Derived, typename Device = Eigen::DefaultDevice> class selfCleaningEvaluator
{
  private:
    using Evaluator = Eigen::TensorEvaluator<const Eigen::TensorForcedEvalOp<const Derived>, Device>;
    Evaluator m_eval;

  public:
    selfCleaningEvaluator(const Evaluator& eval) : m_eval(eval)
    {
    }
    selfCleaningEvaluator(const Eigen::TensorBase<Derived, Eigen::ReadOnlyAccessors>& expr,
                          const Device&                                               device = Device())
        : m_eval(Evaluator(expr.eval(), device))
    {
        m_eval.evalSubExprsIfNeeded(nullptr);
    }

    ~selfCleaningEvaluator()
    {
        // This whole point of this object is to call cleanup automatically on destruct.
        // If there are pending operations to evaluate, m_eval will allocate a buffer to hold a result,
        // which needs to be deallocated.
        m_eval.cleanup();
    }

    constexpr auto
    rank()
    {
        using DimType = typename decltype(m_eval)::Dimensions::Base;
        return DimType{}.size(); // Because Derived::Dimensions is sometimes wrong
    }
    const Evaluator*
    operator->() const
    {
        return &m_eval;
    }
    Evaluator*
    operator->()
    {
        return &m_eval;
    }
    constexpr auto
    map()
    {
        // We inspect m_eval to get the type, rank and dimensions, because it has the resulting tensor,
        // whereas Derived is the tensor type _before_ whatever operation is pending (if any).
        using DimType       = typename decltype(m_eval)::Dimensions::Base;
        constexpr auto rank = DimType{}.size();
        using Scalar        = typename Eigen::internal::remove_const<typename decltype(m_eval)::Scalar>::type;
        return Eigen::TensorMap<Eigen::Tensor<Scalar, rank>>(m_eval.data(), m_eval.dimensions());
    }
};

// Evaluates expressions if needed
template <typename T, typename Device = Eigen::DefaultDevice>
auto
asEval(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>&
                     expr, // An Eigen::TensorBase object (Tensor, TensorMap, TensorExpr... )
       const Device& device = Device())
{ // Override to evaluate on another device, e.g. thread pool or gpu.
    return selfCleaningEvaluator(expr, device);
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Matrix with shape rows/cols
template <typename T, typename sizeType, typename Device = Eigen::DefaultDevice>
auto
MatrixCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const sizeType rows, const sizeType cols,
           const Device& device = Device())
{
    auto tensor    = asEval(expr, device);
    auto tensorMap = tensor.map();
    using Scalar   = typename decltype(tensorMap)::Scalar;
    return static_cast<MatrixType<Scalar>>(Eigen::Map<const MatrixType<Scalar>>(tensorMap.data(), rows, cols));
}

// Converts any Eigen::Tensor (or expression) to an Eigen::Matrix with shape rows/cols
template <typename T, typename Device = Eigen::DefaultDevice>
auto
VectorCast(const Eigen::TensorBase<T, Eigen::ReadOnlyAccessors>& expr, const Device& device = Device())
{
    auto tensor    = asEval(expr, device);
    auto tensorMap = tensor.map();
    auto size      = Eigen::internal::array_prod(tensorMap.dimensions());
    using Scalar   = typename decltype(tensorMap)::Scalar;
    return static_cast<VectorType<Scalar>>(Eigen::Map<const VectorType<Scalar>>(tensorMap.data(), size));
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