#ifndef GMISOLVER_SIMPLEXTRAITS_H
#define GMISOLVER_SIMPLEXTRAITS_H

#include "src/DataModel/MatrixTypes.h"

#include <algorithm>
#include <cmath>
#include <numeric>

template <typename T> struct SimplexTraits {
  static_assert(!std::numeric_limits<T>::is_integer);

  constexpr static T ABSOLUTE_EPSILON = 1e-9;
  constexpr static T PIVOT_ABSOLUTE_EPSILON = 1e-6;
  //  constexpr static T ABSOLUTE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static T RELATIVE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static size_t ULP_COUNT = 4;

  static bool equal(const T &x, const T &y) {
    const T diff = std::fabs(x - y);
    return (diff < ABSOLUTE_EPSILON) ||
           (diff <= RELATIVE_EPSILON * std::fabs(x + y) * ULP_COUNT);
  }

  static bool greater(const T &x, const T &y) {
    return (x > y) && !equal(x, y);
  }

  static bool less(const T &x, const T &y) { return (x < y) && !equal(x, y); }

  static bool isEligibleForPivot(const T &x) {
    return std::fabs(x) > PIVOT_ABSOLUTE_EPSILON;
  }

  constexpr static T RELATIVE_TOLERANCE = 1e-10;
  constexpr static T ABSOLUTE_TOLERANCE = 1e-14;

  static T addSafeNumerical(const T& x, const T& y)
  {
    if (std::fabs(x + y) < RELATIVE_TOLERANCE * std::max(std::fabs(x), std::fabs(y)))
      return T{0.0};

    return x + y;
  }

  static T addNormal(const T& x, const T& y)
  {
    return x + y;
  }

  static T dotProduct(const std::vector<T>& vecNormal, const SparseVector<T>& vecSparse)
  {
    T result{};

    for (const auto& elem : vecSparse._elements)
    {
      if (elem._index >= vecNormal.size())
        break;
      result += vecNormal[elem._index] * elem._data;
    }

    return result;
  }

  static void multiplyByETM(const ElementaryMatrix<T>& etm, std::vector<T>& modifiedVec)
  {
    const T pivotModifiedVecTerm = modifiedVec[etm._pivotIdx];
    for (int i = 0; i < modifiedVec.size(); ++i) {
      if (i == etm._pivotIdx)
        modifiedVec[i] *= etm._pivotingTermInverse;
      else
        modifiedVec[i] -= etm._pivotingTermInverse * etm._vec[i] * pivotModifiedVecTerm;
    }
  }


  static void multiplyByETM(const ElementaryMatrix<T>& etm, Matrix<T>&modifiedMatrix)
  {
    for (int rowIdx = 0; rowIdx < modifiedMatrix.size(); ++rowIdx) {
      if (rowIdx == etm._pivotIdx)
        continue;

      const auto commonCoeff = etm._vec[rowIdx] * etm._pivotingTermInverse;

      for (int colIdx = 0; colIdx < modifiedMatrix[rowIdx].size(); ++colIdx)
        modifiedMatrix[rowIdx][colIdx] -= commonCoeff * modifiedMatrix[etm._pivotIdx][colIdx];
    }

    for (int j = 0; j < modifiedMatrix[etm._pivotIdx].size(); ++j)
      modifiedMatrix[etm._pivotIdx][j] *= etm._pivotingTermInverse;

//    for (int i = 0; i < _rowInfos.size(); ++i) {
//      if (i == leavingRowIdx)
//        continue;
//      const auto commonCoeff = enteringColumn[i] * pivotingTermInverse;
//
//      for (int j = 0; j < _basisMatrixInverse.size(); ++j)
//        _basisMatrixInverse[i][j] -=
//            commonCoeff * _basisMatrixInverse[leavingRowIdx][j];
//    }
//
//    for (int j = 0; j < _basisMatrixInverse.size(); ++j)
//      _basisMatrixInverse[leavingRowIdx][j] *= pivotingTermInverse;
  }
};

template <typename T> struct SimpleComparisonTraits {
  static bool equal(const T &x, const T &y) { return x == y; }

  static bool less(const T &x, const T &y) { return x < y; }

  static bool greater(const T &x, const T &y) { return x > y; }
};

#endif // GMISOLVER_SIMPLEXTRAITS_H
