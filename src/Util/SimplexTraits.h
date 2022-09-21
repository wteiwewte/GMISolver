#ifndef GMISOLVER_SIMPLEXTRAITS_H
#define GMISOLVER_SIMPLEXTRAITS_H

#include "src/DataModel/MatrixTypes.h"
#include "src/Util/SpdlogHeader.h"

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

  static T add(const T& x, const T& y)
  {
    return addNormal(x, y);
  }

  struct PositiveNegativeAdder
  {
    void addValue(const T& val)
    {
      if (val > 0.0)
        _sumPositives = add(_sumPositives, val);
      else
        _sumNegatives = add(_sumNegatives, val);
    }

    T currentSum() const
    {
      return add(_sumPositives, _sumNegatives);
    }

    T _sumPositives = 0.0;
    T _sumNegatives = 0.0;
  };

  struct KahanAdder {
    void addValue(const T& val)
    {
      T t = add(_sum, val);
      if (std::fabs(_sum) >= std::fabs(val))
        _c = add(_c, (_sum - t) + val);
      else
        _c = add(_c, (val - t) + _sum);

      _sum = t;
    }

    T currentSum() const
    {
      return add(_sum, _c);
    }

    T _sum = 0.0;
    T _c = 0.0;
  };

  struct SimpleAdder {
    void addValue(const T& val)
    {
      _sum = add(_sum, val);
    }

    T currentSum() const
    {
      return _sum;
    }

    T _sum = 0.0;
  };

  using Adder = SimpleAdder;

  static T dotProduct(const std::vector<T>& vecNormal, const SparseVector<T>& vecSparse)
  {
    Adder adder;
    for (const auto& elem : vecSparse._elements)
    {
      if (elem._index >= vecNormal.size())
        break;

      adder.addValue(vecNormal[elem._index] * elem._data);
    }

    return adder.currentSum();
  }

  static T dotProduct(const std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    Adder adder;

    for (int i = 0; i < vec1.size(); ++i)
      adder.addValue(vec1[i] * vec2[i]);

    return adder.currentSum();
  }

  static void multiplyByETM(const ElementaryMatrix<T>& etm, std::vector<T>& modifiedVec)
  {
    const T pivotModifiedVecTerm = modifiedVec[etm._pivotRowIdx];
    for (int i = 0; i < modifiedVec.size(); ++i) {
      if (i == etm._pivotRowIdx)
        modifiedVec[i] *= etm._pivotingTermInverse;
      else
        modifiedVec[i] = add(modifiedVec[i], -(etm._pivotingTermInverse * etm._vec[i] * pivotModifiedVecTerm));
    }
  }

  static void multiplyByETMSparse(const SparseElementaryMatrix<T>&sparseEtm, std::vector<T>& modifiedVec)
  {
    const T pivotModifiedVecTerm = modifiedVec[sparseEtm._pivotRowIdx];
    for (const auto&etmElem : sparseEtm._sparseVec._elements) {
      if (etmElem._index == sparseEtm._pivotRowIdx)
        modifiedVec[sparseEtm._pivotRowIdx] *= sparseEtm._pivotingTermInverse;
      else
        modifiedVec[etmElem._index] = add(modifiedVec[etmElem._index], -(sparseEtm._pivotingTermInverse *
                                               etmElem._data * pivotModifiedVecTerm));
    }
  }

  static void multiplyByETMSparse(const SparseElementaryMatrix<T>&sparseEtm, SparseVector<T>& modifiedVec)
  {
//    SPDLOG_INFO("{} {}", modifiedVec._normalVec.size(), sparseEtm._pivotRowIdx);

    const T pivotModifiedVecTerm = modifiedVec._normalVec[sparseEtm._pivotRowIdx];
    if (pivotModifiedVecTerm == T{0.0})
      return;

    auto currentModifiedVecIt = modifiedVec._elements.begin();
    for (const auto&etmElem : sparseEtm._sparseVec._elements) {
      while (currentModifiedVecIt != modifiedVec._elements.end() && currentModifiedVecIt->_index < etmElem._index)
        ++currentModifiedVecIt;

      if (currentModifiedVecIt != modifiedVec._elements.end() && currentModifiedVecIt->_index == etmElem._index)
      {
        if (etmElem._index == sparseEtm._pivotRowIdx)
        {
          currentModifiedVecIt->_data *= sparseEtm._pivotingTermInverse;
          modifiedVec._normalVec[currentModifiedVecIt->_index] = currentModifiedVecIt->_data;
        }
        else
        {
          const auto updatedValue =
              add(currentModifiedVecIt->_data, -sparseEtm._pivotingTermInverse * etmElem._data * pivotModifiedVecTerm);
//          if (std::fabs(updatedValue) < 1e-10)
//            SPDLOG_INFO("{}, {}", add(currentModifiedVecIt->_data, -sparseEtm._pivotingTermInverse * etmElem._data * pivotModifiedVecTerm),
//                        addSafeNumerical(currentModifiedVecIt->_data, -sparseEtm._pivotingTermInverse * etmElem._data * pivotModifiedVecTerm));
          // TODO addSafeNumerical ?
          if (updatedValue == T{0.0})
          {
            modifiedVec._normalVec[etmElem._index] = T{0.0};
            currentModifiedVecIt = modifiedVec._elements.erase(currentModifiedVecIt);
          }
          else
          {
            currentModifiedVecIt->_data = modifiedVec._normalVec[currentModifiedVecIt->_index]
                = updatedValue;
          }
        }
      }
      else
      {
        currentModifiedVecIt = modifiedVec._elements.insert(currentModifiedVecIt,
                                                            {add(0.0, -sparseEtm._pivotingTermInverse * etmElem._data * pivotModifiedVecTerm),
                                                            etmElem._index});
        modifiedVec._normalVec[etmElem._index] = currentModifiedVecIt->_data;
      }
    }
  }


  static T dotProductSparse(const SparseVector<T>& vec1, const SparseVector<T>& vec2)
  {
    Adder adder;
    auto currentVec2It = vec2._elements.begin();

    for (const auto& vec1Elem : vec1._elements)
    {
      while ((currentVec2It != vec2._elements.end()) && (currentVec2It->_index < vec1Elem._index))
        ++currentVec2It;

      if (currentVec2It == vec2._elements.end())
        break;

      if (currentVec2It->_index == vec1Elem._index)
        adder.addValue(vec1Elem._data * currentVec2It->_data);
    }

    return adder.currentSum();
  }


  static void multiplyByETMFromRightSparse(SparseVector<T>& modifiedVec, const SparseElementaryMatrix<T>& etm)
    {
      Adder adder;
      auto currentEtmIt = etm._sparseVec._elements.begin();
      auto firstModifiecVecElemAfterPivotRowIt = modifiedVec._elements.end();

      for (auto modifiecVecElemIt = modifiedVec._elements.begin();
               modifiecVecElemIt != modifiedVec._elements.end(); ++modifiecVecElemIt)
      {
        if (modifiecVecElemIt->_index >= etm._pivotRowIdx &&
                                         firstModifiecVecElemAfterPivotRowIt == modifiedVec._elements.end())
          firstModifiecVecElemAfterPivotRowIt = modifiecVecElemIt;

        auto& modifiedVecElem = *modifiecVecElemIt;
        while ((currentEtmIt != etm._sparseVec._elements.end()) && (currentEtmIt->_index < modifiedVecElem._index))
          ++currentEtmIt;

        if (currentEtmIt == etm._sparseVec._elements.end())
          break;

        if (currentEtmIt->_index == modifiedVecElem._index)
        {
          if (currentEtmIt->_index == etm._pivotRowIdx)
            adder.addValue(modifiedVecElem._data * etm._pivotingTermInverse);
          else
            adder.addValue(-modifiedVecElem._data * currentEtmIt->_data * etm._pivotingTermInverse);
        }
      }

      const auto sum = adder.currentSum();
      if (sum == T{0.0})
      {
        if ((firstModifiecVecElemAfterPivotRowIt != modifiedVec._elements.end()) &&
        (firstModifiecVecElemAfterPivotRowIt->_index == etm._pivotRowIdx))
        {
          modifiedVec._elements.erase(firstModifiecVecElemAfterPivotRowIt);
          modifiedVec._normalVec[etm._pivotRowIdx] = T{0.0};
        }
      }
      else
      {
        if ((firstModifiecVecElemAfterPivotRowIt != modifiedVec._elements.end())  &&
            (firstModifiecVecElemAfterPivotRowIt->_index == etm._pivotRowIdx))
        {
          firstModifiecVecElemAfterPivotRowIt->_data =
              modifiedVec._normalVec[etm._pivotRowIdx] = sum;
        }
        else
        {
          modifiedVec._elements.insert(firstModifiecVecElemAfterPivotRowIt, {sum, etm._pivotRowIdx});
          modifiedVec._normalVec[etm._pivotRowIdx] = sum;
        }
      }
    }

  static void multiplyByETMFromRight(std::vector<T>& modifiedVec, const ElementaryMatrix<T>& etm)
  {
    Adder adder;
    for (int i = 0; i < modifiedVec.size(); ++i) {
      if (i == etm._pivotRowIdx)
        adder.addValue(modifiedVec[i] * etm._pivotingTermInverse);
      else
        adder.addValue(-modifiedVec[i] * etm._vec[i] * etm._pivotingTermInverse);
    }

    modifiedVec[etm._pivotRowIdx] = adder.currentSum();
  }

  static void multiplyByETMFromRightSparse(std::vector<T>& modifiedVec, const SparseElementaryMatrix<T>&sparseEtm)
  {
    Adder adder;
    for (const auto&etmElem : sparseEtm._sparseVec._elements) {
      if (etmElem._index == sparseEtm._pivotRowIdx)
        adder.addValue(modifiedVec[sparseEtm._pivotRowIdx] * sparseEtm._pivotingTermInverse);
      else
        adder.addValue(-modifiedVec[etmElem._index] * etmElem._data * sparseEtm._pivotingTermInverse);
    }

    modifiedVec[sparseEtm._pivotRowIdx] = adder.currentSum();
  }


  static void multiplyByETM(const ElementaryMatrix<T>& etm, Matrix<T>& modifiedMatrix)
  {
    for (int rowIdx = 0; rowIdx < modifiedMatrix.size(); ++rowIdx) {
      if (rowIdx == etm._pivotRowIdx)
        continue;

      const auto commonCoeff = etm._vec[rowIdx] * etm._pivotingTermInverse;

      for (int colIdx = 0; colIdx < modifiedMatrix[rowIdx].size(); ++colIdx)
        modifiedMatrix[rowIdx][colIdx] =
            add(modifiedMatrix[rowIdx][colIdx], -(commonCoeff * modifiedMatrix[etm._pivotRowIdx][colIdx]));
    }

    for (int j = 0; j < modifiedMatrix[etm._pivotRowIdx].size(); ++j)
      modifiedMatrix[etm._pivotRowIdx][j] *= etm._pivotingTermInverse;
  }
};

template <typename T> struct SimpleComparisonTraits {
  static bool equal(const T &x, const T &y) { return x == y; }

  static bool less(const T &x, const T &y) { return x < y; }

  static bool greater(const T &x, const T &y) { return x > y; }
};

#endif // GMISOLVER_SIMPLEXTRAITS_H



