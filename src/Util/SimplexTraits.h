#ifndef GMISOLVER_SIMPLEXTRAITS_H
#define GMISOLVER_SIMPLEXTRAITS_H

#include "src/DataModel/MatrixTypes.h"
#include "src/Util/SpdlogHeader.h"

#include <algorithm>
#include <cmath>
#include <numeric>

template <typename T> struct NumericalTraits {
  static_assert(!std::numeric_limits<T>::is_integer);

  constexpr static T PRIMAL_FEASIBILITY_TOLERANCE = 1e-6;
  constexpr static T DUAL_FEASIBILITY_TOLERANCE = 1e-6;
  constexpr static T OBJECTIVE_MONOTONICITY_TOLERANCE = 1e-6;
  constexpr static T INTEGRALITY_TOLERANCE = 1e-5;

  constexpr static T ABSOLUTE_EPSILON = 1e-9;
  constexpr static T PIVOT_ABSOLUTE_EPSILON = 1e-6;
  //  constexpr static T ABSOLUTE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static T RELATIVE_EPSILON = std::numeric_limits<T>::epsilon();
  constexpr static size_t ULP_COUNT = 4;

  static bool equal(const T &x, const T &y,
                    const T &absoluteTolerance = ABSOLUTE_EPSILON) {
    const T diff = std::fabs(x - y);
    return (diff < absoluteTolerance) ||
           (diff <= RELATIVE_EPSILON * std::fabs(x + y) * ULP_COUNT);
  }

  static bool greater(const T &x, const T &y,
                      const T &absoluteTolerance = ABSOLUTE_EPSILON) {
    return (x > y) && !equal(x, y, absoluteTolerance);
  }

  static bool less(const T &x, const T &y,
                   const T &absoluteTolerance = ABSOLUTE_EPSILON) {
    return (x < y) && !equal(x, y, absoluteTolerance);
  }

  static bool isEligibleForPivot(const T &x) {
    return std::fabs(x) > PIVOT_ABSOLUTE_EPSILON;
  }

  constexpr static T ABSOLUTE_TOLERANCE = 1e-14;

  static T addNormal(const T &x, const T &y) { return x + y; }

  struct NormalAddOp {
    T operator()(const T &x, const T &y) const { return addNormal(x, y); }
  };

  struct SafeNumericalAddOp {
    constexpr static T RELATIVE_TOLERANCE = 1e-10;

    T operator()(const T &x, const T &y) const {
      if (std::fabs(x + y) <
          RELATIVE_TOLERANCE * std::max(std::fabs(x), std::fabs(y)))
        return T{0.0};

      return x + y;
    }
  };

  static T add(const T &x, const T &y) { return addNormal(x, y); }

  template <typename AddOp> struct PositiveNegativeAdder {
    void addValue(const T &val) {
      if (val > 0.0)
        _sumPositives = _addOp(_sumPositives, val);
      else
        _sumNegatives = _addOp(_sumNegatives, val);
    }

    T currentSum() const { return _addOp(_sumPositives, _sumNegatives); }

    const AddOp _addOp{};
    T _sumPositives = 0.0;
    T _sumNegatives = 0.0;
  };

  template <typename AddOp> struct KahanAdder {
    void addValue(const T &val) {
      T t = _addOp(_sum, val);
      if (std::fabs(_sum) >= std::fabs(val))
        _c = _addOp(_c, _addOp(_addOp(_sum, -t), val));
      else
        _c = _addOp(_c, _addOp(_addOp(val, -t), _sum));

      _sum = t;
    }

    T currentSum() const { return _addOp(_sum, _c); }

    const AddOp _addOp{};
    T _sum = 0.0;
    T _c = 0.0;
  };

  template <typename AddOp = NormalAddOp> struct SimpleAdder {
    void addValue(const T &val) { _sum = _addOp(_sum, val); }

    T currentSum() const { return _sum; }

    const AddOp _addOp{};
    T _sum = 0.0;
  };

  static bool isZero(const T &x) { return std::fabs(x) < ABSOLUTE_TOLERANCE; }
};

template <typename T,
          MatrixRepresentationType representationType =
              MatrixRepresentationType::NORMAL,
          typename NumericalTraitsType = NumericalTraits<T>>
struct SimplexTraits {
  using NumericalTraitsT = NumericalTraitsType;
  using CurrentAdder = typename NumericalTraitsT::template SimpleAdder<
      typename NumericalTraitsT::NormalAddOp>;

  constexpr static bool useSparseRepresentationValue =
      representationType == MatrixRepresentationType::SPARSE;
  using ElementaryMatrixT =
      std::conditional_t<useSparseRepresentationValue,
                         SparseElementaryMatrix<T>, ElementaryMatrix<T>>;
  using VectorT = std::conditional_t<useSparseRepresentationValue,
                                     SparseVector<T>, std::vector<T>>;

  template <typename AdderT = CurrentAdder>
  static T dotProduct(const std::vector<T> &vec1, const std::vector<T> &vec2) {
    AdderT adder;

    for (int i = 0; i < vec1.size(); ++i)
      adder.addValue(vec1[i] * vec2[i]);

    return adder.currentSum();
  }

  template <typename AddOp = typename NumericalTraitsT::NormalAddOp>
  static void multiplyByETM(const ElementaryMatrix<T> &etm,
                            std::vector<T> &modifiedVec) {
    const T pivotModifiedVecTerm = modifiedVec[etm._pivotRowIdx];
    for (int i = 0; i < modifiedVec.size(); ++i) {
      if (i == etm._pivotRowIdx)
        modifiedVec[i] *= etm._pivotingTermInverse;
      else
        modifiedVec[i] =
            AddOp{}(modifiedVec[i], -(etm._pivotingTermInverse * etm._vec[i] *
                                      pivotModifiedVecTerm));
    }
  }

  template <typename AdderT = CurrentAdder>
  static void multiplyByETMFromRight(std::vector<T> &modifiedVec,
                                     const ElementaryMatrix<T> &etm) {
    AdderT adder;
    for (int i = 0; i < modifiedVec.size(); ++i) {
      if (i == etm._pivotRowIdx)
        adder.addValue(modifiedVec[i] * etm._pivotingTermInverse);
      else
        adder.addValue(-modifiedVec[i] * etm._vec[i] *
                       etm._pivotingTermInverse);
    }

    modifiedVec[etm._pivotRowIdx] = adder.currentSum();
  }

  template <typename AddOp = typename NumericalTraitsT::NormalAddOp>
  static void multiplyByETM(const ElementaryMatrix<T> &etm,
                            Matrix<T> &modifiedMatrix) {
    for (int rowIdx = 0; rowIdx < modifiedMatrix.size(); ++rowIdx) {
      if (rowIdx == etm._pivotRowIdx)
        continue;

      const auto commonCoeff = etm._vec[rowIdx] * etm._pivotingTermInverse;

      for (int colIdx = 0; colIdx < modifiedMatrix[rowIdx].size(); ++colIdx)
        modifiedMatrix[rowIdx][colIdx] =
            AddOp{}(modifiedMatrix[rowIdx][colIdx],
                    -(commonCoeff * modifiedMatrix[etm._pivotRowIdx][colIdx]));
    }

    for (int j = 0; j < modifiedMatrix[etm._pivotRowIdx].size(); ++j)
      modifiedMatrix[etm._pivotRowIdx][j] *= etm._pivotingTermInverse;
  }

  template <typename AdderT = CurrentAdder>
  static T dotProduct(const std::vector<T> &vecNormal,
                      const SparseVector<T> &vecSparse) {
    AdderT adder;
    for (const auto &elem : vecSparse._indexedValues) {
      if (elem._index >= vecNormal.size())
        break;

      adder.addValue(vecNormal[elem._index] * elem._data);
    }

    return adder.currentSum();
  }

  template <typename AddOp = typename NumericalTraitsT::NormalAddOp>
  static void multiplyByETM(const SparseElementaryMatrix<T> &sparseEtm,
                            std::vector<T> &modifiedVec) {
    const T pivotModifiedVecTerm = modifiedVec[sparseEtm._pivotRowIdx];
    for (const auto &etmElem : sparseEtm._sparseVec._indexedValues) {
      if (etmElem._index == sparseEtm._pivotRowIdx)
        modifiedVec[sparseEtm._pivotRowIdx] *= sparseEtm._pivotingTermInverse;
      else
        modifiedVec[etmElem._index] =
            AddOp{}(modifiedVec[etmElem._index],
                    -(sparseEtm._pivotingTermInverse * etmElem._data *
                      pivotModifiedVecTerm));
    }
  }

  template <typename AddOp = typename NumericalTraitsT::NormalAddOp>
  static void multiplyByETM(const SparseElementaryMatrix<T> &sparseEtm,
                            SparseVector<T> &modifiedVec) {
    const T pivotModifiedVecTerm =
        modifiedVec._normalVec[sparseEtm._pivotRowIdx];
    if (NumericalTraitsT::isZero(pivotModifiedVecTerm))
      return;

    auto currentModifiedVecIt = modifiedVec._indexedValues.begin();
    for (const auto &etmElem : sparseEtm._sparseVec._indexedValues) {
      while (currentModifiedVecIt != modifiedVec._indexedValues.end() &&
             currentModifiedVecIt->_index < etmElem._index)
        ++currentModifiedVecIt;

      if (currentModifiedVecIt != modifiedVec._indexedValues.end() &&
          currentModifiedVecIt->_index == etmElem._index) {
        if (etmElem._index == sparseEtm._pivotRowIdx) {
          currentModifiedVecIt->_data *= sparseEtm._pivotingTermInverse;
          modifiedVec._normalVec[currentModifiedVecIt->_index] =
              currentModifiedVecIt->_data;
        } else {
          const auto updatedValue =
              AddOp{}(currentModifiedVecIt->_data,
                      -sparseEtm._pivotingTermInverse * etmElem._data *
                          pivotModifiedVecTerm);
          if (NumericalTraitsT::isZero(updatedValue)) {
            modifiedVec._normalVec[etmElem._index] = T{0.0};
            currentModifiedVecIt =
                modifiedVec._indexedValues.erase(currentModifiedVecIt);
          } else {
            currentModifiedVecIt->_data =
                modifiedVec._normalVec[currentModifiedVecIt->_index] =
                    updatedValue;
          }
        }
      } else {
        currentModifiedVecIt = modifiedVec._indexedValues.insert(
            currentModifiedVecIt,
            {AddOp{}(0.0, -sparseEtm._pivotingTermInverse * etmElem._data *
                              pivotModifiedVecTerm),
             etmElem._index});
        modifiedVec._normalVec[etmElem._index] = currentModifiedVecIt->_data;
      }
    }
  }

  template <typename AdderT = CurrentAdder>
  static T dotProduct(const SparseVector<T> &vec1,
                      const SparseVector<T> &vec2) {
    AdderT adder;
    //    std::vector<T> addedValues;
    auto currentVec2It = vec2._indexedValues.begin();

    for (const auto &vec1Elem : vec1._indexedValues) {
      while ((currentVec2It != vec2._indexedValues.end()) &&
             (currentVec2It->_index < vec1Elem._index))
        ++currentVec2It;

      if (currentVec2It == vec2._indexedValues.end())
        break;

      if (currentVec2It->_index == vec1Elem._index) {
        adder.addValue(vec1Elem._data * currentVec2It->_data);
        //        addedValues.push_back(vec1Elem._data * currentVec2It->_data);
      }
    }

    //    if (std::fabs(adder.currentSum()) < 1e-10 &&
    //    std::fabs(adder.currentSum()) > 0.0)
    //    {
    //      SPDLOG_INFO("DOT NONZEROS {}", fmt::join(addedValues, ", "));
    //    }
    //      SPDLOG_INFO("DOT PRODUCT SPARSE {}", adder.currentSum());
    return adder.currentSum();
  }

  template <typename AdderT = CurrentAdder>
  static void multiplyByETMFromRight(SparseVector<T> &modifiedVec,
                                     const SparseElementaryMatrix<T> &etm) {
    AdderT adder;
    auto currentEtmIt = etm._sparseVec._indexedValues.begin();
    auto firstModifiecVecElemAfterPivotRowIt = modifiedVec._indexedValues.end();

    for (auto modifiecVecElemIt = modifiedVec._indexedValues.begin();
         modifiecVecElemIt != modifiedVec._indexedValues.end();
         ++modifiecVecElemIt) {
      if (modifiecVecElemIt->_index >= etm._pivotRowIdx &&
          firstModifiecVecElemAfterPivotRowIt ==
              modifiedVec._indexedValues.end())
        firstModifiecVecElemAfterPivotRowIt = modifiecVecElemIt;

      auto &modifiedVecElem = *modifiecVecElemIt;
      while ((currentEtmIt != etm._sparseVec._indexedValues.end()) &&
             (currentEtmIt->_index < modifiedVecElem._index))
        ++currentEtmIt;

      if (currentEtmIt == etm._sparseVec._indexedValues.end())
        break;

      if (currentEtmIt->_index == modifiedVecElem._index) {
        if (currentEtmIt->_index == etm._pivotRowIdx)
          adder.addValue(modifiedVecElem._data * etm._pivotingTermInverse);
        else
          adder.addValue(-modifiedVecElem._data * currentEtmIt->_data *
                         etm._pivotingTermInverse);
      }
    }

    const auto sum = adder.currentSum();
    if (NumericalTraitsT::isZero(sum)) {
      if ((firstModifiecVecElemAfterPivotRowIt !=
           modifiedVec._indexedValues.end()) &&
          (firstModifiecVecElemAfterPivotRowIt->_index == etm._pivotRowIdx)) {
        modifiedVec._indexedValues.erase(firstModifiecVecElemAfterPivotRowIt);
        modifiedVec._normalVec[etm._pivotRowIdx] = T{0.0};
      }
    } else {
      if ((firstModifiecVecElemAfterPivotRowIt !=
           modifiedVec._indexedValues.end()) &&
          (firstModifiecVecElemAfterPivotRowIt->_index == etm._pivotRowIdx)) {
        firstModifiecVecElemAfterPivotRowIt->_data =
            modifiedVec._normalVec[etm._pivotRowIdx] = sum;
      } else {
        modifiedVec._indexedValues.insert(firstModifiecVecElemAfterPivotRowIt,
                                          {sum, etm._pivotRowIdx});
        modifiedVec._normalVec[etm._pivotRowIdx] = sum;
      }
    }
  }

  template <typename AdderT = CurrentAdder>
  static void
  multiplyByETMFromRight(std::vector<T> &modifiedVec,
                         const SparseElementaryMatrix<T> &sparseEtm) {
    AdderT adder;
    for (const auto &etmElem : sparseEtm._sparseVec._indexedValues) {
      if (etmElem._index == sparseEtm._pivotRowIdx)
        adder.addValue(modifiedVec[sparseEtm._pivotRowIdx] *
                       sparseEtm._pivotingTermInverse);
      else
        adder.addValue(-modifiedVec[etmElem._index] * etmElem._data *
                       sparseEtm._pivotingTermInverse);
    }

    modifiedVec[sparseEtm._pivotRowIdx] = adder.currentSum();
  }
};

#endif // GMISOLVER_SIMPLEXTRAITS_H
