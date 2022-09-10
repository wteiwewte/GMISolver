#ifndef GMISOLVER_SIMPLEXTABLEAU_H
#define GMISOLVER_SIMPLEXTABLEAU_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/LinearProgram.h"
#include "src/DataModel/SimplexBasisData.h"

#include <iostream>

template <typename T> class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram,
                 const bool isPrimalSimplex);

  void convertToStandardForm();
  void makeRightHandSidesNonNegative();
  void addBoundsToMatrix();
  void addArtificialVariables();
  void init(const bool isPrimalSimplex);

  std::vector<T> initialPrimalSimplexObjective() const;
  std::vector<T> initialDualSimplexObjective() const;
  void calculateCurrentObjectiveValue();
  void calculateSolution();

  std::string toString() const;
  std::string toStringShort() const;
  std::string toStringShortWithSolution() const;
  std::string toStringLpSolveFormat() const;

  const std::vector<VariableInfo> &getVariableInfos() const {
    return _variableInfos;
  }
  const std::vector<RowInfo> &getRowInfos() const { return _rowInfos; }

private:
  template <typename U, typename ComparisonTraitsT> friend class PrimalSimplex;

  template <typename U, typename ComparisonTraitsT>
  friend class RevisedPrimalSimplexPFI;
  template <typename U, typename ComparisonTraitsT>
  friend class RevisedPrimalSimplexPFIBounds;
  template <typename U, typename ComparisonTraitsT>
  friend class RevisedDualSimplexPFIBounds;

  std::optional<SimplexBasisData>
  createBasisFromArtificialVars() const;

  int basicColumnIdx(const int rowIdx) const { return _simplexBasisData._rowToBasisColumnIdxMap[rowIdx]; }
  bool isColumnAllowedToEnterBasis(const int colIdx);
  std::optional<T> curSatisfiedBound(const int varIdx);

  std::vector<T> computeTableauColumn(const int enteringColumnIdx);
  std::vector<T> computeTableauRow(const int rowIdx);

  void pivot(
      const int rowIdx, const int enteringColumnIdx,
      const std::vector<T> &enteringColumn, const std::vector<T> &pivotRow);
  void updateReducedCosts(const PivotData<T> &pivotData,
                          const std::vector<T> &pivotRow);
  void updateInverseMatrixWithRHS(const PivotData<T> &pivotData,
                                  const std::vector<T> &enteringColumn);

  void initBasisMatrixInverse();
  void initDual();
  void initBoundsForDualSimplex();

  void calculateReducedCostsBasedOnDual();
  void updateBasisData(const PivotData<T> &pivotData);

  const LinearProgram<T> &_initialProgram;
  std::vector<VariableInfo> _variableInfos;
  std::vector<std::optional<T>> _variableLowerBounds;
  std::vector<std::optional<T>> _variableUpperBounds;
  std::set<std::string> _variableLabelSet;
  std::vector<RowInfo> _rowInfos;

  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
  std::vector<T> _initialRightHandSides;
  std::vector<T> _initialObjectiveRow;
  std::vector<T> _objectiveRow;

  std::vector<std::vector<T>> _basisMatrixInverse;
  std::vector<T> _reducedCosts;
  std::vector<T> _y;
  std::vector<T> _x;
  T _objectiveValue{};

  SimplexBasisData _simplexBasisData;

  LPOptimizationResult _result;
};

#endif // GMISOLVER_SIMPLEXTABLEAU_H
