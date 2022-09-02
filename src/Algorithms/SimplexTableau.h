#ifndef GMISOLVER_SIMPLEXTABLEAU_H
#define GMISOLVER_SIMPLEXTABLEAU_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/LinearProgram.h"
#include "src/DataModel/SimplexBasisData.h"

#include <iostream>

template <typename T> class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram);

  void addArtificialVariables();
  void init();

  std::vector<T> artificialObjectiveRow() const;
  void calculateCurrentObjectiveValue();
  void calculateSolution();

  std::string toString() const;
  std::string toStringShort() const;
  std::string toStringLpSolveFormat() const;

  const std::vector<VariableInfo> &getVariableInfos() const {
    return _variableInfos;
  }
  const std::vector<RowInfo> &getRowInfos() const { return _rowInfos; }


private:
  template <typename U, typename ComparisonTraitsT>
  friend class PrimalSimplex;

  template <typename U, typename ComparisonTraitsT>
  friend class RevisedPrimalSimplexPFI;

  std::optional<SimplexBasisData> createBasisFromArtificialVars() const;

  void initBasisMatrixInverse();
  void initDual();

  void calculateReducedCosts();
  void updateBasisData(const PivotData<T>& pivotData);

  const LinearProgram<T>& _initialProgram;
  std::vector<VariableInfo> _variableInfos;
  std::vector<RowInfo> _rowInfos;

  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;
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
