#ifndef GMISOLVER_SIMPLEXTABLEU_H
#define GMISOLVER_SIMPLEXTABLEU_H

#include "SimplexBasisData.h"
#include "src/DataModel/LinearProgram.h"
#include "src/Util/LPPrinter.h"

#include <iostream>

template <typename T> class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram,
                 const std::vector<T> &objectiveRow)
      : _initialProgram(linearProgram), _initialObjectiveRow(objectiveRow),
        _variableInfos(linearProgram._variableInfos),
        _rowInfos(linearProgram._rowInfos),
        _constraintMatrix(linearProgram._constraintMatrix),
        _rightHandSides(linearProgram._rightHandSides) {
    _basisMatrixInverse.resize(linearProgram._rowInfos.size());
    _y.resize(_basisMatrixInverse.size());
    for (int columnIdx = 0; columnIdx < _basisMatrixInverse.size(); ++columnIdx)
      _y[columnIdx] =
          objectiveRow[linearProgram.initialVariableCountInStandardForm() +
                       columnIdx] *
          _constraintMatrix[columnIdx]
                           [linearProgram.initialVariableCountInStandardForm() +
                            columnIdx];

    for (int rowIdx = 0; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
      _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
      _basisMatrixInverse[rowIdx][rowIdx] =
          _constraintMatrix[rowIdx]
                           [linearProgram.initialVariableCountInStandardForm() +
                            rowIdx];
    }

    _reducedCosts.resize(_variableInfos.size());
    for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx) {
      T yAn = objectiveRow[columnIdx];
      for (int i = 0; i < _rowInfos.size(); ++i)
        yAn -= _y[i] * _constraintMatrix[i][columnIdx];

      _reducedCosts[columnIdx] = yAn;
    }
    auto simplexBasisData = linearProgram.createBasisFromArtificialVars();
    if (simplexBasisData.has_value())
      _simplexBasisData = std::move(*simplexBasisData);

    calculateCurrentObjectiveValue();
  }

  std::string toString() const {
    LPPrinter lpPrinter(_variableInfos, _rowInfos);
    lpPrinter.printLineBreak();
    lpPrinter.printVariableInfos(std::cref(_simplexBasisData));
    lpPrinter.printLineBreak();
    lpPrinter.printReducedCostWithObjectiveValue(_reducedCosts,
                                                 _objectiveValue);
    lpPrinter.printMatrixWithRHS(_simplexBasisData._rowToBasisColumnIdxMap,
                                 _constraintMatrix, _rightHandSides);
    lpPrinter.printLineBreak();
    lpPrinter.printInverseBasisWithDual(_basisMatrixInverse);
    lpPrinter.printLineBreak();

    return lpPrinter.toString();
  }

  std::string toStringLpSolveFormat() const {
    LPPrinter lpPrinter(_variableInfos, _rowInfos);
    lpPrinter.printInLpSolveFormat(_constraintMatrix, _rightHandSides);
    return lpPrinter.toString();
  }

  bool primalSimplex() {
    std::optional<int> enteringVarIdx;
    for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
      //            if (_reducedCosts[varIdx] < 0.0 &&
      //            _variableInfos[varIdx]._isBasic)
      //                spdlog::warn
      spdlog::debug("VARIABLE IDX {} {}", varIdx, _reducedCosts[varIdx]);
      if (_reducedCosts[varIdx] < 0.0 &&
          !_simplexBasisData._isBasicColumnIndexBitset[varIdx]) {
        enteringVarIdx = varIdx;

        spdlog::debug("ENTERING VARIABLE IDX {}", *enteringVarIdx);

        std::optional<int> leavingRowIdx;
        for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx)
          if (_constraintMatrix[rowIdx][*enteringVarIdx] > 0.0)
            if (!leavingRowIdx.has_value() ||
                _rightHandSides[rowIdx] *
                        _constraintMatrix[*leavingRowIdx][varIdx] <
                    _rightHandSides[*leavingRowIdx] *
                        _constraintMatrix[rowIdx][varIdx])
              leavingRowIdx = rowIdx;

        if (!leavingRowIdx.has_value()) {
          _result = LPOptimizationResult::UNBOUNDED;
          return true;
        }

        const auto pivotingTermInverse =
            1.0 / _constraintMatrix[*leavingRowIdx][*enteringVarIdx];

        for (int i = 0; i < _rowInfos.size(); ++i) {
          if (i == leavingRowIdx)
            continue;

          const auto commonCoeff =
              _constraintMatrix[i][*enteringVarIdx] * pivotingTermInverse;

          for (int j = 0; j < _variableInfos.size(); ++j)
            _constraintMatrix[i][j] -=
                commonCoeff * _constraintMatrix[*leavingRowIdx][j];

          _rightHandSides[i] -= commonCoeff * _rightHandSides[*leavingRowIdx];
        }

        const auto commonCoeffReducedCost =
            _reducedCosts[*enteringVarIdx] * pivotingTermInverse;
        for (int j = 0; j < _variableInfos.size(); ++j)
          _reducedCosts[j] -=
              commonCoeffReducedCost * _constraintMatrix[*leavingRowIdx][j];
        //        _objectiveValue += commonCoeffReducedCost * _objectiveValue;

        for (int j = 0; j < _variableInfos.size(); ++j)
          _constraintMatrix[*leavingRowIdx][j] *= pivotingTermInverse;

        _rightHandSides[*leavingRowIdx] *= pivotingTermInverse;


        auto &[rowToBasisColumnIdxMap, isBasicColumnIndexBitset] =
            _simplexBasisData;
        spdlog::debug("LEAVING VARIABLE ROW IDX {} COLUMN IDX", *leavingRowIdx,
                      rowToBasisColumnIdxMap[*leavingRowIdx]);

        isBasicColumnIndexBitset[*enteringVarIdx] = true;
        isBasicColumnIndexBitset[rowToBasisColumnIdxMap[*leavingRowIdx]] =
            false;
        rowToBasisColumnIdxMap[*leavingRowIdx] = *enteringVarIdx;
        calculateCurrentObjectiveValue();

        return false;
      }
    }

    _result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

private:
  void calculateCurrentObjectiveValue() {
    _objectiveValue = T{};
    for (int i = 0; i < _rowInfos.size(); ++i)
      _objectiveValue +=
          _rightHandSides[i] *
          _initialObjectiveRow[_simplexBasisData._rowToBasisColumnIdxMap[i]];
  }

  const LinearProgram<T>& _initialProgram;
  const std::vector<T> _initialObjectiveRow;

  std::vector<VariableInfo> _variableInfos;
  std::vector<RowInfo> _rowInfos;
  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;

  std::vector<std::vector<T>> _basisMatrixInverse;
  std::vector<T> _reducedCosts;
  std::vector<T> _y;
  T _objectiveValue{};

  SimplexBasisData _simplexBasisData;

  LPOptimizationResult _result;
};

#endif // GMISOLVER_SIMPLEXTABLEU_H
