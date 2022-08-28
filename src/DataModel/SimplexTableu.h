#ifndef GMISOLVER_SIMPLEXTABLEU_H
#define GMISOLVER_SIMPLEXTABLEU_H

#include "src/DataModel/LinearProgram.h"
#include "src/Util/LPPrinter.h"

#include <iostream>

template <typename T> class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram,
                 const std::vector<T> &objectiveRow)
      : _variableInfos(linearProgram._variableInfos),
        _rowInfos(linearProgram._rowInfos),
        _constraintMatrix(linearProgram._constraintMatrix),
        _rightHandSides(linearProgram._rightHandSides) {
    _constraintMatrix.front() = objectiveRow;

    _basisMatrixInverse.resize(linearProgram._rowInfos.size());
    _basisMatrixInverse[0].resize(_basisMatrixInverse.size());
    _basisMatrixInverse[0][0] = 1;
    for (int columnIdx = 1; columnIdx < _basisMatrixInverse.size(); ++columnIdx)
      _basisMatrixInverse[0][columnIdx] =
          objectiveRow[linearProgram.initialVariableCountInStandardForm() +
                       columnIdx - 1] *
          _constraintMatrix[columnIdx]
                           [linearProgram.initialVariableCountInStandardForm() +
                            columnIdx - 1];

    for (int rowIdx = 1; rowIdx < _basisMatrixInverse.size(); ++rowIdx) {
      _basisMatrixInverse[rowIdx].resize(_basisMatrixInverse.size());
      _basisMatrixInverse[rowIdx][rowIdx] =
          _constraintMatrix[rowIdx]
                           [linearProgram.initialVariableCountInStandardForm() +
                            rowIdx - 1];
    }

    _reducedCosts.resize(_variableInfos.size());
    for (int columnIdx = 0; columnIdx < _variableInfos.size(); ++columnIdx) {
      T yAn = _constraintMatrix[0][columnIdx];
      for (int i = 1; i < _rowInfos.size(); ++i)
        yAn -= _basisMatrixInverse[0][i] * _constraintMatrix[i][columnIdx];

      _reducedCosts[columnIdx] = yAn;
    }
    _constraintMatrix[0] = _reducedCosts;

    auto &d = _rightHandSides[0];
    for (int i = 0; i < _rowInfos.size() - 1; ++i)
      d += _basisMatrixInverse[0][i + 1] * _rightHandSides[i + 1];

    int currentBasicIdx = 1;
    for (int j = 0; j < _variableInfos.size(); ++j) {
      if (_variableInfos[j]._isBasic)
        _rowToBasisColumnIdxMap[currentBasicIdx++] = j;
    }
  }

  std::string toString() const {
    LPPrinter lpPrinter(_variableInfos, _rowInfos);
    lpPrinter.printLineBreak();
    lpPrinter.printVariableInfos();
    lpPrinter.printLineBreak();
    lpPrinter.printMatrixWithRHS(_rowToBasisColumnIdxMap, _constraintMatrix,
                                 _rightHandSides);
    lpPrinter.printLineBreak();
    lpPrinter.printInverseBasisWithDual(_basisMatrixInverse);
    lpPrinter.printLineBreak();
    lpPrinter.printReducedCosts(_reducedCosts);
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
      std::cout << "VARIABLE IDX " << varIdx << ' '
                << _constraintMatrix[0][varIdx] << '\n';
      if (_constraintMatrix[0][varIdx] < 0.0 &&
          !_variableInfos[varIdx]._isBasic) {
        enteringVarIdx = varIdx;
        std::cout << "ENTERING VARIABLE IDX " << *enteringVarIdx << '\n';

        std::optional<int> leavingRowIdx;
        for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx)
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
          if (i != *leavingRowIdx) {
            const auto commonCoeff =
                _constraintMatrix[i][*enteringVarIdx] * pivotingTermInverse;
            for (int j = 0; j < _variableInfos.size(); ++j)
              _constraintMatrix[i][j] -=
                  commonCoeff * _constraintMatrix[*leavingRowIdx][j];

            //                        if (i == 0)
            //                            std::cout << "COMMON COEFF " <<
            //                            commonCoeff << " " <<
            //                            _rightHandSides[*leavingRowIdx] << " "
            //                             << _rightHandSides[i] << '\n';
            // ???
            if (i == 0)
              _rightHandSides[i] +=
                  commonCoeff * _rightHandSides[*leavingRowIdx];
            else
              _rightHandSides[i] -=
                  commonCoeff * _rightHandSides[*leavingRowIdx];
          }
        }

        for (int j = 0; j < _variableInfos.size(); ++j)
          _constraintMatrix[*leavingRowIdx][j] *= pivotingTermInverse;

        _rightHandSides[*leavingRowIdx] *= pivotingTermInverse;

        _variableInfos[*enteringVarIdx]._isBasic = true;
        std::cout << "LEAVING VARIABLE ROW IDX " << *leavingRowIdx << '\n';
        std::cout << "LEAVING VARIABLE COLUMN IDX "
                  << _rowToBasisColumnIdxMap[*leavingRowIdx] << '\n';
        _variableInfos[_rowToBasisColumnIdxMap[*leavingRowIdx]]._isBasic =
            false;
        _rowToBasisColumnIdxMap[*leavingRowIdx] = *enteringVarIdx;
        return false;
      }
    }

    _result = LPOptimizationResult::BOUNDED_AND_FEASIBLE;
    return true;
  }

private:
  std::vector<VariableInfo> _variableInfos;
  std::vector<RowInfo> _rowInfos;
  Matrix<T> _constraintMatrix;
  std::vector<T> _rightHandSides;

  std::vector<std::vector<T>> _basisMatrixInverse;
  std::vector<T> _reducedCosts;

  std::map<int, int> _rowToBasisColumnIdxMap;

  LPOptimizationResult _result;
};

#endif // GMISOLVER_SIMPLEXTABLEU_H
