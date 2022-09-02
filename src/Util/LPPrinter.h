#ifndef GMISOLVER_LPPRINTER_H
#define GMISOLVER_LPPRINTER_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/SimplexBasisData.h"

#include <map>
#include <spdlog/spdlog.h>
// fmt must be included after spdlog - weird bug
#include <fmt/format.h>
#include <functional>
#include <optional>
#include <sstream>
#include <string>

struct LPPrinter {
  constexpr static int CONSTRAINT_SIGN_WIDTH = 6;
  constexpr static int COEFFICIENT_WIDTH = 12;

  LPPrinter(const std::vector<VariableInfo> &variableInfos,
            const std::vector<RowInfo> &rowInfos)
      : _variableInfos(variableInfos), _rowInfos(rowInfos),
        _variableWidths(variableInfos.size()) {
    _totalWidth =
        CONSTRAINT_SIGN_WIDTH + COEFFICIENT_WIDTH + 2 + variableInfos.size();
    for (int variableIdx = 0; variableIdx < variableInfos.size();
         ++variableIdx) {
      _totalWidth +=
          (_variableWidths[variableIdx] = std::max<int>(
               variableInfos[variableIdx]._label.size(), COEFFICIENT_WIDTH));
      _maxVariableWidth =
          std::max(_maxVariableWidth, _variableWidths[variableIdx]);
    }

    _oss << '\n';
  }

  void printVariableInfos(
      const std::optional<std::reference_wrapper<const SimplexBasisData>>
          simplexBasisDataRef) {
    _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx) {
      _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx]._label,
                          _variableWidths[variableIdx]);
    }
    _oss << '\n';

    _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx) {
      _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx].typeStr(),
                          _variableWidths[variableIdx]);
    }
    _oss << '\n';

    if (simplexBasisDataRef.has_value()) {
      _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
      for (int variableIdx = 0; variableIdx < _variableInfos.size();
           ++variableIdx)
        _oss << fmt::format(
            "{:^{}}|",
            simplexBasisDataRef->get()._isBasicColumnIndexBitset[variableIdx]
                ? "BASIC"
                : "NON-BASIC",
            _variableWidths[variableIdx]);
      _oss << '\n';
    }
  }

  template <typename T>
  void printReducedCostWithObjectiveValue(const std::vector<T> &reducedCosts,
                                          const T &objectiveValue) {
    _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx)
      _oss << fmt::format("{:>{}}|",
                          (' ' + std::to_string(reducedCosts[variableIdx])),
                          _variableWidths[variableIdx]);

    _oss << fmt::format("{:^{}}|", "OBJ", CONSTRAINT_SIGN_WIDTH);
    _oss << fmt::format("{:>{}}|\n", (' ' + std::to_string(objectiveValue)),
                        COEFFICIENT_WIDTH);
  }

  template <typename T>
  void printMatrixWithRHS(const std::vector<int> &rowToBasisColumnIdxMap,
                          const Matrix<T> &matrix,
                          const std::vector<T> &rightHandSides) {
    for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
      if (rowIdx >= rowToBasisColumnIdxMap.size())
        _oss << fmt::format("{:^{}}|", "NOT FOUND", _maxVariableWidth);
      else
        _oss << fmt::format(
            "{:^{}}|", _variableInfos[rowToBasisColumnIdxMap[rowIdx]]._label,
            _maxVariableWidth);

      for (int variableIdx = 0; variableIdx < _variableInfos.size();
           ++variableIdx)
        _oss << fmt::format("{:>{}}|",
                            (' ' + std::to_string(matrix[rowIdx][variableIdx])),
                            _variableWidths[variableIdx]);

      _oss << fmt::format("{:^{}}|", rowTypeToStr(_rowInfos[rowIdx]._type),
                          CONSTRAINT_SIGN_WIDTH);
      _oss << fmt::format("{:>{}}|\n",
                          (' ' + std::to_string(rightHandSides[rowIdx])),
                          COEFFICIENT_WIDTH);
    }
  }

  template <typename T>
  void printInverseBasis(const Matrix<T> &basisMatrixInverse) {
    _oss << "BASIS MATRIX INVERSE\n";
    for (int rowIdx = 0; rowIdx < basisMatrixInverse.size(); ++rowIdx) {
      for (int variableIdx = 0; variableIdx < basisMatrixInverse.size();
           ++variableIdx)
        _oss << fmt::format(
            "{:>{}}|",
            (' ' + std::to_string(basisMatrixInverse[rowIdx][variableIdx])),
            COEFFICIENT_WIDTH);
      _oss << '\n';
    }
  }
  template <typename T> void printDual(const std::vector<T> &y) {
    _oss << "DUAL\n";
    for (int i = 0; i < y.size(); ++i) {
      _oss << fmt::format("{:>{}}|", (' ' + std::to_string(y[i])),
                          COEFFICIENT_WIDTH);
    }
    _oss << '\n';
  }

  template <typename T> void printSolution(const std::vector<T> &x) {
    _oss << "SOLUTION\n";
    _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
    for (int variableIdx = 0; variableIdx < _variableInfos.size();
         ++variableIdx) {
      _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx]._label,
                          _variableWidths[variableIdx]);
    }
    _oss << '\n';
    _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
    for (int i = 0; i < x.size(); ++i) {
      _oss << fmt::format("{:>{}}|", (' ' + std::to_string(x[i])),
                          COEFFICIENT_WIDTH);
    }
    _oss << '\n';
  }

  void printLineBreak() {
    _oss.width(_totalWidth);
    _oss.fill('-');
    _oss << '-' << '\n';
    _oss.fill(' ');
  }

  template <typename T>
  void printInLpSolveFormat(const Matrix<T> &matrix,
                            const std::vector<T> &objective,
                            const std::vector<T> &rightHandSides) {
    _oss << "min: ";
    for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
      if (objective[varIdx] != 0.0)
        _oss << fmt::format(" {:+f} {}", objective[varIdx],
                            _variableInfos[varIdx]._label);
    }
    _oss << ";\n";

    for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
      for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx) {
        if (matrix[rowIdx][varIdx] != 0.0)
          _oss << fmt::format(" {:+f} {}", matrix[rowIdx][varIdx],
                              _variableInfos[varIdx]._label);
      }
      _oss << ' ' << rowTypeToStr(_rowInfos[rowIdx]._type);
      _oss << fmt::format(" {:+f}", rightHandSides[rowIdx]) << ";\n";
    }
  }

  std::string toString() const { return _oss.str(); }

  const std::vector<VariableInfo> &_variableInfos;
  const std::vector<RowInfo> &_rowInfos;
  std::vector<int> _variableWidths;
  std::ostringstream _oss;
  int _totalWidth{0};
  int _maxVariableWidth{0};
};

#endif // GMISOLVER_LPPRINTER_H
