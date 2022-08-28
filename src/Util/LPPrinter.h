#ifndef GMISOLVER_LPPRINTER_H
#define GMISOLVER_LPPRINTER_H

#include "src/DataModel/CommonTypes.h"

#include <map>
#include <string>
#include <sstream>
#include <spdlog/spdlog.h>
#include <fmt/format.h>

struct LPPrinter
{
    constexpr static int CONSTRAINT_SIGN_WIDTH = 6;
    constexpr static int COEFFICIENT_WIDTH = 12;

    LPPrinter(const std::vector<VariableInfo>& variableInfos,
                       const std::vector<RowInfo>& rowInfos) :
    _variableInfos(variableInfos),
    _rowInfos(rowInfos),
        _variableWidths(variableInfos.size())
    {
        _totalWidth = CONSTRAINT_SIGN_WIDTH + COEFFICIENT_WIDTH + 2 + variableInfos.size();
        for (int variableIdx = 0; variableIdx < variableInfos.size(); ++ variableIdx)
        {
            _totalWidth += (_variableWidths[variableIdx] = std::max<int>(variableInfos[variableIdx]._label.size(), COEFFICIENT_WIDTH));
            _maxVariableWidth = std::max(_maxVariableWidth, _variableWidths[variableIdx]);
        }
    }

    void printVariableInfos()
    {
        _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
        for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++ variableIdx) {
            _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx]._label, _variableWidths[variableIdx]);
        }
        _oss << '\n';

        _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
        for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++ variableIdx) {
            _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx].typeStr(), _variableWidths[variableIdx]);
        }
        _oss << '\n';

        _oss << fmt::format("{:^{}}|", "", _maxVariableWidth);
        for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++ variableIdx)
            _oss << fmt::format("{:^{}}|", _variableInfos[variableIdx]._isBasic ? "BASIC" : "NON-BASIC", _variableWidths[variableIdx]);
        _oss << '\n';
    }

    template <typename T>
    void printMatrixWithRHS(const std::map<int, int>& rowToBasisColumnIdxMap,
                            const Matrix<T>& matrix, const std::vector<T>& rightHandSides)
    {
        for (int rowIdx = 0; rowIdx < _rowInfos.size(); ++rowIdx) {
            const auto columnBasicIdxIt = rowToBasisColumnIdxMap.find(rowIdx);
            if (columnBasicIdxIt == rowToBasisColumnIdxMap.end())
                _oss << fmt::format("{:^{}}|", "NOT FOUND", _maxVariableWidth);
            else
                _oss << fmt::format("{:^{}}|", _variableInfos[columnBasicIdxIt->second]._label, _maxVariableWidth);

            for (int variableIdx = 0; variableIdx < _variableInfos.size(); ++variableIdx)
                _oss << fmt::format("{:>{}}|", (' ' + std::to_string(matrix[rowIdx][variableIdx])), _variableWidths[variableIdx]);

            _oss << fmt::format("{:^{}}|", rowTypeToStr(_rowInfos[rowIdx]._type), CONSTRAINT_SIGN_WIDTH);
            _oss << fmt::format("{:>{}}|\n", (' ' + std::to_string(rightHandSides[rowIdx])), COEFFICIENT_WIDTH);
        }
    }


    template <typename T>
    void printInverseBasisWithDual(const Matrix<T>& basisMatrixInverse)
    {
        _oss << "INVERSE BASIS WITH DUAL\n";
        for (int rowIdx = 0; rowIdx < basisMatrixInverse.size(); ++rowIdx)
        {
            for (int variableIdx = 0; variableIdx < basisMatrixInverse.size(); ++variableIdx)
                _oss << fmt::format("{:>{}}|", (' ' + std::to_string(basisMatrixInverse[rowIdx][variableIdx])), COEFFICIENT_WIDTH);
            _oss << '\n';
        }
    }

    template <typename T>
    void printReducedCosts(const std::vector<T>& reducedCosts)
    {
        _oss << "REDUCED COSTS\n";
        for (int columnIdx = 0; columnIdx < reducedCosts.size(); ++columnIdx)
            _oss << fmt::format("{:>{}}|", (' ' + std::to_string(reducedCosts[columnIdx])), COEFFICIENT_WIDTH);
        _oss << '\n';
    }

    void printLineBreak()
    {
        _oss.width(_totalWidth);
        _oss.fill('-');
        _oss << '-' << '\n';
        _oss.fill(' ');
    }

    template <typename T>
    void printInLpSolveFormat(const Matrix<T>& matrix, const std::vector<T>& rightHandSides)
    {
        _oss << "min: ";
        for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx)
        {
            if (matrix[0][varIdx] != 0.0)
                _oss << fmt::format(" {:+f} {}", matrix[0][varIdx], _variableInfos[varIdx]._label);
        }
        _oss << ";\n";

        for (int rowIdx = 1; rowIdx < _rowInfos.size(); ++rowIdx)
        {
            for (int varIdx = 0; varIdx < _variableInfos.size(); ++varIdx)
            {
                if (matrix[rowIdx][varIdx] != 0.0)
                    _oss << fmt::format(" {:+f} {}", matrix[rowIdx][varIdx], _variableInfos[varIdx]._label);
            }
            _oss << ' ' << rowTypeToStr(_rowInfos[rowIdx]._type);
            _oss << fmt::format(" {:+f}", rightHandSides[rowIdx]) << ";\n";
        }
    }

    std::string toString() const
    {
        return _oss.str();
    }

    const std::vector<VariableInfo>& _variableInfos;
    const std::vector<RowInfo>& _rowInfos;
    std::vector<int> _variableWidths;
    std::ostringstream _oss;
    int _totalWidth{0};
    int _maxVariableWidth{0};
};

#endif //GMISOLVER_LPPRINTER_H
