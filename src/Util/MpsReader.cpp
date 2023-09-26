#include "src/Util/MpsReader.h"

#include "src/Util/SpdlogHeader.h"

#include <fstream>
#include <sstream>
#include <string_view>

namespace {
constexpr std::string_view MARKER_KEYWORD = "'MARKER'";
constexpr std::string_view INTEGER_SECTION_START_KEYWORD = "'INTORG'";
constexpr std::string_view INTEGER_SECTION_END_KEYWORD = "'INTEND'";

std::vector<std::string> splitString(const std::string &str) {
  std::stringstream stream(str);

  std::vector<std::string> splittedStreamParts;
  while (stream >> splittedStreamParts.emplace_back())
    ;

  std::erase_if(splittedStreamParts, [](const std::string &part) {
    return std::all_of(part.begin(), part.end(), isspace);
  });

  return splittedStreamParts;
}

double convert(const std::string &str) { return std::stod(str); }

template <typename T>
int countBoundsSpecified(const std::vector<std::optional<T>> &bounds) {
  return std::count_if(
      bounds.begin(), bounds.end(),
      [](const std::optional<T> &bound) { return bound.has_value(); });
}
} // namespace

template <typename T>
std::optional<LinearProgram<T>>
MpsReader<T>::read(const std::string &filePath) {
  std::ifstream fileStream(filePath);
  if (!fileStream.is_open()) {
    SPDLOG_DEBUG("Could not open {} file", filePath);
    return std::nullopt;
  }

  LinearProgram<T> linearProgram;
  SectionType currentSection = SectionType::UNDEFINED;
  bool currentSectionIsInteger = false;
  std::string readLine;
  std::string currentWord;

  std::map<std::string, int> rowLabelToRowIdxMap;
  std::map<std::string, int> variableLabelToVariableIdxMap;

  const auto addNewRowInfo = [&](const auto &rowLabelStr,
                                 const auto rowType) -> std::optional<int> {
    const RowInfo newRowInfo{rowLabelStr, rowType};
    const auto [_, inserted] = rowLabelToRowIdxMap.try_emplace(
        rowLabelStr, linearProgram._rowInfos.size());
    if (!inserted) {
      SPDLOG_WARN("Duplicated row label {}", rowLabelStr);
      return false;
    }

    if (rowType == RowType::OBJECTIVE) {
      // Consecutive objective rows are discarded
      if (linearProgram._objectiveInfo._type == RowType::UNKNOWN) {
        linearProgram._objectiveInfo = newRowInfo;
      }
    } else {
      linearProgram._rowInfos.push_back(newRowInfo);
      linearProgram._constraintMatrix.resize(linearProgram._rowInfos.size());
      linearProgram._rightHandSides.resize(linearProgram._rowInfos.size());
    }
    return true;
  };

  const auto tryAddNewVar =
      [&](const auto &variableLabelStr) -> std::optional<int> {
    if (variableLabelStr.find(Constants::SLACK_SUFFIX) != std::string::npos ||
        variableLabelStr.find(Constants::ARTIFICIAL_SUFFIX) !=
            std::string::npos) {
      SPDLOG_WARN("Disallowed variable label {}", variableLabelStr);
      return std::nullopt;
    }

    const auto [variableIt, inserted] =
        variableLabelToVariableIdxMap.try_emplace(
            variableLabelStr, linearProgram._variableInfos.size());
    const auto variableIdx = variableIt->second;
    if (inserted) {
      linearProgram._variableInfos.push_back(
          VariableInfo{variableLabelStr, currentSectionIsInteger
                                             ? VariableType::INTEGER
                                             : VariableType::CONTINUOUS});
      linearProgram._variableLabelSet.insert(variableLabelStr);
      linearProgram._variableLowerBounds.resize(
          linearProgram._variableInfos.size());
      linearProgram._variableUpperBounds.resize(
          linearProgram._variableInfos.size());
      linearProgram._objective.resize(linearProgram._variableInfos.size());
      for (auto &coeffRow : linearProgram._constraintMatrix)
        coeffRow.resize(linearProgram._variableInfos.size());
    }
    return variableIdx;
  };

  while (std::getline(fileStream, readLine)) {
    if (readLine.empty() || readLine[0] == '*')
      continue;
    const auto lineParts = splitString(readLine);

    if (lineParts.empty())
      continue;

    const auto readSectionType = stringToSectionType(lineParts[0]);
    if (readSectionType.has_value())
      currentSection = *readSectionType;

    const bool newSectionLine = lineParts.size() == 1;
    if (newSectionLine) {
      if (!readSectionType.has_value()) {
        SPDLOG_WARN("Unrecognized new section type {}", lineParts[0]);
        return std::nullopt;
      }

      if (*readSectionType == SectionType::END) {
        break;
      }

      continue;
    }

    switch (currentSection) {
    case SectionType::NAME: {
      linearProgram._name = lineParts[1];
      break;
    }
    case SectionType::ROWS: {
      if (lineParts.size() != 2) {
        SPDLOG_WARN("Line {} specifying row should have exactly 2 parts",
                    readLine);
        for (const auto &linePart : lineParts)
          SPDLOG_INFO("LINE PART {}", linePart);
        return std::nullopt;
      }

      const auto &rowTypeStr = lineParts[0];
      const auto &rowLabelStr = lineParts[1];

      const auto readRowType = stringToRowType(rowTypeStr);
      if (!readRowType.has_value()) {
        SPDLOG_WARN("Unrecognized row type {}", rowTypeStr);
        return std::nullopt;
      }

      if (!addNewRowInfo(rowLabelStr, *readRowType))
        return std::nullopt;

      break;
    }
    case SectionType::COLUMNS: {
      if (lineParts.size() != 3 && lineParts.size() != 5) {
        SPDLOG_WARN("Unexpected number of elements in column line {}",
                    readLine);
        return std::nullopt;
      }

      if (lineParts.size() == 3 && lineParts[1] == MARKER_KEYWORD) {
        const auto &integerSectionStr = lineParts[2];
        if (integerSectionStr == INTEGER_SECTION_START_KEYWORD)
          currentSectionIsInteger = true;
        else if (integerSectionStr == INTEGER_SECTION_END_KEYWORD)
          currentSectionIsInteger = false;
        else {
          SPDLOG_WARN("Unrecognized integer section keyword {}",
                      integerSectionStr);
          return std::nullopt;
        }

        break;
      }

      const auto &variableLabelStr = lineParts[0];
      const auto variableIdx = tryAddNewVar(variableLabelStr);
      if (!variableIdx.has_value())
        return std::nullopt;

      const auto updateLpMatrix = [&](const auto &rowLabelStr,
                                      const auto &coefficientValueStr) {
        if (rowLabelStr == linearProgram._objectiveInfo._label) {
          linearProgram._objective[*variableIdx] = convert(coefficientValueStr);
          return true;
        }

        const auto foundRowIt = rowLabelToRowIdxMap.find(rowLabelStr);
        if (foundRowIt == rowLabelToRowIdxMap.end()) {
          SPDLOG_WARN(
              "Row {} given in column section doesn't correspond to any row",
              rowLabelStr);
          return false;
        }

        linearProgram._constraintMatrix[foundRowIt->second][*variableIdx] =
            convert(coefficientValueStr);
        return true;
      };

      if (!updateLpMatrix(lineParts[1], lineParts[2]))
        return std::nullopt;
      if (lineParts.size() == 5) {
        if (!updateLpMatrix(lineParts[3], lineParts[4]))
          return std::nullopt;
      }

      break;
    }
    case SectionType::RHS: {
      if (lineParts.size() != 3 && lineParts.size() != 5) {
        SPDLOG_WARN("Unexpected number of elements in rhs line {}", readLine);
        return std::nullopt;
      }

      const auto updateLpRhs = [&](const auto &rowLabelStr,
                                   const auto &coefficientValueStr) {
        const auto foundRowIt = rowLabelToRowIdxMap.find(rowLabelStr);
        if (foundRowIt == rowLabelToRowIdxMap.end()) {
          SPDLOG_WARN("Row label {} given in column section doesn't "
                      "correspond to any row",
                      rowLabelStr);
          return false;
        }

        linearProgram._rightHandSides.resize(linearProgram._rowInfos.size());
        linearProgram._rightHandSides[foundRowIt->second] =
            convert(coefficientValueStr);
        return true;
      };

      if (!updateLpRhs(lineParts[1], lineParts[2]))
        return std::nullopt;
      if (lineParts.size() == 5) {
        if (!updateLpRhs(lineParts[3], lineParts[4]))
          return std::nullopt;
      }

      break;
    }
    case SectionType::RANGES: {
      if (lineParts.size() != 3 && lineParts.size() != 5) {
        SPDLOG_WARN("Unexpected number of elements in ranges line {}",
                    readLine);
        return std::nullopt;
      }

      const auto addRangeConstraint = [&](const auto &rowLabelStr,
                                          const auto &coefficientValueStr) {
        const auto foundRowIt = rowLabelToRowIdxMap.find(rowLabelStr);
        if (foundRowIt == rowLabelToRowIdxMap.end()) {
          SPDLOG_WARN("Row label {} given in column section doesn't "
                      "correspond to any row",
                      rowLabelStr);
          return false;
        }
        const auto rowIdx = foundRowIt->second;
        const auto rowType = linearProgram._rowInfos[rowIdx]._type;

        std::optional<RowType> rangeConstraintType;
        std::optional<T> rangeRHS;
        switch (rowType) {
        case RowType::GREATER_THAN_OR_EQUAL: {
          rangeConstraintType = RowType::LESS_THAN_OR_EQUAL;
          rangeRHS = linearProgram._rightHandSides[rowIdx] +
                     std::abs(convert(coefficientValueStr));
          break;
        }
        case RowType::LESS_THAN_OR_EQUAL: {
          rangeConstraintType = RowType::GREATER_THAN_OR_EQUAL;
          rangeRHS = linearProgram._rightHandSides[rowIdx] -
                     std::abs(convert(coefficientValueStr));
          break;
        }
        default: {
          SPDLOG_WARN("Row type {} given in ranges section isn't supported",
                      rowTypeToStr(rowType));
          return false;
        }
        }

        if (!rangeConstraintType.has_value() || !rangeRHS.has_value())
          return false;

        if (!addNewRowInfo(rowLabelStr + "_RANGES", *rangeConstraintType))
          return false;

        linearProgram._constraintMatrix.back() =
            linearProgram._constraintMatrix[rowIdx];
        linearProgram._rightHandSides.back() = *rangeRHS;
        return true;
      };

      if (!addRangeConstraint(lineParts[1], lineParts[2]))
        return std::nullopt;
      if (lineParts.size() == 5) {
        if (!addRangeConstraint(lineParts[3], lineParts[4]))
          return std::nullopt;
      }

      break;
    }
    case SectionType::BOUNDS: {
      if (lineParts.size() != 3 && lineParts.size() != 4) {
        SPDLOG_WARN("Unexpected number of elements in bounds line {}",
                    readLine);
        return std::nullopt;
      }

      const auto &boundTypeStr = lineParts[0];
      const auto readBoundType = stringToBoundType(boundTypeStr);
      if (!readBoundType.has_value()) {
        SPDLOG_WARN("Unsupported bound type {}", boundTypeStr);
        return std::nullopt;
      }

      if (!allowedLinesCount(*readBoundType).contains(lineParts.size())) {
        SPDLOG_WARN("Unexpected number of elements in bound line {}", readLine);
        return std::nullopt;
      }

      const auto &variableLabelStr = lineParts[2];
      const auto foundVariableIt =
          variableLabelToVariableIdxMap.find(variableLabelStr);
      if (foundVariableIt == variableLabelToVariableIdxMap.end()) {
        SPDLOG_WARN("Variable label {} given in bounds section doesn't "
                    "correspond to any variable",
                    variableLabelStr);
        return std::nullopt;
      }

      const auto variableIdx = foundVariableIt->second;
      const auto &coefficientValueStr = lineParts[3];

      switch (*readBoundType) {
      // TODO: optimize handling bounds
      case BoundType::LOWER_BOUND: {
        linearProgram._variableLowerBounds[variableIdx] =
            convert(coefficientValueStr);
        break;
      }
      case BoundType::UPPER_BOUND: {
        if (!linearProgram._variableLowerBounds[variableIdx].has_value())
          linearProgram._variableLowerBounds[variableIdx] = 0.0;

        linearProgram._variableUpperBounds[variableIdx] =
            convert(coefficientValueStr);
        break;
      }
      case BoundType::UPPER_BOUND_INTEGER: {
        if (!linearProgram._variableLowerBounds[variableIdx].has_value())
          linearProgram._variableLowerBounds[variableIdx] = 0.0;

        linearProgram._variableUpperBounds[variableIdx] =
            convert(coefficientValueStr);
        linearProgram._variableInfos[variableIdx]._type = VariableType::INTEGER;
        break;
      }
      case BoundType::LOWER_BOUND_INTEGER: {
        linearProgram._variableLowerBounds[variableIdx] =
            convert(coefficientValueStr);
        linearProgram._variableInfos[variableIdx]._type = VariableType::INTEGER;
        break;
      }
      case BoundType::BINARY_VARIABLE: {
        linearProgram._variableLowerBounds[variableIdx] = 0.0;
        linearProgram._variableUpperBounds[variableIdx] = 1.0;
        linearProgram._variableInfos[variableIdx]._type = VariableType::INTEGER;
        break;
      }
      case BoundType::FIXED_VARIABLE: {
        // TODO - maybe opt it
        linearProgram._variableInfos[variableIdx]._isFixed = true;
        linearProgram._variableLowerBounds[variableIdx] =
            linearProgram._variableUpperBounds[variableIdx] =
                convert(coefficientValueStr);
        break;
      }
      case BoundType::FREE_VARIABLE: {
        const auto newMinusVariableIdx =
            tryAddNewVar(variableLabelStr + "_MINUS_VAR");
        if (!newMinusVariableIdx.has_value())
          return std::nullopt;

        linearProgram._objective[*newMinusVariableIdx] =
            -linearProgram._objective[variableIdx];
        for (auto &coeffRow : linearProgram._constraintMatrix)
          coeffRow[*newMinusVariableIdx] = -coeffRow[variableIdx];

        linearProgram._variableLowerBounds[variableIdx] = 0.0;
        linearProgram._variableLowerBounds[*newMinusVariableIdx] = 0.0;
        break;
      }
      case BoundType::LOWER_BOUND_MINUS_INF: {
        linearProgram._variableLowerBounds[variableIdx] = 0.0;
        linearProgram._objective[variableIdx] =
            -linearProgram._objective[variableIdx];
        for (auto &coeffRow : linearProgram._constraintMatrix)
          coeffRow[variableIdx] = -coeffRow[variableIdx];

        break;
      }
      }

      break;
    }
    case SectionType::END: {
      break;
    }
    default: {
      SPDLOG_WARN("Undefined section type");
      return std::nullopt;
    }
    }

    if (currentSection == SectionType::END)
      break;
  }

  if (linearProgram._constraintMatrix.empty() ||
      linearProgram._variableInfos.empty()) {
    SPDLOG_WARN("No data in COLUMNS section");
    return std::nullopt;
  }

  if (linearProgram._rowInfos.empty()) {
    SPDLOG_WARN("No data in ROWS section");
    return std::nullopt;
  }

  if (linearProgram._rowInfos.size() != linearProgram._rightHandSides.size()) {
    SPDLOG_WARN(
        "Number of constraints doesn't match number of right hand sides");
    return std::nullopt;
  }

  if (linearProgram._objectiveInfo._type != RowType::OBJECTIVE) {
    SPDLOG_WARN("Could not find objective row");
    return std::nullopt;
  }

  if (!finalizeBounds(linearProgram))
    return std::nullopt;

  // TODO - add more integrity checks
  return linearProgram;
}

template <typename T>
bool MpsReader<T>::finalizeBounds(LinearProgram<T> &linearProgram) {
  for (int varIdx = 0; varIdx < linearProgram._variableInfos.size(); ++varIdx)
    if (!linearProgram._variableLowerBounds[varIdx].has_value() &&
        !linearProgram._variableUpperBounds[varIdx].has_value())
      linearProgram._variableLowerBounds[varIdx] = 0.0;

  SPDLOG_INFO("PROGRAM NAME {}, VARIABLE COUNT {}, ROW COUNT {}, LOWER BOUNDS "
              "{}, UPPER BOUNDS {}",
              linearProgram._name, linearProgram._variableInfos.size(),
              linearProgram._rowInfos.size(),
              countBoundsSpecified(linearProgram._variableLowerBounds),
              countBoundsSpecified(linearProgram._variableUpperBounds));

  return true;
}

template struct MpsReader<double>;
template struct MpsReader<long double>;
