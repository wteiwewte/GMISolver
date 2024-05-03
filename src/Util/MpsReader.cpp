#include "src/Util/MpsReader.h"

#include "src/Util/CommonFunctions.h"
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

  const auto objectiveVarIdx =
      setupObjectiveRowVar(currentSectionIsInteger, linearProgram);
  if (!objectiveVarIdx.has_value())
    return std::nullopt;

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
      if (!handleRowsSection(readLine, lineParts, linearProgram))
        return std::nullopt;

      break;
    }
    case SectionType::COLUMNS: {
      if (!handleColumnsSection(readLine, lineParts, linearProgram,
                                currentSectionIsInteger))
        return std::nullopt;

      break;
    }
    case SectionType::RHS: {
      if (!handleRHSSection(readLine, lineParts, linearProgram))
        return std::nullopt;

      break;
    }
    case SectionType::RANGES: {
      if (!handleRangesSection(readLine, lineParts, linearProgram))
        return std::nullopt;

      break;
    }
    case SectionType::BOUNDS: {
      if (!handleBoundsSection(readLine, lineParts, linearProgram))
        return std::nullopt;

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
  linearProgram.convertToStandardForm();
  setBoundsForNonFreeVars(linearProgram);
  roundBoundsForIntegerVars(linearProgram);

  if (linearProgram.isPureIP(false) &&
      linearProgram.allCoefficientsAreIntegers()) {
    linearProgram._variableInfos[*objectiveVarIdx]._type =
        VariableType::INTEGER;
  }
  linearProgram.logGeneralInformation();

  return linearProgram;
}

template <typename T>
std::optional<int>
MpsReader<T>::setupObjectiveRowVar(const bool currentSectionIsInteger,
                                   LinearProgram<T> &linearProgram) {
  if (!addNewRowInfo(std::string{"obj_row"}, RowType::EQUALITY, linearProgram))
    return std::nullopt;
  const auto objectiveVarIdx = tryAddNewVar(
      std::string{"obj_value"}, currentSectionIsInteger, linearProgram);
  if (!objectiveVarIdx.has_value())
    return std::nullopt;

  linearProgram
      ._constraintMatrix[OBJECTIVE_ROW_CONSTRAINT_IDX][*objectiveVarIdx] = 1;
  linearProgram._variableInfos[*objectiveVarIdx]._isObjectiveVar = true;
  linearProgram._variableInfos[*objectiveVarIdx]._isFree = true;
  return objectiveVarIdx;
}

template <typename T>
bool MpsReader<T>::handleRowsSection(const std::string &readLine,
                                     const std::vector<std::string> &lineParts,
                                     LinearProgram<T> &linearProgram) {
  if (lineParts.size() != 2) {
    SPDLOG_WARN("Line {} specifying row should have exactly 2 parts", readLine);
    for (const auto &linePart : lineParts)
      SPDLOG_INFO("LINE PART {}", linePart);
    return false;
  }

  const auto &rowTypeStr = lineParts[0];
  const auto &rowLabelStr = lineParts[1];

  const auto readRowType = stringToRowType(rowTypeStr);
  if (!readRowType.has_value()) {
    SPDLOG_WARN("Unrecognized row type {}", rowTypeStr);
    return false;
  }

  return addNewRowInfo(rowLabelStr, *readRowType, linearProgram);
}
template <typename T>
bool MpsReader<T>::handleColumnsSection(
    const std::string &readLine, const std::vector<std::string> &lineParts,
    LinearProgram<T> &linearProgram, bool &currentSectionIsInteger) {
  if (lineParts.size() != 3 && lineParts.size() != 5) {
    SPDLOG_WARN("Unexpected number of elements in column line {}", readLine);
    return false;
  }

  if (lineParts.size() == 3 && lineParts[1] == MARKER_KEYWORD) {
    const auto &integerSectionStr = lineParts[2];
    if (integerSectionStr == INTEGER_SECTION_START_KEYWORD)
      currentSectionIsInteger = true;
    else if (integerSectionStr == INTEGER_SECTION_END_KEYWORD)
      currentSectionIsInteger = false;
    else {
      SPDLOG_WARN("Unrecognized integer section keyword {}", integerSectionStr);
      return false;
    }

    return true;
  }

  const auto &variableLabelStr = lineParts[0];
  const auto variableIdx =
      tryAddNewVar(variableLabelStr, currentSectionIsInteger, linearProgram);
  if (!variableIdx.has_value())
    return false;

  if (!updateLpMatrix(lineParts[1], lineParts[2], *variableIdx, linearProgram))
    return false;
  if (lineParts.size() == 5) {
    if (!updateLpMatrix(lineParts[3], lineParts[4], *variableIdx,
                        linearProgram))
      return false;
  }
  return true;
}
template <typename T>
bool MpsReader<T>::handleRHSSection(const std::string &readLine,
                                    const std::vector<std::string> &lineParts,
                                    LinearProgram<T> &linearProgram) {
  if (lineParts.size() != 3 && lineParts.size() != 5) {
    SPDLOG_WARN("Unexpected number of elements in rhs line {}", readLine);
    return false;
  }

  const auto updateLpRhs = [&](const auto &rowLabelStr,
                               const auto &coefficientValueStr) {
    const auto foundRowIt = _rowLabelToRowIdxMap.find(rowLabelStr);
    if (foundRowIt == _rowLabelToRowIdxMap.end()) {
      SPDLOG_WARN("Row label {} given in rhs section doesn't "
                  "correspond to any row",
                  rowLabelStr);
      return false;
    }

    const auto rowIdx = foundRowIt->second;
    if (!rowIdx.has_value()) {
      SPDLOG_WARN("Discarded non-constraining row {} given in rhs section",
                  rowLabelStr);
      return false;
    }

    linearProgram._rightHandSides.resize(linearProgram._rowInfos.size());
    linearProgram._rightHandSides[*rowIdx] = convert(coefficientValueStr);
    return true;
  };

  if (!updateLpRhs(lineParts[1], lineParts[2]))
    return false;
  if (lineParts.size() == 5) {
    if (!updateLpRhs(lineParts[3], lineParts[4]))
      return false;
  }
  return true;
}
template <typename T>
bool MpsReader<T>::handleRangesSection(
    const std::string &readLine, const std::vector<std::string> &lineParts,
    LinearProgram<T> &linearProgram) {
  if (lineParts.size() != 3 && lineParts.size() != 5) {
    SPDLOG_WARN("Unexpected number of elements in ranges line {}", readLine);
    return false;
  }

  const auto addRangeConstraint = [&](const auto &rowLabelStr,
                                      const auto &coefficientValueStr) {
    const auto foundRowIt = _rowLabelToRowIdxMap.find(rowLabelStr);
    if (foundRowIt == _rowLabelToRowIdxMap.end()) {
      SPDLOG_WARN("Row label {} given in range section doesn't "
                  "correspond to any row",
                  rowLabelStr);
      return false;
    }
    const auto rowIdx = foundRowIt->second;
    if (!rowIdx.has_value()) {
      SPDLOG_WARN("Discarded non-constraining row {} given in range section",
                  rowLabelStr);
      return false;
    }
    const auto rowType = linearProgram._rowInfos[*rowIdx]._type;

    std::optional<RowType> rangeConstraintType;
    std::optional<T> rangeRHS;
    switch (rowType) {
    case RowType::GREATER_THAN_OR_EQUAL: {
      rangeConstraintType = RowType::LESS_THAN_OR_EQUAL;
      rangeRHS = linearProgram._rightHandSides[*rowIdx] +
                 std::abs(convert(coefficientValueStr));
      break;
    }
    case RowType::LESS_THAN_OR_EQUAL: {
      rangeConstraintType = RowType::GREATER_THAN_OR_EQUAL;
      rangeRHS = linearProgram._rightHandSides[*rowIdx] -
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

    if (!addNewRowInfo(rowLabelStr + "_RANGES", *rangeConstraintType,
                       linearProgram))
      return false;

    linearProgram._constraintMatrix.back() =
        linearProgram._constraintMatrix[*rowIdx];
    linearProgram._rightHandSides.back() = *rangeRHS;
    return true;
  };

  if (!addRangeConstraint(lineParts[1], lineParts[2]))
    return false;
  if (lineParts.size() == 5) {
    if (!addRangeConstraint(lineParts[3], lineParts[4]))
      return false;
  }
  return true;
}
template <typename T>
bool MpsReader<T>::handleBoundsSection(
    const std::string &readLine, const std::vector<std::string> &lineParts,
    LinearProgram<T> &linearProgram) {
  if (lineParts.size() != 3 && lineParts.size() != 4) {
    SPDLOG_WARN("Unexpected number of elements in bounds line {}", readLine);
    return false;
  }

  const auto &boundTypeStr = lineParts[0];
  const auto readBoundType = stringToBoundType(boundTypeStr);
  if (!readBoundType.has_value()) {
    SPDLOG_WARN("Unsupported bound type {}", boundTypeStr);
    return false;
  }

  if (!allowedLinesCount(*readBoundType).contains(lineParts.size())) {
    SPDLOG_WARN(
        "Unexpected number of elements ({}) in bound (of type {}) line {}",
        lineParts.size(), boundTypeStr, readLine);
    return false;
  }

  const auto &variableLabelStr = lineParts[2];
  const auto foundVariableIt =
      _variableLabelToVariableIdxMap.find(variableLabelStr);
  if (foundVariableIt == _variableLabelToVariableIdxMap.end()) {
    SPDLOG_WARN("Variable label {} given in bounds section doesn't "
                "correspond to any variable",
                variableLabelStr);
    return false;
  }

  const auto variableIdx = foundVariableIt->second;

  switch (*readBoundType) {
  case BoundType::LOWER_BOUND: {
    if (lineParts.size() != 4) {
      SPDLOG_WARN("Unexpected number of elements in LOWER_BOUND line {}",
                  readLine);
      return false;
    }
    const auto &coefficientValueStr = lineParts[3];
    linearProgram._variableLowerBounds[variableIdx] =
        convert(coefficientValueStr);
    break;
  }
  case BoundType::UPPER_BOUND: {
    if (lineParts.size() != 4) {
      SPDLOG_WARN("Unexpected number of elements in UPPER_BOUND line {}",
                  readLine);
      return false;
    }
    if (!linearProgram._variableLowerBounds[variableIdx].has_value()) {
      linearProgram._variableLowerBounds[variableIdx] = 0.0;
    }

    const auto &coefficientValueStr = lineParts[3];
    linearProgram._variableUpperBounds[variableIdx] =
        convert(coefficientValueStr);
    break;
  }
  case BoundType::UPPER_BOUND_INTEGER: {
    if (lineParts.size() != 4) {
      SPDLOG_WARN(
          "Unexpected number of elements in UPPER_BOUND_INTEGER line {}",
          readLine);
      return false;
    }
    if (!linearProgram._variableLowerBounds[variableIdx].has_value()) {
      linearProgram._variableLowerBounds[variableIdx] = 0.0;
    }

    const auto &coefficientValueStr = lineParts[3];
    linearProgram._variableUpperBounds[variableIdx] =
        convert(coefficientValueStr);
    linearProgram._variableInfos[variableIdx]._type = VariableType::INTEGER;
    break;
  }
  case BoundType::LOWER_BOUND_INTEGER: {
    if (lineParts.size() != 4) {
      SPDLOG_WARN(
          "Unexpected number of elements in LOWER_BOUND_INTEGER line {}",
          readLine);
      return false;
    }
    const auto &coefficientValueStr = lineParts[3];
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
    if (lineParts.size() != 4) {
      SPDLOG_WARN("Unexpected number of elements in FIXED_VARIABLE line {}",
                  readLine);
      return false;
    }
    linearProgram._variableInfos[variableIdx]._isFixed = true;
    const auto &coefficientValueStr = lineParts[3];
    linearProgram._variableLowerBounds[variableIdx] =
        linearProgram._variableUpperBounds[variableIdx] =
            convert(coefficientValueStr);
    break;
  }
  case BoundType::FREE_VARIABLE: {
    linearProgram._variableInfos[variableIdx]._isFree = true;
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
  case BoundType::UPPER_BOUND_PLUS_INF: {
    linearProgram._variableLowerBounds[variableIdx] = 0.0;
    break;
  }
  }
  return true;
}

template <typename T>
bool MpsReader<T>::addNewRowInfo(const std::string &rowLabelStr,
                                 const RowType rowType,
                                 LinearProgram<T> &linearProgram) {
  const RowInfo newRowInfo{rowLabelStr, rowType};
  const bool isObjectiveAlreadySet =
      linearProgram._objectiveInfo._type != RowType::UNKNOWN;
  const bool isRowNeeded =
      (rowType != RowType::OBJECTIVE) || !isObjectiveAlreadySet;
  const auto [_, inserted] = _rowLabelToRowIdxMap.try_emplace(
      rowLabelStr,
      (isRowNeeded ? std::optional<int>(linearProgram._rowInfos.size())
                   : std::nullopt));
  if (!inserted) {
    SPDLOG_WARN("Duplicated row label {}", rowLabelStr);
    return false;
  }

  if (isRowNeeded) {
    if (rowType == RowType::OBJECTIVE) {
      // First objective row is taken into account, consecutive ones are
      // discarded
      linearProgram._objectiveInfo = newRowInfo;
    } else {
      linearProgram._rowInfos.push_back(newRowInfo);
      linearProgram._constraintMatrix.resize(linearProgram._rowInfos.size());
      linearProgram._rightHandSides.resize(linearProgram._rowInfos.size());
    }
  }
  return true;
}
template <typename T>
std::optional<int>
MpsReader<T>::tryAddNewVar(const std::string &variableLabelStr,
                           const bool currentSectionIsInteger,
                           LinearProgram<T> &linearProgram) {
  if (variableLabelStr.find(Constants::SLACK_SUFFIX) != std::string::npos ||
      variableLabelStr.find(Constants::ARTIFICIAL_SUFFIX) !=
          std::string::npos) {
    SPDLOG_WARN("Disallowed variable label {}", variableLabelStr);
    return std::nullopt;
  }

  const auto [variableIt, inserted] =
      _variableLabelToVariableIdxMap.try_emplace(
          variableLabelStr, linearProgram._variableInfos.size());
  const auto variableIdx = variableIt->second;
  if (inserted) {
    linearProgram._variableInfos.push_back(VariableInfo{
        variableLabelStr, currentSectionIsInteger ? VariableType::INTEGER
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
}

template <typename T>
bool MpsReader<T>::updateLpMatrix(const std::string &rowLabelStr,
                                  const std::string &coefficientValueStr,
                                  const int variableIdx,
                                  LinearProgram<T> &linearProgram) {
  if (rowLabelStr == linearProgram._objectiveInfo._label) {
    linearProgram._objective[variableIdx] = convert(coefficientValueStr);
    linearProgram._constraintMatrix[OBJECTIVE_ROW_CONSTRAINT_IDX][variableIdx] =
        -convert(coefficientValueStr);
    return true;
  }

  const auto foundRowIt = _rowLabelToRowIdxMap.find(rowLabelStr);
  if (foundRowIt == _rowLabelToRowIdxMap.end()) {
    SPDLOG_WARN("Row {} given in column section doesn't correspond to any row",
                rowLabelStr);
    return false;
  }

  const auto rowIdx = foundRowIt->second;
  if (!rowIdx.has_value()) {
    return true;
  }

  linearProgram._constraintMatrix[*rowIdx][variableIdx] =
      convert(coefficientValueStr);
  return true;
};

template <typename T>
void MpsReader<T>::setBoundsForNonFreeVars(LinearProgram<T> &linearProgram) {
  for (int varIdx = 0; varIdx < linearProgram._variableInfos.size(); ++varIdx) {
    if (!linearProgram._variableInfos[varIdx]._isFree &&
        !linearProgram._variableLowerBounds[varIdx].has_value() &&
        !linearProgram._variableUpperBounds[varIdx].has_value()) {
      linearProgram._variableLowerBounds[varIdx] = 0.0;

      if (linearProgram._variableInfos[varIdx]._type == VariableType::INTEGER) {
        linearProgram._variableUpperBounds[varIdx] = 1.0;
      }
    }
  }
}
template <typename T>
void MpsReader<T>::roundBoundsForIntegerVars(LinearProgram<T> &linearProgram) {
  for (int varIdx = 0; varIdx < linearProgram._variableInfos.size(); ++varIdx) {
    if (linearProgram._variableInfos[varIdx]._type == VariableType::INTEGER) {
      if (linearProgram._variableLowerBounds[varIdx].has_value()) {
        const auto varLowerBound = *linearProgram._variableLowerBounds[varIdx];
        if (!isInteger(varLowerBound)) {
          SPDLOG_WARN(
              "INTEGER VAR IDX {} HAS NON-INTEGER LOWER BOUND {}, ROUNDING UP",
              varIdx, varLowerBound);
        }
        linearProgram._variableLowerBounds[varIdx] = std::ceil(varLowerBound);
      }

      if (linearProgram._variableUpperBounds[varIdx].has_value()) {
        const auto varUpperBound = *linearProgram._variableUpperBounds[varIdx];
        if (!isInteger(varUpperBound)) {
          SPDLOG_WARN("INTEGER VAR IDX {} HAS NON-INTEGER UPPER BOUND {}, "
                      "ROUNDING DOWN",
                      varIdx, varUpperBound);
        }
        linearProgram._variableUpperBounds[varIdx] = std::floor(varUpperBound);
      }
    }
  }
}

template struct MpsReader<double>;
template struct MpsReader<long double>;
