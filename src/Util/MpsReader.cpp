#include "src/Util/MpsReader.h"

#include <fstream>
#include <spdlog/spdlog.h>
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
std::optional<LinearProgram<T>> MpsReader::read(const std::string &filePath) {
  std::ifstream fileStream(filePath);
  if (!fileStream.is_open()) {
    spdlog::warn("Could not open {} file", filePath);
    return std::nullopt;
  }

  LinearProgram<T> linearProgram;
  SectionType currentSection = SectionType::UNDEFINED;
  bool currentSectionIsInteger = false;
  std::string readLine;
  std::string currentWord;

  std::map<std::string, int> rowLabelToRowIdxMap;
  std::map<std::string, int> variableLabelToVariableIdxMap;

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
      if (!readSectionType.has_value())
        spdlog::warn("Unrecognized new section type {}", lineParts[0]);

      continue;
    }

    switch (currentSection) {
    case SectionType::NAME: {
      linearProgram._name = lineParts[1];
      break;
    }
    case SectionType::ROWS: {
      if (lineParts.size() != 2) {
        spdlog::warn("Line {} specifying row should have exactly 2 parts",
                     readLine);
        for (const auto &linePart : lineParts)
          spdlog::info("LINE PART {}", linePart);
        break;
      }

      const auto &rowTypeStr = lineParts[0];
      const auto &rowLabelStr = lineParts[1];

      const auto readRowType = stringToRowType(rowTypeStr);
      if (!readRowType.has_value()) {
        spdlog::warn("Unrecognized row type {}", rowTypeStr);
        break;
      }

      const RowInfo newRowInfo{rowLabelStr, *readRowType};
      if (*readRowType == RowType::OBJECTIVE)
        linearProgram._objectiveInfo = newRowInfo;
      else
      {
        const auto [_, inserted] = rowLabelToRowIdxMap.try_emplace(
          rowLabelStr, linearProgram._rowInfos.size());
        if (!inserted)
          spdlog::warn("Duplicated row label {}", rowLabelStr);

        linearProgram._rowInfos.push_back(newRowInfo);
      }

      break;
    }
    case SectionType::COLUMNS: {
      if (lineParts.size() != 3 && lineParts.size() != 5) {
        spdlog::warn("Unexpected number of elements in column line {}",
                     readLine);
        break;
      }

      if (lineParts.size() == 3 && lineParts[1] == MARKER_KEYWORD) {
        const auto &integerSectionStr = lineParts[2];
        if (integerSectionStr == INTEGER_SECTION_START_KEYWORD)
          currentSectionIsInteger = true;
        else if (integerSectionStr == INTEGER_SECTION_END_KEYWORD)
          currentSectionIsInteger = false;
        else
          spdlog::warn("Unrecognized integer section keyword {}",
                       integerSectionStr);

        break;
      }

      const auto &variableLabelStr = lineParts[0];

      if (variableLabelStr.find(Constants::SLACK_SUFFIX) != std::string::npos ||
          variableLabelStr.find(Constants::ARTIFICIAL_SUFFIX) !=
              std::string::npos) {
        spdlog::warn("Disallowed variable label {}", variableLabelStr);
        break;
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
        linearProgram._objective.resize(linearProgram._variableInfos.size());
        linearProgram._constraintMatrix.resize(linearProgram._rowInfos.size());
        for (auto &coeffRow : linearProgram._constraintMatrix)
          coeffRow.resize(linearProgram._variableInfos.size());
      }

      const auto updateLpMatrix = [&](const auto &rowLabelStr,
                                      const auto &coefficientValueStr) {
        if (rowLabelStr == linearProgram._objectiveInfo._label)
        {
          linearProgram._objective[variableIdx] = convert(coefficientValueStr);
          return;
        }

        const auto foundRowIt = rowLabelToRowIdxMap.find(rowLabelStr);
        if (foundRowIt == rowLabelToRowIdxMap.end()) {
          spdlog::warn(
              "Row {} given in column section doesn't correspond to any row",
              rowLabelStr);
          return;
        }

        linearProgram._constraintMatrix[foundRowIt->second][variableIdx] =
            convert(coefficientValueStr);
      };

      updateLpMatrix(lineParts[1], lineParts[2]);
      if (lineParts.size() == 5)
        updateLpMatrix(lineParts[3], lineParts[4]);

      break;
    }
    case SectionType::RHS: {
      if (lineParts.size() != 3 && lineParts.size() != 5) {
        spdlog::warn("Unexpected number of elements in rhs line {}", readLine);
        break;
      }

      const auto updateLpRhs = [&](const auto &rowLabelStr,
                                   const auto &coefficientValueStr) {
        const auto foundRowIt = rowLabelToRowIdxMap.find(rowLabelStr);
        if (foundRowIt == rowLabelToRowIdxMap.end()) {
          spdlog::warn("Row label {} given in column section doesn't "
                       "correspond to any row",
                       rowLabelStr);
          return;
        }

        linearProgram._rightHandSides.resize(linearProgram._rowInfos.size());
        linearProgram._rightHandSides[foundRowIt->second] =
            convert(coefficientValueStr);
      };

      updateLpRhs(lineParts[1], lineParts[2]);
      if (lineParts.size() == 5)
        updateLpRhs(lineParts[3], lineParts[4]);

      break;
    }
    case SectionType::BOUNDS: {
      if (lineParts.size() != 4) {
        spdlog::warn("Unexpected number of elements in bounds line {}",
                     readLine);
        break;
      }

      const auto &boundTypeStr = lineParts[0];
      const auto readBoundType = stringToBoundType(boundTypeStr);
      if (!readBoundType.has_value()) {
        spdlog::warn("Unsupported bound type {}", boundTypeStr);
        break;
      }

      const auto &variableStr = lineParts[2];
      const auto foundVariableIt =
          variableLabelToVariableIdxMap.find(variableStr);
      if (foundVariableIt == variableLabelToVariableIdxMap.end()) {
        spdlog::warn("Variable label {} given in bounds section doesn't "
                     "correspond to any variable",
                     variableStr);
        break;
      }

      const auto variableIdx = foundVariableIt->second;
      const auto &coefficientValueStr = lineParts[3];

      linearProgram._constraintMatrix.emplace_back(
          linearProgram._variableInfos.size())[variableIdx] = 1;
      switch (*readBoundType) {
      // TODO: optimize handling bounds
      case BoundType::LOWER_BOUND: {
        linearProgram._rightHandSides.push_back(convert(coefficientValueStr));
        linearProgram._rowInfos.push_back(
            RowInfo{{}, RowType::GREATER_THAN_OR_EQUAL});
        break;
      }
      case BoundType::UPPER_BOUND: {
        linearProgram._rightHandSides.push_back(convert(coefficientValueStr));
        linearProgram._rowInfos.push_back(
            RowInfo{{}, RowType::LESS_THAN_OR_EQUAL});
        break;
      }
      case BoundType::BINARY_VARIABLE: {
        linearProgram._rightHandSides.push_back(1);
        linearProgram._rowInfos.push_back(
            RowInfo{{}, RowType::LESS_THAN_OR_EQUAL});
        linearProgram._variableInfos[variableIdx]._type = VariableType::INTEGER;
        break;
      }
      }

      break;
    }
    case SectionType::END:
      break;
    default: {
      spdlog::warn("Undefined section type");
      break;
    }
    }

    if (currentSection == SectionType::END)
      break;
  }

  if (linearProgram._constraintMatrix.empty() ||
      linearProgram._variableInfos.empty()) {
    spdlog::warn("No data in COLUMNS section");
    return std::nullopt;
  }

  if (linearProgram._rowInfos.empty()) {
    spdlog::warn("No data in ROWS section");
    return std::nullopt;
  }

  if (linearProgram._rowInfos.size() != linearProgram._rightHandSides.size()) {
    spdlog::warn(
        "Number of constraints doesn't match number of right hand sides");
    return std::nullopt;
  }

  if (linearProgram._objectiveInfo._type != RowType::OBJECTIVE) {
    spdlog::warn("Could not find objective row");
    return std::nullopt;
  }

  // TODO - add more integrity checks
  return linearProgram;
}

template std::optional<LinearProgram<double>>
MpsReader::read(const std::string &);