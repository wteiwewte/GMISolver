#include "src/Util/MatchingEdgeFormatReader.h"

#include "src/Util/CommonFunctions.h"
#include "src/Util/SpdlogHeader.h"

#include <fmt/format.h>
#include <fstream>
#include <sstream>
#include <string>

namespace {
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
} // namespace

template <typename T>
std::optional<LinearProgram<T>>
MatchingEdgeFormatReader<T>::read(const std::string &filePath) {
  std::ifstream fileStream(filePath);
  if (!fileStream.is_open()) {
    SPDLOG_DEBUG("Could not open {} file", filePath);
    return std::nullopt;
  }

  LinearProgram<T> linearProgram;
  std::string readLine;
  int numberOfVertices = 0;
  int numberOfEdges = 0;
  int currentEdgeIndex = 0;

  while (std::getline(fileStream, readLine)) {
    if (readLine.empty())
      continue;

    const auto lineParts = splitString(readLine);
    if (readLine[0] == 'p') {
      numberOfVertices = std::stoi(lineParts[2]);
      numberOfEdges = std::stoi(lineParts[3]);
      linearProgram._name =
          fmt::format("MAX_MATCHING_{}_{}", numberOfVertices, numberOfEdges);
      linearProgram._variableInfos.resize(numberOfEdges);
      linearProgram._constraintMatrix.resize(numberOfVertices);
      for (auto &row : linearProgram._constraintMatrix)
        row.resize(numberOfEdges, 0.0);

      linearProgram._variableLowerBounds.resize(numberOfEdges, 0.0);
      linearProgram._variableUpperBounds.resize(numberOfEdges, 1.0);
      linearProgram._objective.resize(numberOfEdges);
      linearProgram._rightHandSides.resize(numberOfVertices, 1.0);
      linearProgram._rowInfos.resize(
          numberOfVertices, RowInfo{._type = RowType::LESS_THAN_OR_EQUAL});

    } else if (lineParts.size() == 4) {
      int v = std::stoi(lineParts[1]);
      int w = std::stoi(lineParts[2]);
      int cost = std::stoi(lineParts[3]);
      linearProgram._variableInfos[currentEdgeIndex] = VariableInfo{
          ._label = fmt::format("{}-{}", v, w), ._type = VariableType::INTEGER};

      linearProgram._constraintMatrix[v - 1][currentEdgeIndex] = 1;
      linearProgram._constraintMatrix[w - 1][currentEdgeIndex] = 1;
      linearProgram._objective[currentEdgeIndex] = -cost;

      ++currentEdgeIndex;
    }
  }

  linearProgram.convertAllConstraintsToEquality();
  linearProgram.logGeneralInformation();
  SPDLOG_DEBUG(linearProgram.toString());

  return linearProgram;
}

template struct MatchingEdgeFormatReader<double>;
template struct MatchingEdgeFormatReader<long double>;