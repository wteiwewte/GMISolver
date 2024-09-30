#ifndef GMISOLVER_MATCHINGEDGEFORMATREADER_H
#define GMISOLVER_MATCHINGEDGEFORMATREADER_H

#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/LinearProgram.h"

#include <optional>

template <typename T> struct MatchingEdgeFormatReader {
  std::optional<LinearProgram<T>> read(const std::string &filePath);
};

#endif // GMISOLVER_MATCHINGEDGEFORMATREADER_H
