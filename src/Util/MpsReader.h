#ifndef GMISOLVER_MPSREADER_H
#define GMISOLVER_MPSREADER_H

#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/LinearProgram.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

struct MpsReader {
public:
  template <typename T>
  static std::optional<LinearProgram<T>> read(const std::string &filePath);
};

#endif // GMISOLVER_MPSREADER_H
