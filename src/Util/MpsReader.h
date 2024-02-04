#ifndef GMISOLVER_MPSREADER_H
#define GMISOLVER_MPSREADER_H

#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/LinearProgram.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

template <typename T> struct MpsReader {
public:
  static std::optional<LinearProgram<T>> read(const std::string &filePath);

private:
  static void setBoundsForNonFreeVars(LinearProgram<T> &linearProgram);
  static void roundBoundsForIntegerVars(LinearProgram<T> &linearProgram);
};

#endif // GMISOLVER_MPSREADER_H
