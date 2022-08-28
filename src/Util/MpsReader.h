#ifndef GMISOLVER_MPSREADER_H
#define GMISOLVER_MPSREADER_H

#include "src/DataModel/EnumTypes.h"
#include "src/DataModel/LinearProgram.h"

#include <map>
#include <string>
#include <vector>
#include <optional>

struct MpsReader
{
public:
    template <typename T>
    static std::optional<LinearProgram<T>> read(const std::string& filePath);
};

#endif //GMISOLVER_MPSREADER_H
