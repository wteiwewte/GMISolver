#ifndef GMISOLVER_SIMPLEXBASISDATA_H
#define GMISOLVER_SIMPLEXBASISDATA_H

#include <iostream>
#include <map>
#include <spdlog/spdlog.h>

struct SimplexBasisData {
  std::vector<int> _rowToBasisColumnIdxMap;
  std::vector<bool> _isBasicColumnIndexBitset;
  std::vector<bool> _isColumnAtLowerBoundBitset;
  std::vector<bool> _isColumnAtUpperBoundBitset;
};

#endif // GMISOLVER_SIMPLEXBASISDATA_H
