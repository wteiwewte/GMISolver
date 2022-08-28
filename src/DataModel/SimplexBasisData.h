#ifndef GMISOLVER_SIMPLEXBASISDATA_H
#define GMISOLVER_SIMPLEXBASISDATA_H


#include <iostream>
#include <map>
#include <spdlog/spdlog.h>

struct SimplexBasisData {
  std::map<int, int> _rowToBasisColumnIdxMap;
  std::vector<bool> _isBasicColumnIndexBitset;
};


#endif // GMISOLVER_SIMPLEXBASISDATA_H
