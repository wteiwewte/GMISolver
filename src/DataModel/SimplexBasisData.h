#ifndef GMISOLVER_SIMPLEXBASISDATA_H
#define GMISOLVER_SIMPLEXBASISDATA_H

#include <vector>

#include <boost/dynamic_bitset.hpp>

struct SimplexBasisData {
  std::vector<int> _rowToBasisColumnIdxMap;
  boost::dynamic_bitset<> _isBasicColumnIndexBitset;
  boost::dynamic_bitset<> _isColumnAtLowerBoundBitset;
  boost::dynamic_bitset<> _isColumnAtUpperBoundBitset;
};

#endif // GMISOLVER_SIMPLEXBASISDATA_H
