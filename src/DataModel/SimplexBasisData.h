#ifndef GMISOLVER_SIMPLEXBASISDATA_H
#define GMISOLVER_SIMPLEXBASISDATA_H

#include "src/Util/SpdlogHeader.h"

#include <vector>

#include <boost/dynamic_bitset.hpp>

struct SimplexBasisData {
  void resizeVarCount(const size_t newVarCount) {
    _isBasicColumnIndexBitset.resize(newVarCount);
    _isColumnAtLowerBoundBitset.resize(newVarCount);
    _isColumnAtUpperBoundBitset.resize(newVarCount);
  }

  void restoreMapping() {
    std::vector<int> unmappedRowIndices;
    std::set<int> mappedColumnIndices;
    for (int rowIdx = 0; rowIdx < _rowToBasisColumnIdxMap.size(); ++rowIdx) {
      if (!_isBasicColumnIndexBitset[_rowToBasisColumnIdxMap[rowIdx]]) {
        unmappedRowIndices.push_back(rowIdx);
      } else {
        mappedColumnIndices.insert(_rowToBasisColumnIdxMap[rowIdx]);
      }
    }
    SPDLOG_INFO(
        "UNMAPPED ROW COUNT {} MAPPED COLUMN COUNT {} ALL COLUMN COUNT {}",
        unmappedRowIndices.size(), mappedColumnIndices.size(),
        _isBasicColumnIndexBitset.size());
    if (unmappedRowIndices.empty())
      return;

    int currentUnmappedRowIdx = 0;
    for (int colIdx = 0; colIdx < _isBasicColumnIndexBitset.size(); ++colIdx) {
      if (_isBasicColumnIndexBitset[colIdx] &&
          !mappedColumnIndices.contains(colIdx)) {
        _rowToBasisColumnIdxMap[unmappedRowIndices[currentUnmappedRowIdx++]] =
            colIdx;
      }
    }
  }

  std::vector<int> _rowToBasisColumnIdxMap;
  boost::dynamic_bitset<> _isBasicColumnIndexBitset;
  boost::dynamic_bitset<> _isColumnAtLowerBoundBitset;
  boost::dynamic_bitset<> _isColumnAtUpperBoundBitset;
};

#endif // GMISOLVER_SIMPLEXBASISDATA_H
