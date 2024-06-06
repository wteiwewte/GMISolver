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

  std::vector<int>
  fixMappingAfterRemoval(const std::vector<bool> &shouldVarBeRemoved,
                         const std::vector<bool> &shouldRowBeRemoved) {
    std::vector<int> oldRowIdxToNewRowIdx(shouldRowBeRemoved.size());
    int currentNewRowIdx = 0;
    for (int oldRowIdx = 0; oldRowIdx < shouldRowBeRemoved.size();
         ++oldRowIdx) {
      if (!shouldRowBeRemoved[oldRowIdx])
        oldRowIdxToNewRowIdx[oldRowIdx] = currentNewRowIdx++;
    }
    std::vector<int> oldVarIdxToNewVarIdx(shouldVarBeRemoved.size());
    int currentNewVarIdx = 0;
    for (int oldVarIdx = 0; oldVarIdx < shouldVarBeRemoved.size();
         ++oldVarIdx) {
      if (!shouldVarBeRemoved[oldVarIdx])
        oldVarIdxToNewVarIdx[oldVarIdx] = currentNewVarIdx++;
    }

    std::vector<int> newVarIndicesToBeRemapped;
    for (int rowIdx = 0; rowIdx < shouldRowBeRemoved.size(); ++rowIdx) {
      const auto basicVarIdx = _rowToBasisColumnIdxMap[rowIdx];
      if (shouldRowBeRemoved[rowIdx] && !shouldVarBeRemoved[basicVarIdx])
        newVarIndicesToBeRemapped.push_back(oldVarIdxToNewVarIdx[basicVarIdx]);
    }
    std::vector<int> newRowToNewBasisColumnIdxMap(currentNewRowIdx);
    int currentVarIndexToBeRemapped = 0;
    for (int oldRowIdx = 0; oldRowIdx < shouldRowBeRemoved.size();
         ++oldRowIdx) {
      if (!shouldRowBeRemoved[oldRowIdx]) {
        const int newRowIdx = oldRowIdxToNewRowIdx[oldRowIdx];
        const int oldBasicVarIdx = _rowToBasisColumnIdxMap[oldRowIdx];
        if (!shouldVarBeRemoved[oldBasicVarIdx]) {
          newRowToNewBasisColumnIdxMap[newRowIdx] =
              oldVarIdxToNewVarIdx[oldBasicVarIdx];
        } else {
          newRowToNewBasisColumnIdxMap[newRowIdx] =
              newVarIndicesToBeRemapped[currentVarIndexToBeRemapped++];
        }
      }
    }
    return newRowToNewBasisColumnIdxMap;
  }

  std::vector<int> _rowToBasisColumnIdxMap;
  boost::dynamic_bitset<> _isBasicColumnIndexBitset;
  boost::dynamic_bitset<> _isColumnAtLowerBoundBitset;
  boost::dynamic_bitset<> _isColumnAtUpperBoundBitset;
  std::optional<boost::dynamic_bitset<>> _isColumnEligibleBitset;
};

#endif // GMISOLVER_SIMPLEXBASISDATA_H
