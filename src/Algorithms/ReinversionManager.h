#ifndef GMISOLVER_REINVERSIONMANAGER_H
#define GMISOLVER_REINVERSIONMANAGER_H

#include "src/Util/SimplexTraits.h"

template <typename T, typename SimplexTraitsT> class SimplexTableau;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class ReinversionManager {
public:
  ReinversionManager(SimplexTableau<T, SimplexTraitsT> &simplexTableau,
                     const int32_t reinversionFrequency);

  bool tryReinverse();
  bool reinverse();

  int reinversionFrequency() const { return _reinversionFrequency; }

private:
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;

  bool reinverseBasis();
  bool reinverseBasisExplicit();
  bool reinverseBasisPFI();
  bool reinverseBasisPFISparse();

  void updateTableau();

  template <typename VectorType>
  std::optional<int>
  findPivotColumnFirstEligible(const std::vector<VectorType> &basisColumns,
                               const std::vector<bool> &isUnusedColumn,
                               const int rowIdx) const;
  template <typename VectorType>
  std::optional<int>
  findPivotColumnMaxAbsValue(const std::vector<VectorType> &basisColumns,
                             const std::vector<bool> &isUnusedColumn,
                             const int rowIdx) const;

  SimplexTableau<T, SimplexTraitsT> &_simplexTableau;
  const int32_t _reinversionFrequency;
  int32_t _iterCount{1};
};

#endif // GMISOLVER_REINVERSIONMANAGER_H
