#ifndef GMISOLVER_SIMPLEXTABLEAU_H
#define GMISOLVER_SIMPLEXTABLEAU_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/LinearProgram.h"
#include "src/DataModel/MatrixTypes.h"
#include "src/DataModel/SimplexBasisData.h"
#include "src/Util/SimplexTraits.h"

#include <set>

#include <boost/dynamic_bitset.hpp>

template <typename T, typename SimplexTraitsT> class DualSimplexGomory;
template <typename T, typename SimplexTraitsT> class PrimalSimplexGomory;
template <typename T, typename SimplexTraitsT> class LexicographicOptimizer;
template <typename T, typename SimplexTraitsT> class PrimalSimplex;
template <typename T, typename SimplexTraitsT> class PrimalSimplexPhaseOne;
template <typename T, typename SimplexTraitsT> class DualSimplex;
template <typename T, typename SimplexTraitsT> class DualSimplexPhaseOne;
template <typename T, typename SimplexTraitsT> class PhaseOneUtilities;
template <typename T, typename SimplexTraitsT> class ReinversionManager;
template <typename T, typename SimplexTraitsT> class SimplexTableauResizer;
template <typename T, typename SimplexTraitsT> class SimplexValidator;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram,
                 const SimplexType simplexType,
                 const SimplexTableauType simplexTableauType);

  std::vector<T> originalObjective() const;
  void calculateCurrentObjectiveValue(const SimplexPhase simplexPhase);
  void calculateSolution();

  std::string toString() const;
  std::string toStringObjectiveValue() const;
  std::string toStringSolution() const;
  std::string
  toStringSolutionWithDual(const LinearProgram<T> &dualProgram) const;
  std::string toStringLpSolveFormat() const;

  const std::string &getName() const { return _initialProgram.getName(); }

private:
  friend class DualSimplexGomory<T, SimplexTraitsT>;
  friend class PrimalSimplexGomory<T, SimplexTraitsT>;
  friend class LexicographicOptimizer<T, SimplexTraitsT>;
  friend class PrimalSimplex<T, SimplexTraitsT>;
  friend class PrimalSimplexPhaseOne<T, SimplexTraitsT>;
  friend class DualSimplex<T, SimplexTraitsT>;
  friend class DualSimplexPhaseOne<T, SimplexTraitsT>;
  friend class PhaseOneUtilities<T, SimplexTraitsT>;
  friend class ReinversionManager<T, SimplexTraitsT>;
  friend class SimplexTableauResizer<T, SimplexTraitsT>;
  friend class SimplexValidator<T, SimplexTraitsT>;

  using ElementaryMatrixT = typename SimplexTraitsT::ElementaryMatrixT;
  using VectorT = typename SimplexTraitsT::VectorT;
  using NumericalTraitsT = typename SimplexTraitsT::NumericalTraitsT;
  using SimpleAdderSafe = typename NumericalTraitsT::template SimpleAdder<
      typename NumericalTraitsT::SafeNumericalAddOp>;
  using PositiveNegativeAdderSafe =
      typename NumericalTraitsT::template PositiveNegativeAdder<
          typename NumericalTraitsT::SafeNumericalAddOp>;
  using KahanAdderSafe = typename NumericalTraitsT::template KahanAdder<
      typename NumericalTraitsT::SafeNumericalAddOp>;
  using KahanAdderNormal = typename NumericalTraitsT::template KahanAdder<
      typename NumericalTraitsT::NormalAddOp>;

  int basicColumnIdx(const int rowIdx) const {
    return _simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
  }
  bool isColumnEligibleToEnterBasis(const int colIdx);
  std::optional<T> curSatisfiedBound(const int varIdx);

  VectorT computeTableauColumnGeneric(const int colIdx);
  std::vector<T> retrieveTableauColumn(const int colIdx);
  std::vector<T> computeTableauColumnExplicit(const int colIdx);
  std::vector<T> computeTableauColumnPFI(const int colIdx);
  SparseVector<T> computeTableauColumnPFISparse(const int colIdx);

  std::shared_ptr<VectorT> computeTableauRowGeneric(const int rowIdx);
  std::shared_ptr<std::vector<T>> retrieveTableauRow(const int rowIdx);
  std::shared_ptr<std::vector<T>> computeTableauRowExplicit(const int rowIdx);
  std::shared_ptr<std::vector<T>> computeTableauRowPFI(const int rowIdx);
  std::shared_ptr<SparseVector<T>> computeTableauRowPFISparse(const int rowIdx);

  void updateReducedCostsGeneric(const PivotData<T> &pivotData,
                                 const VectorT &pivotRow);
  void updateDualGeneric(const PivotData<T> &pivotData,
                         const VectorT &pivotRow);
  void updateInverseMatrixWithRHSGeneric(const PivotData<T> &pivotData,
                                         const VectorT &enteringColumn);
  void pivotGeneric(const int rowIdx, const int enteringColumnIdx,
                    const VectorT &enteringColumn, const VectorT &pivotRow);
  void pivotImplicitBoundsGeneric(const int pivotRowIdx,
                                  const int enteringColumnIdx,
                                  const VectorT &enteringColumn,
                                  const VectorT &pivotRow,
                                  const bool leavingVarBecomesLowerBound);
  void initTableau();
  void initFullTableau();
  void initBasisMatrixInverse();
  void calculateDual();
  void calculateDualExplicit();
  void calculateDualPFI();
  void calculateDualPFISparse();
  void initMatrixRepresentations();

  void calculateReducedCostsBasedOnDual();
  void calculateReducedCostsFullTableau();
  void calculateRHS();
  void calculateRHSFullTableau();
  void calculateRHSBasisInverse();
  void calculateRHSPFI();
  void calculateRHSPFISparse();
  void calculateRHSWithoutInverse();
  void updateBasisData(const PivotData<T> &pivotData);
  void setObjective(const std::vector<T> &newObjective,
                    const SimplexPhase simplexPhase);
  void updateTableauForNewObjective(const SimplexPhase simplexPhase);
  void multiplyByBasisMatrixLeftInverseUsingPFI(std::vector<T> &vec);
  void multiplyByBasisMatrixLeftInverseUsingPFISparse(SparseVector<T> &vec);
  void
  multiplyByBasisMatrixLeftInverseUsingPFISparseNormal(std::vector<T> &vec);
  void multiplyByBasisMatrixRightInverseUsingPFI(std::vector<T> &vec);
  void multiplyByBasisMatrixRightInverseUsingPFISparse(SparseVector<T> &vec);
  void
  multiplyByBasisMatrixRightInverseUsingPFISparseNormal(std::vector<T> &vec);

  boost::dynamic_bitset<> getIsVariableFreeBitset() const;

  const LinearProgram<T> &_initialProgram;
  std::vector<VariableInfo> _variableInfos;
  std::vector<std::optional<T>> _variableLowerBounds;
  std::vector<std::optional<T>> _variableUpperBounds;
  std::set<std::string> _variableLabelSet;
  std::vector<RowInfo> _rowInfos;

  Matrix<T> _constraintMatrix;
  MatrixRepresentation<T> _constraintMatrixNormalForm;
  SparseMatrixRepresentation<T> _constraintMatrixSparseForm;

  std::vector<T> _rightHandSides;
  std::vector<T> _initialRightHandSides;
  std::vector<T> _objectiveRow;

  SimplexTableauType _simplexTableauType;
  std::vector<std::vector<T>> _basisMatrixInverse;
  std::vector<ElementaryMatrix<T>> _pfiEtms;
  std::vector<SparseElementaryMatrix<T>> _sparsePfiEtms;
  Matrix<T> _fullTableau;

  std::vector<T> _reducedCosts;
  std::vector<T> _y;
  std::vector<T> _x;
  T _objectiveValue{};

  SimplexBasisData _simplexBasisData;

  LPOptimizationResult _result;
};

#endif // GMISOLVER_SIMPLEXTABLEAU_H
