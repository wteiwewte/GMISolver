#ifndef GMISOLVER_SIMPLEXTABLEAU_H
#define GMISOLVER_SIMPLEXTABLEAU_H

#include "src/DataModel/CommonTypes.h"
#include "src/DataModel/MatrixTypes.h"
#include "src/DataModel/SimplexBasisData.h"
#include "src/Util/SimplexTraits.h"

#include <set>

template <typename T> class LinearProgram;
template <typename T, typename SimplexTraitsT>
class RevisedPrimalSimplexPFIBounds;
template <typename T, typename SimplexTraitsT>
class RevisedPrimalSimplexPFIBounds;
template <typename T, typename SimplexTraitsT>
class RevisedDualSimplexPFIBounds;

template <typename T, typename SimplexTraitsT = SimplexTraits<T>>
class SimplexTableau {
public:
  SimplexTableau(const LinearProgram<T> &linearProgram,
                 const bool isPrimalSimplex,
                 const bool useProductFormOfInverse);

  void convertToStandardForm();
  void makeRightHandSidesNonNegative();
  void addBoundsToMatrix();
  void addArtificialVariables();
  void init(const bool isPrimalSimplex);

  std::vector<T> artificialObjective() const;
  std::vector<T> originalObjective() const;
  void calculateCurrentObjectiveValue();
  void calculateSolution();

  std::string toString() const;
  std::string toStringObjectiveValue() const;
  std::string toStringSolution() const;
  std::string toStringLpSolveFormat() const;

  const std::vector<VariableInfo> &getVariableInfos() const {
    return _variableInfos;
  }
  const std::vector<RowInfo> &getRowInfos() const { return _rowInfos; }

private:
  friend class RevisedPrimalSimplexPFIBounds<T, SimplexTraitsT>;
  friend class RevisedDualSimplexPFIBounds<T, SimplexTraitsT>;

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

  std::optional<SimplexBasisData> createBasisFromArtificialVars() const;

  int basicColumnIdx(const int rowIdx) const {
    return _simplexBasisData._rowToBasisColumnIdxMap[rowIdx];
  }
  bool isColumnAllowedToEnterBasis(const int colIdx);
  std::optional<T> curSatisfiedBound(const int varIdx);

  auto computeTableauColumnGeneric(const int colIdx) {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return computeTableauColumnPFISparse(colIdx);
    else
      return computeTableauColumn(colIdx);
  }

  std::vector<T> computeTableauColumn(const int colIdx);
  std::vector<T> computeTableauColumnExplicit(const int colIdx);
  std::vector<T> computeTableauColumnPFI(const int colIdx);
  SparseVector<T> computeTableauColumnPFISparse(const int colIdx);
  std::vector<T> computeTableauRow(const int rowIdx);
  std::vector<T> computeTableauRowExplicit(const int rowIdx);
  std::vector<T> computeTableauRowPFI(const int rowIdx);
  SparseVector<T> computeTableauRowPFISparse(const int rowIdx);

  auto computeTableauRowGeneric(const int rowIdx) {
    if constexpr (SimplexTraitsT::useSparseRepresentationValue)
      return computeTableauRowPFISparse(rowIdx);
    else
      return computeTableauRow(rowIdx);
  }

  void updateReducedCostsGeneric(const PivotData<T> &pivotData,
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

  void initBasisMatrixInverse();
  void calculateDual();
  void calculateDualExplicit();
  void calculateDualPFI();
  void calculateDualPFISparse();
  void initBoundsForDualSimplex();
  void initMatrixRepresentations();

  void calculateReducedCostsBasedOnDual();
  void calculateRHS();
  void calculateRHSExplicit();
  void calculateRHSPFI();
  void calculateRHSPFISparse();
  void calculateRHSWithoutInverse();
  void updateBasisData(const PivotData<T> &pivotData);
  bool reinversion();
  bool reinversionExplicit();
  bool reinversionPFI();
  bool reinversionPFISparse();
  void setObjective(const std::vector<T> &newObjective);
  void multiplyByBasisMatrixLeftInverseUsingPFI(std::vector<T> &vec);
  void multiplyByBasisMatrixLeftInverseUsingPFISparse(SparseVector<T> &vec);
  void
  multiplyByBasisMatrixLeftInverseUsingPFISparseNormal(std::vector<T> &vec);
  void multiplyByBasisMatrixRightInverseUsingPFI(std::vector<T> &vec);
  void multiplyByBasisMatrixRightInverseUsingPFISparse(SparseVector<T> &vec);
  void
  multiplyByBasisMatrixRightInverseUsingPFISparseNormal(std::vector<T> &vec);

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
  std::vector<T> _initialObjectiveRow;
  std::vector<T> _objectiveRow;

  bool _useProductFormOfInverse;
  std::vector<std::vector<T>> _basisMatrixInverse;
  std::vector<ElementaryMatrix<T>> _pfiEtms;
  std::vector<SparseElementaryMatrix<T>> _sparsePfiEtms;
  std::vector<T> _reducedCosts;
  std::vector<T> _y;
  std::vector<T> _x;
  T _objectiveValue{};

  SimplexBasisData _simplexBasisData;

  LPOptimizationResult _result;
};

#endif // GMISOLVER_SIMPLEXTABLEAU_H
