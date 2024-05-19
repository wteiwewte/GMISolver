add_library(RevisedDualSimplexPFIBounds DualSimplex.cpp DualSimplex.h)
target_link_libraries(
  RevisedDualSimplexPFIBounds
  CommonTypes
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader
  Time)

add_library(RevisedPrimalSimplexPFIBounds PrimalSimplex.cpp PrimalSimplex.h)
target_link_libraries(
  RevisedPrimalSimplexPFIBounds
  CommonTypes
  LinearProgram
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexTableauResizer
  SimplexValidator
  SpdlogHeader
  Time)

add_library(LexicographicOptimizer LexicographicOptimizer.cpp
                                   LexicographicOptimizer.h)
target_link_libraries(
  LexicographicOptimizer CommonTypes LexReoptStatistics SimplexTraits
  SimplexTableau RevisedPrimalSimplexPFIBounds)

add_library(DualSimplexGomory DualSimplexGomory.h DualSimplexGomory.cpp)
target_link_libraries(
  DualSimplexGomory
  CommonTypes
  LexicographicOptimizer
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexTableauResizer
  RevisedDualSimplexPFIBounds
  Time)

add_library(SimplexTableauResizer SimplexTableauResizer.cpp
                                  SimplexTableauResizer.h)
target_link_libraries(SimplexTableauResizer ReinversionManager SimplexTraits
                      SimplexTableau)

add_library(ReinversionManager ReinversionManager.cpp ReinversionManager.h)
target_link_libraries(ReinversionManager SimplexTraits SimplexTableau)

add_library(SimplexValidator SimplexValidator.cpp SimplexValidator.h)
target_link_libraries(SimplexValidator SimplexTableau LPOptStatistics expected)

add_library(SimplexTableau SimplexTableau.cpp SimplexTableau.h)
target_link_libraries(SimplexTableau CommonTypes LinearProgram MatrixTypes
                      SimplexBasisData SimplexTraits)
