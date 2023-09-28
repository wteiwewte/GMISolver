add_library(RevisedDualSimplexPFIBounds RevisedDualSimplexPFIBounds.cpp
                                        RevisedDualSimplexPFIBounds.h)
target_link_libraries(
  RevisedDualSimplexPFIBounds
  CommonTypes
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader)

add_library(RevisedPrimalSimplexPFIBounds RevisedPrimalSimplexPFIBounds.cpp
                                          RevisedPrimalSimplexPFIBounds.h)
target_link_libraries(
  RevisedPrimalSimplexPFIBounds
  CommonTypes
  LinearProgram
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader)

add_library(LexicographicOptimizer LexicographicOptimizer.cpp
                                   LexicographicOptimizer.h)
target_link_libraries(
  LexicographicOptimizer CommonTypes LPOptStatistics SimplexTraits
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
  RevisedDualSimplexPFIBounds)

add_library(ReinversionManager ReinversionManager.cpp ReinversionManager.h)
target_link_libraries(ReinversionManager SimplexTraits SimplexTableau)

add_library(SimplexValidator SimplexValidator.cpp SimplexValidator.h)
target_link_libraries(SimplexValidator SimplexTableau LPOptStatistics expected)

add_library(SimplexTableau SimplexTableau.cpp SimplexTableau.h)
target_link_libraries(SimplexTableau CommonTypes LinearProgram MatrixTypes
                      SimplexBasisData SimplexTraits)
