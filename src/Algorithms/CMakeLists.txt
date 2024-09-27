add_library(DualSimplex DualSimplex.cpp DualSimplex.h)
target_link_libraries(
  DualSimplex
  CommonTypes
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader
  Time)

add_library(DualSimplexPhaseOne DualSimplexPhaseOne.cpp DualSimplexPhaseOne.h)
target_link_libraries(
  DualSimplexPhaseOne
  CommonTypes
  DualSimplex
  LPOptStatistics
  PhaseOneUtilities
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader
  Time)

add_library(PhaseOneUtilities PhaseOneUtilities.cpp PhaseOneUtilities.h)
target_link_libraries(PhaseOneUtilities CommonTypes SimplexBasisData
                      SimplexTraits SimplexTableau)

add_library(PrimalSimplex PrimalSimplex.cpp PrimalSimplex.h)
target_link_libraries(
  PrimalSimplex
  CommonTypes
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexValidator
  SpdlogHeader
  Time)

add_library(PrimalSimplexGomory PrimalSimplexGomory.h PrimalSimplexGomory.cpp
                                PrimalSimplexGomory.h)
target_link_libraries(
  PrimalSimplexGomory
  CommonTypes
  LexicographicOptimizer
  LPOptStatistics
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexTableauResizer
  PrimalSimplex
  Time)

add_library(PrimalSimplexPhaseOne PrimalSimplexPhaseOne.cpp
                                  PrimalSimplexPhaseOne.h)
target_link_libraries(
  PrimalSimplexPhaseOne
  CommonTypes
  LinearProgram
  LPOptStatistics
  PhaseOneUtilities
  PrimalSimplex
  ReinversionManager
  SimplexTraits
  SimplexTableau
  SimplexTableauResizer
  SimplexValidator
  SpdlogHeader
  Time)

add_library(LexicographicOptimizer LexicographicOptimizer.cpp
                                   LexicographicOptimizer.h)
target_link_libraries(LexicographicOptimizer CommonTypes LexReoptStatistics
                      SimplexTraits SimplexTableau PrimalSimplex)

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
  DualSimplex
  Time)

add_library(SimplexTableauResizer SimplexTableauResizer.cpp
                                  SimplexTableauResizer.h)
target_link_libraries(SimplexTableauResizer ReinversionManager SimplexTraits
                      SimplexTableau)

add_library(ReinversionManager ReinversionManager.cpp ReinversionManager.h)
target_link_libraries(ReinversionManager SimplexTraits SimplexTableau)

add_library(SimplexValidator SimplexValidator.cpp SimplexValidator.h)
target_link_libraries(SimplexValidator SimplexTableau LPOptStatistics
                      Overloaded expected)

add_library(SimplexTableau SimplexTableau.cpp SimplexTableau.h)
target_link_libraries(SimplexTableau CommonTypes LinearProgram MatrixTypes
                      SimplexBasisData SimplexTraits)
