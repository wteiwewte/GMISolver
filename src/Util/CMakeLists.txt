add_library(CommonFunctions INTERFACE CommonFunctions.h)

add_library(MpsReader MpsReader.cpp MpsReader.h)
target_link_libraries(MpsReader CommonFunctions EnumTypes LinearProgram
                      SpdlogHeader)

add_library(MatchingEdgeFormatReader MatchingEdgeFormatReader.cpp
                                     MatchingEdgeFormatReader.h)
target_link_libraries(MatchingEdgeFormatReader CommonFunctions EnumTypes
                      LinearProgram SpdlogHeader)

add_library(LPPrinter INTERFACE LPPrinter.h)
target_link_libraries(
  LPPrinter INTERFACE CommonTypes MatrixTypes SimplexBasisData SpdlogHeader
                      fmt::fmt-header-only)

add_library(SimplexTraits INTERFACE SimplexTraits.h)
target_link_libraries(SimplexTraits INTERFACE EnumTypes MatrixTypes
                                              SpdlogHeader)

add_library(SpdlogHeader INTERFACE SpdlogHeader.h)
target_link_libraries(SpdlogHeader INTERFACE spdlog::spdlog)

add_library(LexReoptStatistics INTERFACE LexReoptStatistics.h)
target_link_libraries(LexReoptStatistics INTERFACE EnumTypes LPOptStatistics)

add_library(LPOptStatistics INTERFACE LPOptStatistics.h)
target_link_libraries(LPOptStatistics INTERFACE EnumTypes)

add_library(IPOptStatistics INTERFACE IPOptStatistics.h)
target_link_libraries(IPOptStatistics INTERFACE LexReoptStatistics
                                                LPOptStatistics)

add_library(OptStatisticsPrinter INTERFACE OptStatisticsPrinter.h)
target_link_libraries(
  OptStatisticsPrinter INTERFACE LPOptStatistics LexReoptStatistics
                                 IPOptStatistics fmt::fmt-header-only)
add_library(Overloaded INTERFACE Overloaded.h)

add_library(GurobiOptimizer INTERFACE GurobiOptimizer.h)
target_link_libraries(
  GurobiOptimizer INTERFACE LPOptStatistics SpdlogHeader ${GUROBI_CXX_LIBRARY}
                            ${GUROBI_LIBRARY})
target_include_directories(GurobiOptimizer INTERFACE "${GUROBI_INCLUDE_DIRS}")

add_library(Time INTERFACE Time.h)

set_target_properties(
  LPPrinter
  SimplexTraits
  SpdlogHeader
  LPOptStatistics
  IPOptStatistics
  OptStatisticsPrinter
  GurobiOptimizer
  PROPERTIES LINKER_LANGUAGE CXX)
