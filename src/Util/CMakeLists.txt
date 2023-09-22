add_library(MpsReader MpsReader.cpp MpsReader.h)
target_link_libraries(MpsReader EnumTypes LinearProgram SpdlogHeader)

add_library(LPPrinter INTERFACE LPPrinter.h)
target_link_libraries(LPPrinter INTERFACE CommonTypes MatrixTypes SimplexBasisData SpdlogHeader fmt::fmt-header-only)

add_library(SimplexTraits INTERFACE SimplexTraits.h)
target_link_libraries(SimplexTraits INTERFACE MatrixTypes SpdlogHeader)

add_library(SpdlogHeader INTERFACE SpdlogHeader.h)
target_link_libraries(SpdlogHeader INTERFACE spdlog::spdlog)

add_library(LPOptStatistics INTERFACE LPOptStatistics.h)
target_link_libraries(LPOptStatistics INTERFACE EnumTypes)

add_library(IPOptStatistics INTERFACE IPOptStatistics.h)
target_link_libraries(IPOptStatistics INTERFACE LPOptStatistics)

add_library(LPOptStatisticsPrinter INTERFACE LPOptStatisticsPrinter.h)
target_link_libraries(LPOptStatisticsPrinter INTERFACE LPOptStatistics fmt::fmt-header-only)

add_library(GurobiOptimizer INTERFACE GurobiOptimizer.h)
target_link_libraries(GurobiOptimizer INTERFACE LPOptStatistics SpdlogHeader ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARY})
target_include_directories(GurobiOptimizer INTERFACE "${GUROBI_INCLUDE_DIRS}")

set_target_properties(LPPrinter SimplexTraits SpdlogHeader LPOptStatistics IPOptStatistics LPOptStatisticsPrinter GurobiOptimizer
        PROPERTIES LINKER_LANGUAGE CXX)