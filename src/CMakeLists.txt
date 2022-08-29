add_library(Algorithms Algorithms/PrimalSimplex.cpp Algorithms/PrimalSimplex.h)
target_link_libraries(Algorithms DataModel Util)

add_library(DataModel DataModel/LinearProgram.cpp DataModel/LinearProgram.h DataModel/EnumTypes.h DataModel/EnumTypes.cpp Algorithms/SimplexTableau.cpp Algorithms/SimplexTableau.h DataModel/CommonTypes.h DataModel/CommonConstants.h DataModel/SimplexBasisData.h)
target_link_libraries(DataModel spdlog::spdlog fmt::fmt-header-only)

add_library(Util Util/MpsReader.cpp Util/MpsReader.h Util/LPPrinter.cpp Util/LPPrinter.h Util/ComparisonTraits.h)
target_link_libraries(Util DataModel)
