add_library(Algorithms Algorithms/Simplex.cpp Algorithms/Simplex.h)
target_link_libraries(Algorithms DataModel)

add_library(DataModel DataModel/LinearProgram.cpp DataModel/LinearProgram.h DataModel/EnumTypes.h DataModel/EnumTypes.cpp DataModel/SimplexTableu.cpp DataModel/SimplexTableu.h DataModel/CommonTypes.h DataModel/CommonConstants.h)
target_link_libraries(DataModel spdlog::spdlog fmt::fmt-header-only)

add_library(Util Util/MpsReader.cpp Util/MpsReader.h Util/LPPrinter.cpp Util/LPPrinter.h)
target_link_libraries(Util DataModel)
