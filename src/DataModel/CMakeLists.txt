add_library(LinearProgram LinearProgram.cpp LinearProgram.h)
target_link_libraries(LinearProgram CommonConstants CommonTypes EnumTypes MatrixTypes LPPrinter)

add_library(EnumTypes EnumTypes.h EnumTypes.cpp)
target_link_libraries(EnumTypes absl::flags)

add_library(CommonTypes INTERFACE CommonTypes.h)
target_link_libraries(CommonTypes INTERFACE EnumTypes)

add_library(CommonConstants INTERFACE CommonConstants.h)
add_library(SimplexBasisData INTERFACE SimplexBasisData.h)
add_library(MatrixTypes INTERFACE MatrixTypes.h)

set_target_properties(CommonTypes CommonConstants MatrixTypes SimplexBasisData PROPERTIES LINKER_LANGUAGE CXX)