add_library(arg STATIC)
target_link_libraries(arg PRIVATE compiler-flags)
target_sources(arg PRIVATE FILE_SET arg TYPE CXX_MODULES FILES arg.cpp)
