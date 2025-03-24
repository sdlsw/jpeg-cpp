add_library(msg STATIC)
target_link_libraries(msg PRIVATE compiler-flags)
target_sources(
    msg
    PRIVATE
    FILE_SET msg TYPE CXX_MODULES FILES msg.cpp loglevel_static.cpp
)
