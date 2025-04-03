add_library(jpeg STATIC)
target_link_libraries(jpeg PRIVATE compiler-flags msg)
target_sources(
    jpeg
    PRIVATE
    FILE_SET jpeg TYPE CXX_MODULES FILES
        jpeg.cpp
        data.cpp
        tables.cpp
        decode.cpp
        encode.cpp
        file.cpp
        codingbase.cpp
)
