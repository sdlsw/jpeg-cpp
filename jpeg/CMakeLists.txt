add_library(jpeg STATIC)
target_link_libraries(jpeg PRIVATE compiler-flags msg)
target_sources(
    jpeg
    PRIVATE
    FILE_SET jpeg TYPE CXX_MODULES FILES
        jpeg.cpp
        concepts.cpp
        data_image.cpp
        data_jpeg.cpp
        data_matrix.cpp
        data_tables.cpp
        tables.cpp
        decode.cpp
        encode.cpp
        file_bmp.cpp
        file_jpeg.cpp
        codingbase.cpp
)
