cmake_minimum_required(VERSION 3.10)

set(SRC_FILES
        "punch.c"
        "subio.c"
)

set(HDR_FILES
        "punch.h"
)

add_executable(punch ${SRC_FILES})
target_compile_options(punch PUBLIC -Wall -Wextra -Werror)
target_link_libraries(punch m)
