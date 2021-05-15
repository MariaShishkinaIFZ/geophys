cmake_minimum_required(VERSION 3.10)

set(SRC_FILES
        "3d-seis-w19-2.cpp"
        "haar.cpp"
        "randgen.cpp"
        "spvec.cpp"
        "subio.cpp"
        "submem.cpp"
        "bicsb.cpp"
        "nrutil.cpp"
        "csc.cpp"
        "lsqr.cpp"
)

set(SRC_C_FILES
        "subio.c"
        "submem.c"
)

set(HDR_FILES
        "cilk_util.h"
        "csc.h"
        "lsqr.h"
        "mortoncompare.h"
        "utility.h"
        "datatypes.h"
        "haar.h"
        "randgen.h"
        "spvec.h"
        "subio.h"
        "submem.h"
        "bicsb.h"
        "gmtio.h"
        "matmul.h"
        "nrutil.h"
        "triple.h"
)

set(PARAMETER_FILES
        "sample.par"
)

add_executable(inv-19-2 ${SRC_FILES})
target_compile_options(inv-19-2 PUBLIC -Wall -Wextra -Werror)
target_link_libraries(inv-19-2 m)
