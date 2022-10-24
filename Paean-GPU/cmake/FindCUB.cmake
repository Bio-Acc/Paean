# FindCUB

#  CUB_INCLUDE_DIR - the CUB include directory
#  CUB_VERSION -     Version of CUB in the form "major.minor.patch".

find_path(CUB_INCLUDE_DIR
        HINTS
	/usr/local/cuda/include
	${CUDA_INCLUDE_DIRS}
	${CUB_DIR}
        NAMES cub/cub.cuh
        DOC "CUB headers"
)

# Find CUB version
file(STRINGS ${CUB_INCLUDE_DIR}/cub/version.cuh
        version
        REGEX "#define CUB_VERSION[ \t]+([0-9x]+)"
        )
string(REGEX REPLACE
        "#define CUB_VERSION[ \t]+"
        ""
        version
        "${version}"
        )

string(REGEX MATCH "^[0-9]" major ${version})
string(REGEX REPLACE "^${major}00" "" version "${version}")
string(REGEX MATCH "^[0-9]" minor ${version})
string(REGEX REPLACE "^${minor}0" "" version "${version}")
set(CUB_VERSION "${major}.${minor}.${version}")

# Check for required components
include( FindPackageHandleStandardArgs )
find_package_handle_standard_args(CUB
        REQUIRED_VARS CUB_INCLUDE_DIR
        VERSION_VAR CUB_VERSION
        )