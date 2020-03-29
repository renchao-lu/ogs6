######################
### Find tools     ###
######################

string(REPLACE ".windows.1" "" GIT_VERSION_STRING ${GIT_VERSION_STRING})
if(${GIT_VERSION_STRING} VERSION_LESS ${ogs.minimum_version.git})
    message(FATAL_ERROR "Git version ${ogs.minimum_version.git} is required. \
        Found version ${GIT_VERSION_STRING}.")
endif()

find_package(Doxygen OPTIONAL_COMPONENTS dot)

# Find gnu profiler gprof
find_program(GPROF_PATH gprof DOC "GNU profiler gprof" QUIET)

find_program(CPPCHECK_TOOL_PATH cppcheck)

# Find bash itself ...
find_program(BASH_TOOL_PATH bash
    HINTS ${GITHUB_BIN_DIR} DOC "The bash executable")

# Dumpbin is a windows dependency analaysis tool required for packaging.
# Variable has to be named gp_cmd to override the outdated find routines
# of the GetPrerequisites CMake-module.
if(WIN32)
    include(MSVCPaths)
    find_program(gp_cmd dumpbin DOC "Windows dependency analysis tool"
        PATHS ${MSVC_INSTALL_PATHS} PATH_SUFFIXES VC/bin)
    if(gp_cmd)
        get_filename_component(dir ${gp_cmd} PATH)
        set(ENV{PATH} "${dir}/../../../Common7/IDE;$ENV{PATH}")
    endif()
endif()

find_program(CURL_TOOL_PATH curl DOC "The curl-tool")

find_program(S3CMD_TOOL_PATH s3cmd DOC "S3cmd tool for uploading to Amazon S3")

find_program(CCACHE_TOOL_PATH ccache)

# Tools for web
find_program(VTKJS_CONVERTER vtkDataConverter
    PATHS ${PROJECT_SOURCE_DIR}/web/node_modules/.bin)
find_program(HUGO hugo)
find_program(NPM npm)
find_program(YARN yarn)
find_program(PIP pip)
find_program(PANDOC_CITEPROC pandoc-citeproc)

find_program(MODULE_CMD modulecmd
    PATHS /usr/local/modules/3.2.10-1/Modules/3.2.10/bin)

######################
### Find libraries ###
######################
find_package(Boost ${ogs.minimum_version.boost} REQUIRED)

set(VTK_COMPONENTS vtkIOXML)
if(OGS_BUILD_GUI)
    set(VTK_COMPONENTS ${VTK_COMPONENTS}
        vtkIOImage vtkIOLegacy vtkIOExport
        vtkIOExportOpenGL2 vtkInteractionStyle vtkInteractionWidgets
        vtkGUISupportQt vtkRenderingOpenGL2 vtkRenderingContextOpenGL2
        vtkFiltersTexture vtkRenderingCore vtkFiltersParallel
    )
endif()
if(OGS_USE_MPI)
    set(VTK_COMPONENTS ${VTK_COMPONENTS} vtkIOParallelXML vtkParallelMPI)
endif()
find_package(VTK ${ogs.minimum_version.vtk} REQUIRED COMPONENTS ${VTK_COMPONENTS})
include(${VTK_USE_FILE})

find_package(Eigen3 ${ogs.minimum_version.eigen} REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

# Do not search for libs if this option is set
if(OGS_NO_EXTERNAL_LIBS)
    return()
endif() # OGS_NO_EXTERNAL_LIBS

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(Metis QUIET)

## Qt5 library ##
if(OGS_BUILD_GUI)
    set(QT_MODULES Gui Widgets Xml XmlPatterns)
    if(OGS_USE_CONAN AND UNIX AND NOT APPLE)
        set(QT_MODULES ${QT_MODULES} X11Extras)
    endif()
    find_package(Qt5 ${ogs.minimum_version.qt} REQUIRED ${QT_MODULES})
    cmake_policy(SET CMP0020 NEW)
endif()

if(OGS_USE_NETCDF)
    set(NETCDF_ROOT ${CONAN_NETCDF-C_ROOT})
    set(NETCDF_CXX_ROOT ${CONAN_NETCDF-CXX_ROOT})
    find_package(NetCDF REQUIRED)
    if(OGS_USE_CONAN AND MSVC)
        set(HDF5_DIR ${CONAN_HDF5_ROOT}/cmake/hdf5)
        find_package(HDF5 REQUIRED COMPONENTS C HL NO_MODULE)
    else()
        find_package(HDF5 REQUIRED COMPONENTS C HL)
    endif()
    add_compile_options(-DOGS_USE_NETCDF)
endif()

# lapack
find_package(LAPACK QUIET)

## geotiff ##
find_package(LibGeoTiff)
if(GEOTIFF_FOUND)
    add_definitions(-DGEOTIFF_FOUND)
endif() # GEOTIFF_FOUND

## lis ##
if(OGS_USE_LIS)
    find_package( LIS REQUIRED )
endif()

if(OGS_USE_MKL)
    find_package( MKL REQUIRED )
endif()

if(OGS_USE_PETSC)
    message(STATUS "Configuring for PETSc")

    option(FORCE_PETSC_EXECUTABLE_RUNS
        "Force CMake to accept a given PETSc configuration" ON)

    # Force CMake to accept a given PETSc configuration in case the failure of
    # MPI tests. This may cause the compilation broken.
    if(FORCE_PETSC_EXECUTABLE_RUNS)
        set(PETSC_EXECUTABLE_RUNS YES)
    endif()

    find_package(PETSc ${ogs.minimum_version.petsc} REQUIRED)

    include_directories(SYSTEM ${PETSC_INCLUDES})

    add_definitions(-DPETSC_VERSION_NUMBER=PETSC_VERSION_MAJOR*1000+PETSC_VERSION_MINOR*10)

endif()

find_package(OpenSSL)

## Check MPI package
if(OGS_USE_MPI)
    find_package(MPI REQUIRED)
endif()

find_package(Shapelib)
if(Shapelib_FOUND)
    include_directories(SYSTEM ${Shapelib_INCLUDE_DIRS})
elseif(OGS_BUILD_GUI)
    message(FATAL_ERROR "Shapelib not found but it is required for OGS_BUILD_GUI!")
endif()

## Sundials cvode ode-solver library
if(OGS_USE_CVODE)
    find_package(CVODE REQUIRED)
    add_definitions(-DCVODE_FOUND)
endif()
