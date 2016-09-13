#TOP LEVEL DIRECTORY


cmake_minimum_required (VERSION 2.8)
project (SimAdapt__Proj)

# path to libraries
set (APPLICATIONS /Users/rogsoares/applications)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_FLAGS "-DPARALLEL -D_SEEKFORBUGS_ -D__ERROR_ANALYSIS_DEBUG__ -D__ADAPTATION_DEBUG__ -Wall")
#set(CMAKE_CXX_FLAGS "-DPARALLEL -O2")  

# FMDB
set (FMDB_INCLUDE "${APPLICATIONS}/FMDB-2011/include" CACHE PATH "Directory where FMDB header files can be found")
set (FMDB_LIB "${APPLICATIONS}/FMDB-2011/lib" CACHE PATH "Directory where the FMDB library can be found")
include_directories(${FMDB_INCLUDE})
link_directories(${FMDB_LIB})
link_libraries(FMDB-O)

# PARMETIS
set (PARMETIS_INCLUDE "${APPLICATIONS}/ParMetis-3.1/include" CACHE PATH "Directory where FMDB header files can be found")
set (PARMETIS_LIB "${APPLICATIONS}/ParMetis-3.1/lib" CACHE PATH "Directory where the FMDB library can be found")
include_directories(${PARMETIS_INCLUDE})
link_directories(${PARMETIS_LIB})
link_libraries(parmetis metis)

# AUTOPACK
set (AUTOPACK_INCLUDE "${APPLICATIONS}/autopack/include" CACHE PATH "Directory where FMDB header files can be found")
set (AUTOPACK_LIB "${APPLICATIONS}/autopack/lib" CACHE PATH "Directory where the FMDB library can be found")
include_directories(${AUTOPACK_INCLUDE})
link_directories(${AUTOPACK_LIB})
link_libraries(autopack-O)

# PETSC
set(PETSC_DIR "${APPLICATIONS}/petsc")
set(PETSC_INCLUDE "${PETSC_DIR}/include" "${PETSC_DIR}/include/private" CACHE PATH "Directory where PETSc header files can be found")
set(PETSC_LIB "${PETSC_DIR}/lib" CACHE PATH "Directory where the PETSc library can be found" )
include_directories(${PETSC_INCLUDE})
link_directories(${PETSC_LIB})
link_libraries(petscmat f2clapack petscksp petsc petscvec petscdm petscsnes petscts f2cblas petsccontrib)


# PROJECT
include_directories(
                    "${PROJECT_SOURCE_DIR}/include" 
                    "${PROJECT_SOURCE_DIR}/src"
                    "${PROJECT_SOURCE_DIR}/src/pre-processor"
                    "${PROJECT_SOURCE_DIR}/src/elliptic"
                    "${PROJECT_SOURCE_DIR}/src/adaptation"
                    "${PROJECT_SOURCE_DIR}/src/adaptation/adaptative-remeshing"
                    "${PROJECT_SOURCE_DIR}/src/adaptation/h-refinement"
                    "${PROJECT_SOURCE_DIR}/src/high-order"
                    "${PROJECT_SOURCE_DIR}/src/hyperbolic"
                    "${PROJECT_SOURCE_DIR}/src/tools"
                    "${PROJECT_SOURCE_DIR}/src/interpolation"
                    "${PROJECT_SOURCE_DIR}/src/error-analysis"
                    "${PROJECT_SOURCE_DIR}/src/SIMULATION_core"
                    )
                                    
add_subdirectory(src/SIMULATION_core)
add_subdirectory(src/pre-processor)
add_subdirectory(src/elliptic)
add_subdirectory(src/hyperbolic)
add_subdirectory(src/error-analysis)
add_subdirectory(src/adaptation)
add_subdirectory(src/adaptation/adaptative-remeshing)
add_subdirectory(src/adaptation/h-refinement)
add_subdirectory(src/high-order)
add_subdirectory(src/tools)
add_subdirectory(src/interpolation)
add_subdirectory(src)
