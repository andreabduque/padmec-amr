# src CMakeLists.txt
# ==================================

set(prepro_sources ${preprocessors_STAT_SRCS})
add_library(preprocessorslib ${prepro_sources})

set(elliptic_sources ${elliptic_STAT_SRCS})
add_library(ellipticlib ${elliptic_sources})

set(hyperbolic_sources ${hyperbolic_STAT_SRCS})
add_library(hyperboliclib ${hyperbolic_sources})

set(erroranalysis_sources ${erroranalysis_STAT_SRCS})
add_library(erroranalysislib ${erroranalysis_sources})

#set(adaptation_sources ${adaptation_STAT_SRCS})
#add_library(adaptationlib ${adaptation_sources})

set(highorder_sources ${highorder_STAT_SRCS})
add_library(highorderlib ${highorder_sources})

set(interpolation_sources ${interpolation_STAT_SRCS})
add_library(interpolationlib ${interpolation_sources})

set(tools_sources ${tools_STAT_SRCS})
add_library(toolslib ${tools_sources})

set(simulationcore_sources ${simulationcore_STAT_SRCS})
add_library(simulationcorelib ${simulationcore_sources})

set(hrefinement_sources ${hrefinement_STAT_SRCS})
add_library(hrefinementlib ${hrefinement_sources})

set(adaptativeremeshing_sources ${adaptativeremeshing_STAT_SRCS})
add_library(adaptativeremeshinglib ${adaptativeremeshing_sources})

file(GLOB sources *.cpp)
set(SRCS ${sources} ${SRCS})

file(GLOB headers *.h)
set(HDRS ${headers} ${HDRS})

# ============================================================================================================
# Create main executable
# ============================================================================================================
add_executable(PADMEC_AMR main.cpp ${HDRS})
target_link_libraries(PADMEC_AMR preprocessorslib ellipticlib hyperboliclib hrefinementlib adaptativeremeshinglib erroranalysislib highorderlib interpolationlib toolslib simulationcorelib)


# ============================================================================================================
# Create adaptation executable
# ============================================================================================================
# add_executable(SimAdapt__Proj__testingAdaptation main.cpp ${HDRS})
# target_link_libraries(SimAdapt__Proj__testingAdaptation preprocessorslib ellipticlib hyperboliclib adaptationlib erroranalysislib highorderlib interpolationlib toolslib simulationcorelib)

