# simulationcore CMakeLists.txt
# ==================================


file(GLOB headers *.h)
file(GLOB sources *.cpp)
set(simulationcore_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

