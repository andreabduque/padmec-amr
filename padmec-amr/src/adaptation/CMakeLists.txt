# multigrid CMakeLists.txt
# ==================================


file(GLOB headers *.h)
file(GLOB sources *.cpp)
set(adaptation_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

