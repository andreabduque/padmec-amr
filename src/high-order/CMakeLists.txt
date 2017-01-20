# multigrid CMakeLists.txt
# ==================================


file(GLOB headers *.h)
file(GLOB sources *.cpp)
set(highorder_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

