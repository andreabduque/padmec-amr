# multigrid CMakeLists.txt
# ==================================


file(GLOB headers *.h)
file(GLOB sources *.cpp)
set(preprocessors_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

