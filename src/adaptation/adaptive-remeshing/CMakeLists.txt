# multigrid CMakeLists.txt
# ==================================


file(GLOB headers *.h)
file(GLOB sources *.cpp)
set(adaptativeremeshing_STAT_SRCS ${headers} ${sources} PARENT_SCOPE)

