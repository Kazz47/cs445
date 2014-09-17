include_directories(${Boost_INCLUDE_DIRS} ${GLOG_INCLUDE_DIRS})

# add_library(slog OBJECT ${SOURCES})
# add_library(slog_static STATIC $<TARGET_OBJECTS:soundex>)
# add_library(slog_shared SHARED $<TARGET_OBJECTS:soundex>)

add_executable(kg_hw1
    hw1)

target_link_libraries(kg_hw1 ${Boost_SYSTEM_LIBRARY} ${GLOG_LIBRARIES})


