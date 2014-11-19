include_directories(${Boost_INCLUDE_DIRS} ${GLOG_INCLUDE_DIRS})

add_executable(kg_hw1 hw1)
target_link_libraries(kg_hw1 ${Boost_SYSTEM_LIBRARY} ${GLOG_LIBRARIES})

add_executable(kg_hw2 hw2)
target_link_libraries(kg_hw2 ${Boost_SYSTEM_LIBRARY} ${GLOG_LIBRARIES})

add_executable(snowflake snowflake)
target_link_libraries(snowflake ${Boost_SYSTEM_LIBRARY} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${GLOG_LIBRARIES})

add_executable(kg_dna dna_sim)
target_link_libraries(kg_dna ${Boost_SYSTEM_LIBRARY} ${GLOG_LIBRARIES})

add_executable(kg_nbody nbody_sim)
target_link_libraries(kg_nbody ${GLOG_LIBRARIES})

add_executable(kg_fifo fifo_sim)
target_link_libraries(kg_fifo ${GLOG_LIBRARIES})

