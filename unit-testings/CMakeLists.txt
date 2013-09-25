# Francois Keith, Mehdi Benallegue JRL, CNRS/AIST
# Creation: 2013
# History:
#
# Copyright CNRS/AIST



# prefix and suffix each element of list by ${prefix}elemnt${suffix}
macro(ADDPREFIX newlist prefix list_name)
    # create empty list - necessary?
    SET(${newlist})

     # prefix and suffix elements
    foreach(l ${${list_name}})
      list(APPEND ${newlist} ${prefix}${l} )
    endforeach()

endmacro(ADDPREFIX)


INCLUDE(CTest)

ADDPREFIX(${PROJECT_NAME}_ABSOLUTE_HEADERS  "${CMAKE_SOURCE_DIR}/"  ${PROJECT_NAME}_HEADERS)

SET(${PROJECT_NAME}_test_CXXFLAGS "-I${PROJECT_SOURCE_DIR}/include/ ${${PROJECT_NAME}_CXXFLAGS}")


 include_directories(${Boost_INCLUDE_DIR} ${PROJECT_INCLUDE_DIR})
 
 ADD_EXECUTABLE(kalman-testing test-kalman-filter.cpp)
ADD_EXECUTABLE(kalman-testing-template test-kalman-filter-template.cpp)


 
 target_link_libraries(kalman-testing ${PROJECT_NAME} ${Boost_LIBRARIES})
 target_link_libraries(kalman-testing-template ${Boost_LIBRARIES})
 
 ADD_DEPENDENCIES(kalman-testing ${PROJECT_NAME}) 

 ADD_TEST(kalman-testing kalman-testing)
 ADD_TEST(kalman-testing kalman-testing-template)
