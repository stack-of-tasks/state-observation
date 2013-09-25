
# provide path to dependency libraries

ADD_EXECUTABLE(test-kalman-filter
  test-kalman-filter.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

LINK_DIRECTORIES(${${PROJECT_NAME}_BINARY_DIR}/src)
TARGET_LINK_LIBRARIES(test-kalman-filter
  ${Boost_LIBRARIES} ${PROJECT_NAME})



ADD_TEST(test-kalman-filter test-kalman-filter)

ADD_EXECUTABLE(test-kalman-filter-template
  test-kalman-filter-template.cpp)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

TARGET_LINK_LIBRARIES(test-kalman-filter ${Boost_LIBRARIES})

ADD_TEST(test-kalman-filter-template test-kalman-filter-template)
