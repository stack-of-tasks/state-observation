
# provide path to dependency libraries

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)
LINK_DIRECTORIES(${${PROJECT_NAME}_BINARY_DIR}/src)

ADD_EXECUTABLE(test-kalman-filter test-kalman-filter.cpp)
ADD_EXECUTABLE(test-kalman-filter-template test-kalman-filter-template.cpp)
ADD_EXECUTABLE(imu-test imu-test.cpp ${${PROJECT_NAME}_ABSOLUTE_HEADERS})
ADD_EXECUTABLE(HRP2LAAS-IMU-Observation HRP2LAAS-IMU-Observation.cpp ${${PROJECT_NAME}_ABSOLUTE_HEADERS})
ADD_EXECUTABLE(flex-estimation-test flex-estimation-test.cpp ${${PROJECT_NAME}_ABSOLUTE_HEADERS})
ADD_EXECUTABLE(real-offline-flexibility-estimator real-offline-flexibility-estimator.cpp ${${PROJECT_NAME}_ABSOLUTE_HEADERS})

TARGET_LINK_LIBRARIES(test-kalman-filter  ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(test-kalman-filter ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(test-kalman-filter-template ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(imu-test ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(HRP2LAAS-IMU-Observation ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(flex-estimation-test ${Boost_LIBRARIES} ${PROJECT_NAME})
TARGET_LINK_LIBRARIES(real-offline-flexibility-estimator ${Boost_LIBRARIES} ${PROJECT_NAME})

ADD_TEST(test-kalman-filter test-kalman-filter)
ADD_TEST(test-kalman-filter-template test-kalman-filter-template)
ADD_TEST(imu-test imu-test)
ADD_TEST(flex-estimation-test flex-estimation-test)
ADD_TEST(HRP2LAAS-IMU-Observation HRP2LAAS-IMU-Observation)
#ADD_TEST(real-offline-flexibility-estimator real-offline-flexibility-estimator)
