# Copyright 2008, 2009, 2010, 2013 Olivier Stasse, Mehdi Benallegue JRL, CNRS/AIST
#
# This file is part of state-observation.
# state-observation is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# state-observation is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Lesser Public License for more details.
# You should have received a copy of the GNU Lesser General Public License
# along with state-observation.  If not, see <http://www.gnu.org/licenses/>.

# prefix and suffix each element of list by ${prefix}elemnt${suffix}
macro(ADDPREFIX newlist prefix list_name)
    # create empty list - necessary?
    SET(${newlist})

     # prefix and suffix elements
    foreach(l ${${list_name}})
      list(APPEND ${newlist} ${prefix}${l} )
    endforeach()

endmacro(ADDPREFIX)

ADDPREFIX(${PROJECT_NAME}_ABSOLUTE_HEADERS  "${CMAKE_SOURCE_DIR}/"  ${PROJECT_NAME}_HEADERS)

SET(LIBRARY_NAME ${PROJECT_NAME})

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/include)

ADD_LIBRARY(${LIBRARY_NAME}
  SHARED
  kalman-filter-base.cpp  
  observer-base.cpp
  extended-kalman-filter.cpp
  kalman-filter.cpp
  zero-delay-observer.cpp
  sensor-base.cpp
  algebraic-sensor.cpp
  rotation-velocity.cpp
  linear-acceleration.cpp
  definitions.cpp
  accelerometer-gyrometer.cpp
  probability-law-simulation.cpp
  gaussian-white-noise.cpp
  dynamical-system-functor-base.cpp
  dynamical-system-simulator.cpp
  ${${PROJECT_NAME}_ABSOLUTE_HEADERS}
  )

SET_TARGET_PROPERTIES(${LIBRARY_NAME}
  PROPERTIES
  SOVERSION ${PROJECT_VERSION}
  INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)

INSTALL(TARGETS ${LIBRARY_NAME}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)


