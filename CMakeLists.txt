# Copyright 2013, CNRS
#
# Author: Mehdi Benallegue
#
# This file is part of state-observer.
# state-observer is free software: you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# state-observer is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# state-observer. If not, see <http://www.gnu.org/licenses/>.

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/eigen.cmake)

SET(PROJECT_NAME state-observer)
SET(PROJECT_DESCRIPTION
  "General implementation of observers.")
SET(PROJECT_URL "")

SET(CUSTOM_HEADER_DIR "${PROJECT_NAME}")

SET(${PROJECT_NAME}_HEADERS
  include/state-observer/extended-kalman-filter.hpp
  include/state-observer/kalman-filter-base.hpp
  include/state-observer/kalman-filter.hpp
  include/state-observer/observer-base.hpp
  include/state-observer/zero-delay-observer.hpp
  include/state-observer/template/extended-kalman-filter.hpp
  include/state-observer/template/extended-kalman-filter.hxx
  include/state-observer/template/kalman-filter-base.hpp
  include/state-observer/template/kalman-filter-base.hxx
  include/state-observer/template/kalman-filter.hpp
  include/state-observer/template/kalman-filter.hxx
  include/state-observer/template/observer-base.hpp
  include/state-observer/template/observer-base.hxx
  include/state-observer/template/zero-delay-observer.hpp
  include/state-observer/template/zero-delay-observer.hxx
)
SET(CXX_DISABLE_WERROR TRUE)

SETUP_PROJECT()

SEARCH_FOR_BOOST()
SEARCH_FOR_EIGEN()

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(unit-testings)
SETUP_PROJECT_FINALIZE()
