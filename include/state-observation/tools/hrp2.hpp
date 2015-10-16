/**
 * \file     hrp2.hpp
 * \author   Alexis Mifsud
 * \date     2014
 * \brief    Definitions of Hrp2 constants.
 *
 * \details
 *
 *
 */



#ifndef HRP2CONSTANTS
#define HRP2CONSTANTS

//#define STATEOBSERVATION_VERBOUS_CONSTRUCTORS

#include <vector>
#include <deque>

#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
#   include <iostream>
#endif

#include <boost/assert.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace stateObservation
{

    namespace hrp2
    {
        /// mass of the robot
        const double m=56.8;

        /// stifness and damping
        const double linKe=40000;
        const double angKe=400;
        const double linKv=600;
        const double angKv=10;
    }

}

#endif //HRP2CONSTANTS
