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


//Tracking NaN
#include <fenv.h>

namespace stateObservation
{

    namespace hrp2
    {
        /// ATTENTION: spécifique à hrp2: mettre ca ds un fichier a part plus tard.
        /// mass of the robot
        const double m=56.8679920;//59.8;//59.8; // 58 ds la doc
        const double H=1.54;
        const double R=0.31;

//        const Vector3 linKe(53200,53200,53200);
//        const Vector3 angKe(510,510,510);
//        const Vector3 linKv(200,200,200);
//        const Vector3 angKv(20,20,20);

        const double linKe=40000;//53200;
        const double angKe=600;//510;
        const double linKv=600; //200
        const double angKv=60;
    }

}

#endif //HRP2CONSTANTS
