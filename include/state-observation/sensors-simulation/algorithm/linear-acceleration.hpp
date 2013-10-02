/**
 * \file      linear-acceleration.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 * \details
 *
 *
 */



#ifndef SENSORALGORITHMSLINEARACCELERATIONHPP
#define SENSORALGORITHMSLINEARACCELERATIONHPP

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    namespace algorithm
    {
        /**
         * \class  LinearAcceleration
         * \brief
         *
         */

        class LinearAcceleration
        {
        public:
            Vector3 accelerationMeasure(const Vector3 & acceleration, const Matrix3 & orientation) const;


        protected:

        };
    }

}

#endif //SENSORALGORITHMSLINEARACCELERATIONHPP
