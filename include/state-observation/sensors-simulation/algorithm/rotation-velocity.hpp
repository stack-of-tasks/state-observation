/**
 * \file      rotation-velocity.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 * \details
 *
 *
 */



#ifndef SENSORALGORITHMSROTATIONVELOCITYHPP
#define SENSORALGORITHMSROTATIONVELOCITYHPP

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
     namespace algorithm
    {
        /**
         * \class  RotationVelocity
         * \brief
         *
         */

        class RotationVelocity
        {
        public:
            Vector3 rotationVelocityMeasure(const Vector3 & rotationVector, const Matrix3 & orientation) const;


        protected:

        };
    }

}

#endif //SENSORALGORITHMSROTATIONVELOCITYHPP
