/**
 * \file      accelerometer-gyrometer.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief
 *
 *
 */



#ifndef SIMULATIONACCELEROMETERGYROMETERSENSORHPP
#define SIMULATIONACCELEROMETERGYROMETERSENSORHPP

#include <Eigen/Core>
#include <boost/assert.hpp>

#include <state-observation/sensors-simulation/algebraic-sensor.hpp>

#include <state-observation/sensors-simulation/algorithm/linear-acceleration.hpp>

#include <state-observation/sensors-simulation/algorithm/rotation-velocity.hpp>

namespace stateObservation
{
    /**
     * \class  ObserverBase
     * \brief  The base class for observers.
     *
     * \details
     *
     */

    class AccelerometerGyrometer : public AlgebraicSensor,
        protected algorithm::LinearAcceleration,
        protected algorithm::RotationVelocity
    {
    public:
        virtual ~AccelerometerGyrometer(){}

        virtual unsigned getStateSize() const;

        virtual unsigned getMeasurementSize() const;


    protected:
        virtual Vector computeNoiselessMeasurement_();

        static const int stateSize_= 10;

        static const int measurementSize_=6;

    };

}

#endif //SIMULATIONACCELEROMETERGYROMETERSENSORHPP
