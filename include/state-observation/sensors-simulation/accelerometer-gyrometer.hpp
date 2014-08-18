/**
 * \file      accelerometer-gyrometer.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief     Implements the accelerometer-gyrometer inertial measuremen
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
     * \class  AccelerometerGyrometer
     * \brief  Implements the accelerometer-gyrometer measurements
     *
     *
     *
     * \details
     *
     */

    class AccelerometerGyrometer : public AlgebraicSensor,
        protected algorithm::LinearAcceleration,
        protected algorithm::RotationVelocity
    {
    public:
        AccelerometerGyrometer()
        {
#ifdef STATEOBSERVATION_VERBOUS_CONSTRUCTORS
           std::cout<<std::endl<<"AccelerometerGyrometer Constructor"<<std::endl;
#endif //STATEOBSERVATION_VERBOUS_CONSTRUCTOR

       }

        ///Virtual destructor
        virtual ~AccelerometerGyrometer(){}

        ///Gets the state vector Size
        virtual unsigned getStateSize() const;

        ///Gets the measurements vector size
        virtual unsigned getMeasurementSize() const;


    protected:
        virtual Vector computeNoiselessMeasurement_();

        static const int stateSize_= 10;

        static const int measurementSize_=6;

    };

}

#endif //SIMULATIONACCELEROMETERGYROMETERSENSORHPP
