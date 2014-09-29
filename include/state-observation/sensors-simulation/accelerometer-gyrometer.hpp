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
        AccelerometerGyrometer();

        ///Virtual destructor
        virtual ~AccelerometerGyrometer(){}

        ///Gets the state vector Size
        virtual unsigned getStateSize() const;

        ///Gets the measurements vector size
        virtual unsigned getMeasurementSize() const;

        void setMatrixMode(bool matrixMode);


    protected:
        virtual Vector computeNoiselessMeasurement_();

        Matrix3 r_;
        Vector3 acc_;
        Vector3 omega_;
        Vector output_;

        bool matrixMode_;

        static const int stateSize_= 10;
        static const int stateSizeMatrix_= 15;

        static const int measurementSize_=6;

    };

}

#endif //SIMULATIONACCELEROMETERGYROMETERSENSORHPP
