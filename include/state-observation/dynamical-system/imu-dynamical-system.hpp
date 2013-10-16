/**
 * \file      imu-dynamical-system.hpp
 * \author    Mehdi Benallegue
 * \date       2013
 * \brief   The file describes the dynamical system defined by an inertial
 *          measurement unit (IMU) fixed on a rigid body.
 *
 * \details
 *
 *
 */

#ifndef IMU-DYNAMICAL-SYSTEM_HPP
#define IMU-DYNAMICAL-SYSTEM_HPP

#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/noise-base.hpp>

namespace stateObservation
{

     /**
    * \class  IMUDynamicalSystem
    * \brief  The class is an implementation of the dynamical system defined by
    *         an inertial measurement unit (IMU) fixed on a rigid body. The state
    *         is the position velocity and acceleration and the orientaion and rotation
    *         velocity and acceleration. The sensors are the accelerometer and the gyrometer
    *
    *
    */
    class IMUDynamicalSystem : public DynamicalSystemFunctorBase,
        private algorithm::RigidBodyKinematics
    {
    public:
        ///The constructor
        IMUDynamicalSystem();

        ///The virtual destructor
        virtual ~IMUDynamicalSystem();

        ///Description of the state dynamics
        virtual Vector stateDynamics
        (const Vector& x, const Vector& u, unsigned k);

        ///Description of the sensor's dynamics
        virtual Vector measureDynamics
        (const Vector& x, const Vector& u, unsigned k);

        ///Sets a noise which disturbs the state dynamics
        virtual void setProcessNoise( NoiseBase * );
        ///Removes the process noise
        virtual void resetProcessNoise();
        ///Gets the process noise
        virtual NoiseBase * getProcessNoise() const;

        ///Sets a noise which disturbs the measurements
        virtual void setMeasurementNoise( NoiseBase * );
        ///Removes the measurement noise
        virtual void resetMeasurementNoise();
        ///Gets a pointer on the measurement noise
        virtual NoiseBase * getMeasurementNoise() const;

        ///Set the period of the time discretization
        virtual void setSamplingPeriod(double dt);

        ///Gets the state size
        virtual unsigned getStateSize();
        ///Gets the input size
        virtual unsigned getInputSize();
        ///Gets the measurement size
        virtual unsigned getMeasurementSize();

    protected:
        AccelerometerGyrometer sensor_;

        NoiseBase * processNoise_;

        double dt_;

        Vector3 orientationVector_;
        Quaternion quaternion_;

        Quaternion computeQuaternion_(const Vector3 & x);

        static const unsigned stateSize_=18;
        static const unsigned inputSize_=6;
        static const unsigned measurementSize_=6;


    private:
    };
}
#endif // IMU-DYNAMICAL-SYSTEM_HPP
