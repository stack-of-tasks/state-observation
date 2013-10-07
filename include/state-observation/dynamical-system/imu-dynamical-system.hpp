#ifndef IMU-DYNAMICAL-SYSTEM_HPP
#define IMU-DYNAMICAL-SYSTEM_HPP

#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/noise/noise-base.hpp>

namespace stateObservation
{
    class IMUDynamicalSystem : public DynamicalSystemFunctorBase,
        private algorithm::RigidBodyKinematics
    {
    public:
        IMUDynamicalSystem();
        virtual ~IMUDynamicalSystem();

        virtual Vector stateDynamics
        (const Vector& x, const Vector& u, unsigned k);

        virtual Vector measureDynamics
        (const Vector& x, const Vector& u, unsigned k);

        virtual void setProcessNoise( NoiseBase * );
        virtual void resetProcessNoise();

        virtual void setMeasurementNoise( NoiseBase * );
        virtual void resetMeasurementNoise();

        virtual void setSamplingPeriod(double dt);

        virtual unsigned getStateSize();
        virtual unsigned getInputSize();
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
