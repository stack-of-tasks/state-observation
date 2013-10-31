#ifndef FIXED-CONTACTS-IMU-DYNAMICS-FUNCTOR_HPP
#define FIXED-CONTACTS-IMU-DYNAMICS-FUNCTOR_HPP

#include <vector>

#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/noise/noise-base.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

    class IMUFixedContactDynamicalSystem :
        public stateObservation::DynamicalSystemFunctorBase,
        private stateObservation::algorithm::RigidBodyKinematics
    {
    public:
        IMUFixedContactDynamicalSystem();
        virtual ~IMUFixedContactDynamicalSystem();

        ///Description of the state dynamics
        virtual stateObservation::Vector stateDynamics
        (const stateObservation::Vector& x, const stateObservation::Vector& u,
            unsigned k);

        ///Description of the sensor's dynamics
        virtual stateObservation::Vector measureDynamics
        (const stateObservation::Vector& x, const stateObservation::Vector& u,
            unsigned k);

        ///Sets a noise which disturbs the state dynamics
        virtual void setProcessNoise( stateObservation::NoiseBase * );

        ///Removes the process noise
        virtual void resetProcessNoise();

        ///Gets the process noise
        virtual stateObservation::NoiseBase * getProcessNoise() const;

        ///Sets a noise which disturbs the measurements
        virtual void setMeasurementNoise( stateObservation::NoiseBase * );

        ///Removes the measurement noise
        virtual void resetMeasurementNoise();

        ///Gets a pointer on the measurement noise
        virtual stateObservation::NoiseBase * getMeasurementNoise() const;

        ///Set the period of the time discretization
        virtual void setSamplingPeriod(double dt);

        ///Gets the state size
        virtual unsigned getStateSize();
        ///Gets the input size
        virtual unsigned getInputSize();
        ///Gets the measurement size
        virtual unsigned getMeasurementSize();

        ///Sets the number of contacts
        virtual void setContactsNumber(unsigned);

        ///Sets the position of the contact number i
        virtual void setContactPosition(unsigned i, const Vector3 & position);

    protected:

        stateObservation::AccelerometerGyrometer sensor_;

        stateObservation::NoiseBase * processNoise_;

        double dt_;

        Vector3 orientationVector_;
        Quaternion quaternion_;

        Quaternion computeQuaternion_(const Vector3 & x);

        static const unsigned stateSize_=18;
        static const unsigned inputSize_=15;
        static const unsigned measurementSizeBase_=6;

        unsigned measurementSize_;

        std::vector <Vector3,Eigen::aligned_allocator<Vector3> >
            contactPositions_;

    private:
    };
}
}

#endif // FIXED-CONTACTS-IMU-DYNAMICS-FUNCTOR_HPP
