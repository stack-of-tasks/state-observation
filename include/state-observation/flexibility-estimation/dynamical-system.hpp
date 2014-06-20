/*
 * dynamical-system.hpp
 *
 *  Created on: 19 mai 2014
 *      Author: alexis
 */

#ifndef DYNAMICAL_SYSTEM_HPP_
#define DYNAMICAL_SYSTEM_HPP_


#include <vector>

#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>
#include <state-observation/noise/noise-base.hpp>
#include <state-observation/sensors-simulation/accelerometer-gyrometer.hpp>
#include <state-observation/dynamical-system/algorithm/rigid-body-kinematics.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

    /**
    * \class  DynamicalSystem
    * \brief  This class describes the dynamics of a robot's flexibility
    *         this dynamics is the simplest possible system, the flexibility
    *         is expressed as a rotation against the contact positions with no
    *         other hypothesis than that the contact points are at constant position
    *
    */
    class 	DynamicalSystem :
        public stateObservation::DynamicalSystemFunctorBase,
        private stateObservation::algorithm::RigidBodyKinematics
    {
    public:
        ///constructor
        explicit DynamicalSystem(double dt);

        ///virtual destructor
        virtual ~DynamicalSystem();

        // stabilization of the acceleration linear
        virtual Vector3 stabilizeAccelerationLinear
    	(Vector3 , Vector3);

        // stabilization of the acceleration angular
        virtual Vector3 stabilizeAccelerationAngular
        (Vector3 , Vector3);

        virtual Vector3 computeFc(const stateObservation::Vector& x);

        virtual Vector3 computeTc(const stateObservation::Vector& x);



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

        Vector3Unaligned orientationVector_;
        QuaternionUnaligned quaternion_;

        Quaternion computeQuaternion_(const Vector3 & x);

        static const unsigned stateSize_=18;
        static const unsigned inputSize_=15;
        static const unsigned measurementSizeBase_=6;

        unsigned measurementSize_;

        std::vector <Vector3,Eigen::aligned_allocator<Vector3> >
            contactPositions_;

    private:

    public:
    };
}
}



#endif /* DYNAMICAL_SYSTEM_HPP_ */
