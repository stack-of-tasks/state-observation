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
#include <state-observation/tools/hrp2.hpp>


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
    class 	IMUElasticLocalFrameDynamicalSystem :
        public stateObservation::DynamicalSystemFunctorBase,
        private stateObservation::algorithm::RigidBodyKinematics
    {
    public:
        struct input
        {
        ///indexes of the different components of a vector of the input state
          static const unsigned posCom = 0;
          static const unsigned velCom = 3;
          static const unsigned accCom = 6;
          static const unsigned inertia = 9;
          static const unsigned angMoment = 15;
          static const unsigned dotInertia = 18;
          static const unsigned dotAngMoment = 24;
          static const unsigned posIMU = 27;
          static const unsigned oriIMU = 30;
          static const unsigned linVelIMU = 33;
          static const unsigned angVelIMU = 36;
          static const unsigned linAccIMU = 39;
          static const unsigned contacts = 42;

        };

        ///constructor
        explicit IMUElasticLocalFrameDynamicalSystem(double dt);

        ///virtual destructor
        virtual ~IMUElasticLocalFrameDynamicalSystem();

        void test();

        // computation of the acceleration linear
        virtual Vector3 computeAccelerations
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& Inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              const Vector3& position, const Vector3& linVelocity,
              Vector3& linearAcceleration,  const Vector3 &oriVector ,
              const Matrix3& orientation, const Vector3& angularVel,
              Vector3& angularAcceleration);

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

        ///Sets the input size
        virtual void setInputSize(unsigned i);

        ///Gets the contact number
        ///virtual

        ///Gets the contacts position

        ///Gets the measurement size
        virtual unsigned getMeasurementSize();

        ///Sets the number of contacts
        virtual void setContactsNumber(unsigned);

        ///Sets the position of the contact number i
        virtual void setContactPosition(unsigned i, const Vector3 & position);

        ///Gets the position of the contact number i
        virtual Vector3 getContactPosition(unsigned i);

        ///Gets the nimber of contacts
        unsigned getContactsNumber(void) const;

        virtual void setContactModelNumber(unsigned nb);

        virtual void computeElastFeetContactForcesAndMoments
                              (const IndexedMatrixArray& contactPosArray,
                               const IndexedMatrixArray& contactOriArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector3& forces, Vector3& moments);

        virtual void computeElastPendulumForcesAndMoments
                              (const IndexedMatrixArray& PrArray,
                               const IndexedMatrixArray& PeArray,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector3& forces, Vector3& moments);

        virtual void computeForcesAndMoments
                              (const IndexedMatrixArray& position1,
                               const IndexedMatrixArray& position2,
                               const Vector3& position, const Vector3& linVelocity,
                               const Vector3& oriVector, const Matrix3& orientation,
                               const Vector3& angVel,
                               Vector3& forces, Vector3& moments);

        virtual Vector getForcesAndMoments();

        virtual void iterateDynamicsEuler
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& Inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              Vector3& position, Vector3& linVelocity, Vector3& linearAcceleration,
              Vector3 &oriVector, Vector3& angularVel, Vector3& angularAcceleration,
              double dt
              );

        virtual void iterateDynamicsRK4
             (const Vector3& positionCom, const Vector3& velocityCom,
              const Vector3& accelerationCom, const Vector3& AngMomentum,
              const Vector3& dotAngMomentum,
              const Matrix3& Inertia, const Matrix3& dotInertia,
              const IndexedMatrixArray& contactPos,
              const IndexedMatrixArray& contactOri,
              Vector3& position, Vector3& linVelocity, Vector3& linearAcceleration,
              Vector3 &oriVector, Vector3& angularVel, Vector3& angularAcceleration,
              double dt
              );


        virtual void setKfe(const Matrix3 & m);
        virtual void setKfv(const Matrix3 & m);
        virtual void setKte(const Matrix3 & m);
        virtual void setKtv(const Matrix3 & m);

        virtual void setRobotMass(double d);

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    protected:

        Matrix3 curRotation0;
        Vector3 orientationVector0;
        Matrix3 curRotation1;
        Vector3 orientationVector1;
        Matrix3 curRotation2;
        Vector3 orientationVector2;
        Matrix3 curRotation3;
        Vector3 orientationVector3;

        inline Vector3& orientationVector(int i)
        {
            if (i==0)
              return orientationVector0;
            if (i==1)
              return orientationVector1;
            if (i==2)
              return orientationVector2;

            return orientationVector3;
        }

        inline Matrix3& curRotation(int i)
        {
            if (i==0)
              return curRotation0;
            if (i==1)
              return curRotation1;
            if (i==2)
              return curRotation2;

            return curRotation3;
        }

        stateObservation::AccelerometerGyrometer sensor_;

        stateObservation::NoiseBase * processNoise_;

        double dt_;

        double robotMass_;
        double robotMassInv_;

        Matrix3& computeRotation_(const Vector3 & x, int i);

        static const unsigned stateSize_=18;
        static const unsigned inputSizeBase_=42;
        unsigned inputSize_;
        static const unsigned measurementSizeBase_=6;

        unsigned nbContacts_;
        unsigned contactModel_;

        Vector fc_;
        Vector tc_;

        unsigned measurementSize_;

        std::vector <Vector3,Eigen::aligned_allocator<Vector3> > contactPositions_;


        Matrix3 Kfe_, Kte_, Kfv_, Ktv_;



        unsigned kcurrent_;


        struct Optimization
        {


          Vector3 positionFlex;
          Vector3 velocityFlex;
          Vector3 accelerationFlex;
          Vector3 orientationFlexV;
          Vector3 angularVelocityFlex;
          Vector3 angularAccelerationFlex;

          Matrix3 rFlex;

          Vector xk1;

          Vector3 positionCom;
          Vector3 velocityCom;
          Vector3 accelerationCom;
          Vector3 AngMomentum;
          Vector3 dotAngMomentum;

          Vector3 positionControl;
          Vector3 velocityControl;
          Vector3 accelerationControl;
          Vector3 orientationControlV;
          Vector3 angularVelocityControl;

          Matrix3 rControl;

          Matrix3 rimu;

          IndexedMatrixArray contactPosV;
          IndexedMatrixArray contactOriV;


          Matrix3 inertia;
          Matrix3 dotInertia;

          EIGEN_MAKE_ALIGNED_OPERATOR_NEW

          //optimization of orientation transformation between vector3 to rotation matrix



        } op_;

    private:


    public:
    };
}
}



#endif /* DYNAMICAL_SYSTEM_HPP_ */
