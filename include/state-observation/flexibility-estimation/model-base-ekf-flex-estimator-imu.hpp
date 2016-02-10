/**
 * \file      model-base-ekf-flex-estimator-imu.hpp
 * \author    Mehdi Benallegue
 * \date      2013
 * \brief     Declares the class of the estimation of the flexibility using an
 *            extended Kalman filter and a fixed contact hypothesis
 *
 * \details
 *
 *
 */


#ifndef FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H
#define FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H

#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
//#include <state-observation/flexibility-estimation/stable-imu-fixed-contact-dynamical-system.hpp>
#include <state-observation/flexibility-estimation/imu-elastic-local-frame-dynamical-system.hpp>
//#include <state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

    /**
    * \class  ModelBaseEKFFlexEstimatorIMU
    * \brief  This class implements the flexibility estimation of a robot with
    *         the hypothesis that the contact positions do not move. This constraint
    *         is expressed using fictious measurements but the interface is transparent
    *         to this assumption, the state is expressed using classical representation
    *         of position, velocity, acceleration, orientation (using (theta x mu) representation)
    *         angular velocity (omega) and acceleration (omega dot)
    *
    */

    class ModelBaseEKFFlexEstimatorIMU :
                                            public EKFFlexibilityEstimatorBase,
                                            private boost::noncopyable
    {
    public:

        struct contactModel
        {
          ///indexes of the different components of a vector of the input state
          static const unsigned elasticContact= IMUElasticLocalFrameDynamicalSystem::contactModel::elasticContact;
          static const unsigned pendulum= IMUElasticLocalFrameDynamicalSystem::contactModel::pendulum;

        };

        ///The constructor, it requires the value of the time discretization period
        explicit ModelBaseEKFFlexEstimatorIMU( double dt=0.005 );

        ///Virtual destructor
        virtual ~ModelBaseEKFFlexEstimatorIMU();

        ///Sets the number of contacts can be changed online
        void setContactsNumber(unsigned i);

        void setContactModel(unsigned nb);

        /// Sets the value of the next sensor measurement y_{k+1}
        virtual void setMeasurement(const Vector & y);

        ///Sets the process covariance matrice
        virtual void setProcessNoiseCovariance(const Matrix & Q);

        ///Sets the measurements covariance matrice
        virtual void setMeasurementNoiseCovariance(const Matrix & R);

        ///gets the covariance matrices for the process noises
        virtual Matrix getProcessNoiseCovariance() const ;

        ///gets the covariance matrices for the sensor noises
        virtual Matrix getMeasurementNoiseCovariance() const ;

        virtual Vector getForcesAndMoments();

        // get state covariance
        stateObservation::Vector getStateCovariance() const
        {
            stateObservation::Matrix P(ekf_.getStateCovariance());
            stateObservation::Vector Pvec(ekf_.getStateSize());
            for(int i=0;i<ekf_.getStateSize();++i) Pvec(i)=P(i,i);
            return Pvec;
        }

        ///Sets a value of the flexibility x_k provided from another source
        /// can be used for initialization of the estimator
        virtual void setFlexibilityGuess(const Matrix & x);

        /// Gets an estimation of the flexibility in the form of a homogeneous matrix
        virtual Matrix4 getFlexibility();

        /// Gets an estimation of the flexibility in the form of a state vector \hat{x_{k+1}}
        virtual const Vector& getFlexibilityVector();

        virtual double& getComputeFlexibilityTime();


        virtual unsigned getMeasurementSize() const ;

        virtual unsigned getStateSize() const ;

        virtual unsigned getInputSize() const ;


        ///sets to whether or not the force mesurements are taken into account
        virtual void setWithForcesMeasurements(bool);

        virtual void setWithComBias(bool b);

        virtual bool getWithComBias()
        {
            return withComBias_;
        }

        virtual void setForceVariance(double d);

        /// sets the sampling period
        virtual void setSamplingPeriod(double);

        /// Enable or disable the estimation
        void setOn(bool & b);

        virtual void setKfe(const Matrix3 & m);
        virtual void setKfv(const Matrix3 & m);
        virtual void setKte(const Matrix3 & m);
        virtual void setKtv(const Matrix3 & m);

        ///Resets the covariance matrices to their original values
        virtual void resetCovarianceMatrices();

        virtual void setMass(double m);

        virtual void setAngularAccelerationLimit(const Vector3 & v);
        virtual void setLinearAccelerationLimit(const Vector3 & v);

    protected:

        virtual void updateCovarianceMatrix_();

        IMUElasticLocalFrameDynamicalSystem functor_;

        Vector x_;

        Matrix R_,Q_;

        static const unsigned measurementSizeBase_=42;
        static const unsigned inputSizeBase_=42;
        unsigned inputSize_;
        const unsigned stateSize_;

        double dt_;//sampling period
        bool on_;
        double computeFlexibilityTime_;

        double forceVariance_;//force sensor variance

        bool useFTSensors_;
        bool withComBias_;

        Vector3 limitAngularAcceleration_;
        Vector3 limitLinearAcceleration_;

    private:
    };

}
}
#endif // FLEXBILITYESTMATOR_MODELBASEEKFFLEXIBILITYESTIMATOR_IMU_H
