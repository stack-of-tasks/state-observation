/**
 * \file      ekf-flexibility-estimator-base.hpp
 * \author    Mehdi Benallegue
 * \date      2013
 * \brief     Declare the class of the flexibility estimation using the extended
 *            Kalman Filter.
 *
 * \details
 *
 *
 */

#ifndef FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
#define FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H

#include <boost/utility.hpp>

#include <state-observation/observer/extended-kalman-filter.hpp>

#include <state-observation/flexibility-estimation/flexibility-estimator-base.hpp>


namespace stateObservation
{
namespace flexibilityEstimation
{
   /**
    * \class  EKFFlexibilityEstimatorBase
    * \brief  This class is the base class of the flexibility estimators that
    *         use an extended Kalman Filter. Several methods require to be overloaded
    *         to derive an implementation from this base class.
    *
    */

    class EKFFlexibilityEstimatorBase:
                public FlexibilityEstimatorBase
    {
    public:


        /// The constructor.
        ///  \li stateSize : size of the state vector
        ///  \li measurementSize : size of the measurements vector
        ///  \li inputSize : size of the input vector
        ///  \li dx gives the derivation step for a finite differences derivation method

        EKFFlexibilityEstimatorBase(unsigned stateSize,
                                    unsigned measurementSize,
                                    unsigned inputSize,
                                    const Vector & dx );


        ///virtual destructor
        virtual ~EKFFlexibilityEstimatorBase();

        ///Sets a value of the flexibility x_k provided from another source
        /// can be used for initialization of the estimator
        /// This is a pure virtual function that requires to be overloaded in
        /// implementation
        virtual void setFlexibilityGuess(const Matrix &)=0;

        ///Sets the covariance matrix of the flexibility Guess
        virtual void setFlexibilityGuessCovariance(const Matrix & P);

        ///Sets the covariance matrices for the process and the sensor noises
        /// \li Q process noise
        /// \li R sensor noise
        virtual void setNoiseCovariances(const Matrix & Q, const Matrix & R);

        ///gets the covariance matrices for the process noises
        virtual Matrix getProcessNoiseCovariance() const ;

        ///gets the covariance matrices for the sensor noises
        virtual Matrix getMeasurementNoiseCovariance() const ;

        /// Sets the value of the next sensor measurement y_{k+1}
        virtual void setMeasurement(const Vector & y);

        /// Sets the value of the next input for the state process dynamics
        /// i.e. : gives u_k such that x_{k+1} = f(x_k,u_k)
        virtual void setInput(const Vector & u);

        /// Sets the value of the next  measurement
        /// i.e. : gives u_{k+1} such that y_{k+1}=h(x_{k+1},u_{k+1})
        virtual void setMeasurementInput(const Vector & u);

        /// Gets an estimation of the flexibility in the form of a state vector \hat{x_{k+1}}
        virtual Vector getFlexibilityVector();

        /// Gets an estimation of the flexibility in the form of a homogeneous matrix
        virtual Matrix4 getFlexibility()=0;

        /// Gets a const reference on the extended Kalman filter
        virtual const stateObservation::ExtendedKalmanFilter & getEKF() const;

        /// Gets a reference on the extended Kalman filter
        virtual stateObservation::ExtendedKalmanFilter & getEKF();

        /// Gets the state size
        /// this method is pure virtual and reauires to be overloaded in implementation
        virtual unsigned getStateSize() const =0;

        /// Gets the measurements size
        /// this method is pure virtual and reauires to be overloaded in implementation
        virtual unsigned getMeasurementSize() const =0;

        /// Gets the input size
        /// this method is pure virtual and reauires to be overloaded in implementation
        virtual unsigned getInputSize() const =0;

        /// Gets a simulation of the
        virtual Vector getSimulatedMeasurement();

        ///Resets the covariance matrices to their original values
        virtual void resetCovarianceMatrices()=0;

        ///Get the last vector of inovation of the Kalman filter
        virtual Vector getInovation();

        ///Get the simulated measurement of the predicted state
        virtual Vector getPredictedMeaurement();


    protected:
        virtual void setJacobians(const Matrix & A, const Matrix & C);

        virtual void useFiniteDifferencesJacobians(Vector dx);

        stateObservation::ExtendedKalmanFilter ekf_;

        bool finiteDifferencesJacobians_;

        Vector dx_;

        Vector lastX_;

        unsigned k_;

    private:

    };
}
}
#endif // FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
