#ifndef FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
#define FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H

#include <boost/utility.hpp>

#include <state-observation/observer/extended-kalman-filter.hpp>

#include <state-observation/flexibility-estimation/flexibility-estimator-base.hpp>


namespace stateObservation
{
namespace flexibilityEstimation
{
    class EKFFlexibilityEstimatorBase:
                public FlexibilityEstimatorBase
    {
    public:

        EKFFlexibilityEstimatorBase(unsigned stateSize,
                                    unsigned measurementSize,
                                    unsigned inputSize,
                                    const Vector & dx,
                                    double d=0.005);

        virtual ~EKFFlexibilityEstimatorBase();

        virtual void setFlexibilityGuess(const Matrix &)=0;

        virtual void setFlexibilityGuessCovariance(const Matrix & P);

        virtual void setNoiseCovariances(const Matrix & Q, const Matrix & R);

        virtual void setMeasurement(const Vector & y);

        virtual void setInput(const Vector & u);

        virtual void setMeasurementInput(const Vector & u);

        virtual Vector getFlexibilityVector();

        virtual Matrix4 getFlexibility()=0;

        virtual const stateObservation::ExtendedKalmanFilter & getEKF() const;

        virtual stateObservation::ExtendedKalmanFilter & getEKF();

        virtual unsigned getStateSize() const =0;

        virtual unsigned getMeasurementSize() const =0;

        virtual unsigned getInputSize() const =0;


    protected:
        virtual void setJacobians(const Matrix & A, const Matrix & C);

        virtual void useFiniteDifferencesJacobians(Vector dx);

        stateObservation::ExtendedKalmanFilter ekf_;

        bool finiteDifferencesJacobians_;

        Vector dx_;

        unsigned k_;

    private:

    };
}
}
#endif // FLEXIBILITYESTIMATION_EKFFLEXIBILITYESTIMATORBASE_H
