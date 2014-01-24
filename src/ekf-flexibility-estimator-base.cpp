#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
namespace stateObservation
{
namespace flexibilityEstimation
{
    EKFFlexibilityEstimatorBase::EKFFlexibilityEstimatorBase(unsigned stateSize,
                                    unsigned measurementSize,
                                    unsigned inputSize,
                                    const Vector & dx):
        FlexibilityEstimatorBase(),
        ekf_(stateSize,measurementSize,inputSize),
        k_(0),
        finiteDifferencesJacobians_(true),
        dx_(dx)
    {
        //ctor
    }

    EKFFlexibilityEstimatorBase::~EKFFlexibilityEstimatorBase()
    {
        //dtor
    }

    void EKFFlexibilityEstimatorBase::setFlexibilityGuessCovariance
                                                            (const Matrix & P)
    {
        ekf_.setStateCovariance(P);
    }

    void EKFFlexibilityEstimatorBase::setMeasurement(const Vector & y)
    {
        ekf_.setMeasurement(y,k_+1);
    }

    void EKFFlexibilityEstimatorBase::setNoiseCovariances
                                     (const Matrix & Q, const Matrix & R)
    {
        ekf_.setQ(Q);
        ekf_.setR(R);
    }

    void EKFFlexibilityEstimatorBase::useFiniteDifferencesJacobians(Vector x)
    {
        finiteDifferencesJacobians_= true ;
        dx_ = x;
    }

    void EKFFlexibilityEstimatorBase::setJacobians
                                          (const Matrix & A, const Matrix & C)
    {
        finiteDifferencesJacobians_= false ;
        ekf_.setA(A);
        ekf_.setC(C);
    }

    void EKFFlexibilityEstimatorBase::setInput(const Vector & u)
    {
        ekf_.setInput(u,k_);
    }

    void EKFFlexibilityEstimatorBase::setMeasurementInput(const Vector & u)
    {
        ekf_.setInput(u,k_+1);
    }

    Vector EKFFlexibilityEstimatorBase::getFlexibilityVector()
    {
        if (ekf_.getMeasurementsNumber()>0)
        {
            k_=ekf_.getMeasurementTime();

            unsigned i;
            for (i=ekf_.getCurrentTime()+1; i<=k_; ++i)
            {
                    ekf_.setA(ekf_.getAMatrixFD(dx_));
                    ekf_.setC(ekf_.getCMatrixFD(dx_));

                    ekf_.getEstimatedState(i);
            }
        }

        Vector x(ekf_.getEstimatedState(k_));

        if (x==x)//detect NaN values
        {
            return lastX_=x;
        }
        else //delete NaN values
        {
            ekf_.setState(lastX_,k_);
            resetCovarianceMatrices();
            return lastX_;
        }
    }

    const stateObservation::ExtendedKalmanFilter &
        EKFFlexibilityEstimatorBase::getEKF() const
    {
        return ekf_;
    }

    stateObservation::ExtendedKalmanFilter & EKFFlexibilityEstimatorBase::getEKF()
    {
        return ekf_;
    }

    Vector EKFFlexibilityEstimatorBase::getSimulatedMeasurement()
    {

        getFlexibilityVector();

        return ekf_.getSimulatedMeasurement(ekf_.getCurrentTime());
    }

    Vector EKFFlexibilityEstimatorBase::getInovation()
    {
        return ekf_.getInovation();
    }

    Vector EKFFlexibilityEstimatorBase::getPredictedMeaurement()
    {
        return ekf_.getPredictedMeaurement();
    }

}
}
