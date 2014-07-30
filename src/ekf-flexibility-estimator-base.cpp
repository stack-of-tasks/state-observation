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

    void EKFFlexibilityEstimatorBase::setFlexibilityCovariance
                                                            (const Matrix & P)
    {
        ekf_.setStateCovariance(P);
    }

    Matrix EKFFlexibilityEstimatorBase::getFlexibilityCovariance
                                                            () const
    {
        return ekf_.getStateCovariance();
    }

    void EKFFlexibilityEstimatorBase::setMeasurement(const Vector & y)
    {
        ekf_.setMeasurement(y,k_+1);
    }

    void EKFFlexibilityEstimatorBase::setProcessNoiseCovariance(const Matrix & Q)
    {
        ekf_.setQ(Q);
    }

    void EKFFlexibilityEstimatorBase::setMeasurementNoiseCovariance(const Matrix & R)
    {
        ekf_.setR(R);
    }

    ///gets the covariance matrices for the process noises
    Matrix EKFFlexibilityEstimatorBase::getProcessNoiseCovariance() const
    {
        return ekf_.getQ();
    }

    ///gets the covariance matrices for the sensor noises
    Matrix EKFFlexibilityEstimatorBase::getMeasurementNoiseCovariance() const
    {
        return ekf_.getR();
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
            std::cout << "k " << k_ << std::endl;

            unsigned i;
            for (i=ekf_.getCurrentTime()+1; i<=k_; ++i)
            {
                    ekf_.setA(ekf_.getAMatrixFD(dx_));
                    ekf_.setC(ekf_.getCMatrixFD(dx_));
                    ekf_.getEstimatedState(i);
            }
            Vector x(ekf_.getEstimatedState(k_));

            if (x==x)//detect NaN values
            {
                lastX_=x;
            }
            else //delete NaN values
            {
                ekf_.setState(lastX_,k_);
                std::cout << "\n\n\n\n\n Need to reset covariance matrix \n\n\n\n\n" << std::endl;
                resetCovarianceMatrices();
            }
        }
        return lastX_;
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

    Vector EKFFlexibilityEstimatorBase::getPrediction()
    {
        return ekf_.getPrediction();
    }

}
}
