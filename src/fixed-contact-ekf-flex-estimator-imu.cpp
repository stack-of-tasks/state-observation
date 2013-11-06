#include <state-observation/flexibility-estimation/fixed-contact-ekf-flex-estimator-imu.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>


const double initialVirtualMeasurementCovariance=1.e-10;

namespace stateObservation
{
namespace flexibilityEstimation
{
    FixedContactEKFFlexEstimatorIMU::FixedContactEKFFlexEstimatorIMU(double dt):
        EKFFlexibilityEstimatorBase
            (stateSizeConst_,measurementSizeConst_,inputSizeConst_,
                            Matrix::Constant(getStateSize(),1,1.0e-8),dt),
        virtualMeasurementCovariance_(initialVirtualMeasurementCovariance)
    {
        ekf_.setDirectInputStateFeedthrough(false);

        ekf_.setMeasureSize(functor_.getMeasurementSize());

        R_=ekf_.getRmatrixIdentity();
        R_=R_*1.e-8;
        ekf_.setR(R_);

        Q_=ekf_.getQmatrixIdentity();
        Q_=Q_*1.e-4;
        Q_.block(6,6,3,3)=Matrix3::Identity()*1.e-2;
        Q_.block(15,15,3,3)=Matrix3::Identity()*1.e-2;
        ekf_.setQ(Q_);

        Vector x0=ekf_.stateVectorZero();
        ekf_.setState(x0,0);

        ekf_.setStateCovariance(Q_);

        ekf_.setFunctor(& functor_);
    }

    FixedContactEKFFlexEstimatorIMU::~FixedContactEKFFlexEstimatorIMU()
    {
        //dtor
    }

    void FixedContactEKFFlexEstimatorIMU::setContactsNumber(unsigned i)
    {
        finiteDifferencesJacobians_=true;
        functor_.setContactsNumber(i);
        ekf_.setMeasureSize(functor_.getMeasurementSize());

    }

    void FixedContactEKFFlexEstimatorIMU::setContactPosition
                                            (unsigned i, Vector3 position)
    {
        functor_.setContactPosition(i,position);
    }

    void FixedContactEKFFlexEstimatorIMU::setMeasurement(const Vector & y)
    {
        Vector y2 = ekf_.measureVectorZero();
        y2.head(getMeasurementSize()) = y;
        ekf_.setMeasurement(y2,k_+1);
        updateCovarianceMatrix_();
    }

    void FixedContactEKFFlexEstimatorIMU::setVirtualMeasurementsCovariance
                                                                    (double c_)
    {
        virtualMeasurementCovariance_=c_;
    }

    void FixedContactEKFFlexEstimatorIMU::setFlexibilityGuess(const Matrix & x)
    {
        bool bstate =ekf_.checkStateVector(x);
        bool b6= (x.rows()==6 && x.cols()==1);
        bool bhomogeneous = (x.rows()==4 && x.cols()==4);

        BOOST_ASSERT((bstate||b6||bhomogeneous) &&
                "ERROR: The flexibility state has incorrect size \
                    must be 18x1 vector, 6x1 vector or 4x4 matrix");

        Vector x0=x;

        if (bstate)
        {
            ekf_.setState(x0,k_);
        }
        else
        {
            if (bhomogeneous)
                x0=tools::homogeneousMatrixToVector6(x);

            Vector x_s = ekf_.stateVectorZero();

            x_s.head(3)=x0.head(3);

            x_s.segment(9,3)=x0.tail(3);

            ekf_.setState(x_s,k_);

            ekf_.setQ(Q_);
        }
    }

    void FixedContactEKFFlexEstimatorIMU::setNoiseCovariances
                                            (const Matrix & Q, const Matrix & R)
    {
        BOOST_ASSERT(R.rows()==getMeasurementSize() &&
                     R.cols()==getMeasurementSize() &&
                    "ERROR: The measurement noise covariance matrix Q has \
                        incorrect size");

        R_=R;
        updateCovarianceMatrix_();

        Q_=Q;
        ekf_.setQ(Q_);

    }

    void FixedContactEKFFlexEstimatorIMU::updateCovarianceMatrix_()
    {
        Matrix R = ekf_.getRmatrixIdentity();
        R.block(0,0, getMeasurementSize() , getMeasurementSize())=R_;

        for (int i= getMeasurementSize() ; i<ekf_.getMeasureSize() ; ++i)
        {
            R(i,i)=virtualMeasurementCovariance_;
        }

        ekf_.setR(R);
    }

    unsigned FixedContactEKFFlexEstimatorIMU::getStateSize() const
    {
        return stateSizeConst_;
    }

    unsigned FixedContactEKFFlexEstimatorIMU::getInputSize() const
    {
        return inputSizeConst_;
    }

    unsigned FixedContactEKFFlexEstimatorIMU::getMeasurementSize() const
    {
        return measurementSizeConst_;
    }

    Matrix4 FixedContactEKFFlexEstimatorIMU::getFlexibility()
    {
        Vector v (getFlexibilityVector());
        Vector v2 (Matrix::Zero(6,1));

        v2.head(3) = v.head(3);
        v2.tail(3) = v.segment(9,3);

        return tools::vector6ToHomogeneousMatrix(v2);

    }

    Vector FixedContactEKFFlexEstimatorIMU::getFlexibilityVector()
    {
        Vector v(EKFFlexibilityEstimatorBase::getFlexibilityVector());

        ///regulate the part of orientation vector in the state vector
        v.segment(9,3)=tools::regulateOrientationVector(v.segment(9,3));

        ekf_.setState(v,ekf_.getCurrentTime());

        return v;
    }

}
}
