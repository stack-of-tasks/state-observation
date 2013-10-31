#ifndef FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_IMU_H
#define FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_IMU_H

#include <state-observation/flexibility-estimation/ekf-flexibility-estimator-base.hpp>
#include <state-observation/flexibility-estimation/imu-fixed-contact-dynamical-system.hpp>

namespace stateObservation
{
namespace flexibilityEstimation
{

    class FixedContactEKFFlexEstimatorIMU :
                                            public EKFFlexibilityEstimatorBase,
                                            private boost::noncopyable
    {
    public:
        FixedContactEKFFlexEstimatorIMU( double d=0.005 );

        virtual ~FixedContactEKFFlexEstimatorIMU();

        void setContactsNumber(unsigned i);

        void setContactPosition(unsigned i, Vector3 position);

        virtual void setMeasurement(const Vector & y);

        virtual void setVirtualMeasurementsCovariance(double c_);

        virtual void setNoiseCovariances(const Matrix & Q, const Matrix & R);

        virtual void setFlexibilityGuess(const Matrix & x);

        virtual Matrix4 getFlexibility();

        virtual unsigned getMeasurementSize() const ;

        virtual unsigned getStateSize() const ;

        virtual unsigned getInputSize() const ;


    protected:



        virtual void updateCovarianceMatrix_();

        IMUFixedContactDynamicalSystem functor_;

        double virtualMeasurementCovariance_;

        Matrix R_,Q_;

        static const unsigned stateSizeConst_=18;
        static const unsigned measurementSizeConst_=6;
        static const unsigned inputSizeConst_=15;

    private:
    };

}
}
#endif // FLEXBILITYESTMATOR_FIXEDCONTACTEKFFLEXIBILITYESTIMATOR_H
