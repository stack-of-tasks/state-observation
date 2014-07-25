#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

const double initialVirtualMeasurementCovariance=1.e-10;

const double dxFactor = 1.0e-8;

namespace stateObservation
{
namespace flexibilityEstimation
{
    ModelBaseEKFFlexEstimatorIMU::ModelBaseEKFFlexEstimatorIMU(double dt):
        EKFFlexibilityEstimatorBase
            (stateSizeConst_,measurementSizeConst_,inputSizeBase_,
                            Matrix::Constant(getStateSize(),1,dxFactor)),
        virtualMeasurementCovariance_(initialVirtualMeasurementCovariance),
        functor_(dt)
    {

        //ekf_.setDirectInputStateFeedthrough(false);

        ObserverBase::InputVector u;
        u.resize(42);
        u <<    0.0145673,
                0.00153601,
                0.807688,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                9.59691,
                8.39051,
                1.73881,
                -0.00189139,
                0.146455,
                -0.0378545,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                1.1174,
                0.0,
                0.0,
                0.0,
                0.227298,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                22.7298,
                0.0,
                0.0,
                0.0;

        //setInput(u);

        ekf_.setMeasureSize(functor_.getMeasurementSize());

        ModelBaseEKFFlexEstimatorIMU::resetCovarianceMatrices();

        Vector x0=ekf_.stateVectorZero();
        x0(2)=-0.010835;

        lastX_=x0;

        ekf_.setState(x0,0);

        ekf_.setStateCovariance(Q_);

        ekf_.setFunctor(& functor_);

        on_=true;

    }



    ModelBaseEKFFlexEstimatorIMU::~ModelBaseEKFFlexEstimatorIMU()
    {
        //dtor
    }

    void ModelBaseEKFFlexEstimatorIMU::resetCovarianceMatrices()
    {
        std::cout << "\n\n\n ============> RESET COVARIANCE MATRIX <=============== \n\n\n" << std::endl;

        R_=Matrix::Identity(getMeasurementSize(),getMeasurementSize());
        R_.block(0,0,3,3)=Matrix3::Identity()*1.e-6;//accelerometer
        R_.block(3,3,3,3)=Matrix3::Identity()*1.e-6;//gyrometer

        updateCovarianceMatrix_();

        Q_=ekf_.getQmatrixIdentity();
        Q_=Q_*1.e-8;
        Q_.block(kine::linVel,kine::linVel,3,3)=Matrix3::Identity()*1.e-4;
        Q_.block(kine::angVel,kine::angVel,3,3)=Matrix3::Identity()*1.e-4;
        Q_.block(kine::linAcc,kine::linAcc,3,3)=Matrix3::Identity()*1.e-2;
        Q_.block(kine::angAcc,kine::angAcc,3,3)=Matrix3::Identity()*1.e-2;

        ekf_.setQ(Q_);

        Matrix P0 (ekf_.getQmatrixIdentity());
        P0=P0*1e-2;
        P0.block(kine::linVel,kine::linVel,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(kine::angVel,kine::angVel,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(kine::linAcc,kine::linAcc,3,3)=Matrix3::Identity()*1.e-2;
        P0.block(kine::angAcc,kine::angAcc,3,3)=Matrix3::Identity()*1.e-2;

        ekf_.setStateCovariance(P0);

    }

    void ModelBaseEKFFlexEstimatorIMU::setContactsNumber(unsigned i)
    {
        finiteDifferencesJacobians_=true;
        functor_.setContactsNumber(i);
        updateCovarianceMatrix_();

    }


    void ModelBaseEKFFlexEstimatorIMU::setMeasurement(const Vector & y)
    {
        BOOST_ASSERT((getMeasurementSize()==y.size()) &&
                "ERROR: The measurement vector has incorrect size");


        Vector y2 = ekf_.measureVectorZero();
        y2.head(getMeasurementSize()) = y;
        ekf_.setMeasurement(y2,k_+1);

    }

    void ModelBaseEKFFlexEstimatorIMU::setVirtualMeasurementsCovariance
                                                                    (double c)
    {
        virtualMeasurementCovariance_=c;
        updateCovarianceMatrix_();
    }

    double ModelBaseEKFFlexEstimatorIMU::getVirtualMeasurementsCovariance() const
    {
        return virtualMeasurementCovariance_;
    }

    void ModelBaseEKFFlexEstimatorIMU::setFlexibilityGuess(const Matrix & x)
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
                x0=kine::homogeneousMatrixToVector6(x);

            Vector x_s = ekf_.stateVectorZero();

            x_s.segment(kine::pos,3)=x0.head(3);

            x_s.segment(kine::ori,3)=x0.tail(3);

            ekf_.setState(x_s,k_);

            ekf_.setQ(Q_);
        }
    }

    void ModelBaseEKFFlexEstimatorIMU::setMeasurementNoiseCovariance
                                            (const Matrix & R)
    {
        BOOST_ASSERT(R.rows()==getMeasurementSize() &&
                     R.cols()==getMeasurementSize() &&
                    "ERROR: The measurement noise covariance matrix R has \
                        incorrect size");

        R_=R;
        updateCovarianceMatrix_();
    }


    void ModelBaseEKFFlexEstimatorIMU::setProcessNoiseCovariance
                                            (const Matrix & Q)
    {
        Q_=Q;
        ekf_.setQ(Q_);
    }

    Matrix ModelBaseEKFFlexEstimatorIMU::getProcessNoiseCovariance() const
    {
        return Q_;
    }


    Matrix ModelBaseEKFFlexEstimatorIMU::getMeasurementNoiseCovariance() const
    {
        return R_;
    }

    void ModelBaseEKFFlexEstimatorIMU::updateCovarianceMatrix_()
    {
        Matrix R = ekf_.getRmatrixIdentity();
        R.block(0,0, getMeasurementSize() , getMeasurementSize())=R_;

        for (int i= getMeasurementSize() ; i<ekf_.getMeasureSize() ; ++i)
        {
            R(i,i)=virtualMeasurementCovariance_;
        }

        ekf_.setR(R);
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getStateSize() const
    {
        return stateSizeConst_;
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getInputSize() const
    {
        return inputSize_;
    }

    void ModelBaseEKFFlexEstimatorIMU::setInputSize(unsigned i)
    {

        stateObservation::Vector saveu, newu;
        unsigned saveInputSize;
        int v, vmax;

        saveInputSize=saveu.size();
        saveu=ekf_.getInput(ekf_.getInputTime());

        newu=Vector::Zero(i,1);

//        vmax=min(saveu.size(),newu.size());
//        for(v=0;v<vmax;++v)
//        {
//            newu(v)=saveu(v);
//        }

        inputSize_=i;
        ekf_.setInputSize(i);
        functor_.setInputSize(i);

        setInput(newu);

    }


    unsigned ModelBaseEKFFlexEstimatorIMU::getMeasurementSize() const
    {
        return measurementSizeConst_;
    }

    Matrix4 ModelBaseEKFFlexEstimatorIMU::getFlexibility()
    {

        Vector v (getFlexibilityVector());

        Vector6 v2;
        v2.head(3) = v.segment(kine::pos,3);
        v2.tail(3) = v.segment(kine::ori,3);

        return kine::vector6ToHomogeneousMatrix(v2);

    }

    Vector ModelBaseEKFFlexEstimatorIMU::getFlexibilityVector()
    {

        if (ekf_.getMeasurementsNumber()>0)
        {
            if(on_==true)
            {
                lastX_ =EKFFlexibilityEstimatorBase::getFlexibilityVector();
            }
            else
            {
               // vector v(getEKF().getStateSize());
                lastX_.setZero();
            }

            ///regulate the part of orientation vector in the state vector
            lastX_.segment(kine::ori,3)=kine::regulateOrientationVector(lastX_.segment(kine::ori,3));
            ekf_.setState(lastX_,ekf_.getCurrentTime());

        }

        return lastX_;
    }

    void ModelBaseEKFFlexEstimatorIMU::setSamplingPeriod(double dt)
    {
        dt_=dt;
    }

    /// Enable or disable the estimation
    void ModelBaseEKFFlexEstimatorIMU::setOn(bool & b)
    {
        on_=b;
    }


}
}
