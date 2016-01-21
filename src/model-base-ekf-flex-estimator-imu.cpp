#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

#include <iostream>

const double initialVirtualMeasurementCovariance=1.e-10;

const double dxFactor = 1.0e-8;

namespace stateObservation
{
namespace flexibilityEstimation
{
    ModelBaseEKFFlexEstimatorIMU::ModelBaseEKFFlexEstimatorIMU(double dt):
        EKFFlexibilityEstimatorBase
            (stateSizeBase_,measurementSizeBase_,inputSizeBase_,
                            Matrix::Constant(getStateSize(),1,dxFactor)),
        functor_(dt),
        forceVariance_(1e-6),
        stateSize_(stateSizeBase_)
    {
        ekf_.setMeasureSize(functor_.getMeasurementSize());

        withComBias_=false;
        ekf_.setStateSize(stateSizeBase_);
        ModelBaseEKFFlexEstimatorIMU::resetCovarianceMatrices();

        Vector dx( Matrix::Constant(getStateSize(),1,dxFactor));//thanks Justin

        dx.segment(kine::ori,3).fill(1e-4) ;
        dx.segment(kine::angVel,3).fill(1e-4) ;

        ModelBaseEKFFlexEstimatorIMU::useFiniteDifferencesJacobians(dx);

        Vector x0(ekf_.stateVectorZero());
        x0(2)=-0.010835;

        lastX_=x0;
        ekf_.setState(x0,0);

        ekf_.setStateCovariance(Q_);

        ekf_.setFunctor(& functor_);

        on_=true;

        useFTSensors_= false;
    }



    ModelBaseEKFFlexEstimatorIMU::~ModelBaseEKFFlexEstimatorIMU()
    {
        //dtor
    }

    void ModelBaseEKFFlexEstimatorIMU::resetCovarianceMatrices()
    {

//            std::cout << "\n\n\n ============> RESET COVARIANCE MATRIX <=============== \n\n\n" << std::endl;

            R_=Matrix::Identity(getMeasurementSize(),getMeasurementSize());
            R_.block(0,0,3,3)=Matrix3::Identity()*1.e-6;//accelerometer
            R_.block(3,3,3,3)=Matrix3::Identity()*1.e-6;//gyrometer

            updateCovarianceMatrix_();

            Q_=ekf_.getQmatrixIdentity();
            Q_=Q_*1.e-8;
            Q_.block(kine::linVel,kine::linVel,3,3)=Matrix3::Identity()*1.e-8;
            Q_.block(kine::angVel,kine::angVel,3,3)=Matrix3::Identity()*1.e-8;
            Q_.block(kine::linAcc,kine::linAcc,3,3)=Matrix3::Identity()*1.e-4;
            Q_.block(kine::angAcc,kine::angAcc,3,3)=Matrix3::Identity()*1.e-4;
            if(withComBias_)  Q_.block(kine::comBias,kine::comBias,2,2)=Matrix3::Identity().block(0,0,2,2)*2.5e-8;

            ekf_.setQ(Q_);

            Matrix P0 (ekf_.getQmatrixIdentity());
            P0=P0*1e-2;
            P0.block(kine::linVel,kine::linVel,3,3)=Matrix3::Identity()*1.e-2;
            P0.block(kine::angVel,kine::angVel,3,3)=Matrix3::Identity()*1.e-2;
            P0.block(kine::linAcc,kine::linAcc,3,3)=Matrix3::Identity()*1.e-2;
            P0.block(kine::angAcc,kine::angAcc,3,3)=Matrix3::Identity()*1.e-2;
            if(withComBias_) P0.block(kine::comBias,kine::comBias,2,2)=Matrix3::Identity().block(0,0,2,2)*1.e-2;

            ekf_.setStateCovariance(P0);

    }

    void ModelBaseEKFFlexEstimatorIMU::setContactsNumber(unsigned i)
    {
        finiteDifferencesJacobians_=true;
        functor_.setContactsNumber(i);

        unsigned usize = functor_.getInputSize();

        if (inputSize_!=usize)
        {
          stateObservation::Vector saveu, newu;
          int v, vmin;

          saveu=ekf_.getInput(ekf_.getInputTime());

          newu=Vector::Zero(usize,1);

          vmin=std::min(saveu.size(),newu.size());
          for(v=0;v<vmin;++v)
          {
              newu(v)=saveu(v);
          }

          inputSize_=usize;
          ekf_.setInputSize(usize);

          setInput(newu);
        }


        ekf_.setInputSize(usize);

        if (useFTSensors_)
        {
          ekf_.setMeasureSize(functor_.getMeasurementSize());
          updateCovarianceMatrix_();
        }
    }

    void ModelBaseEKFFlexEstimatorIMU::setContactModel(unsigned nb)
    {
        functor_.setContactModel(nb);
    }


    void ModelBaseEKFFlexEstimatorIMU::setMeasurement(const Vector & y)
    {
        BOOST_ASSERT((getMeasurementSize()==unsigned(y.size())) &&
                "ERROR: The measurement vector has incorrect size");


        ekf_.setMeasurement(y,k_+1);

    }



    void ModelBaseEKFFlexEstimatorIMU::setFlexibilityGuess(const Matrix & x)
    {
        bool bstate =ekf_.checkStateVector(x);
        bool b6= (x.rows()==6 && x.cols()==1);
        bool bhomogeneous = (x.rows()==4 && x.cols()==4);

        (void)b6;//avoid warning

        BOOST_ASSERT((bstate||b6||bhomogeneous) &&
                "ERROR: The flexibility state has incorrect size \
                    must be 20x1 vector, 6x1 vector or 4x4 matrix");

        Vector x0 (x);

        if (bstate)
        {
            ekf_.setState(x0,k_);
        }
        else
        {
            if (bhomogeneous)
                x0=kine::homogeneousMatrixToVector6(x);

            Vector x_s = ekf_.stateVectorZero();

            x_s.segment(kine::pos,3)=x0.head<3> ();

            x_s.segment(kine::ori,3)=x0.tail<3> ();

            ekf_.setState(x_s,k_);

            ekf_.setQ(Q_);
        }
    }

    void ModelBaseEKFFlexEstimatorIMU::setMeasurementNoiseCovariance
                                            (const Matrix & R)
    {
        BOOST_ASSERT(unsigned(R.rows())==6 &&
                     unsigned(R.cols())==6 &&
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

    Vector ModelBaseEKFFlexEstimatorIMU::getForcesAndMoments()
    {
        return functor_.getForcesAndMoments(getFlexibilityVector(),ekf_.getInput(ekf_.getInputTime()));
    }

    void ModelBaseEKFFlexEstimatorIMU::updateCovarianceMatrix_()
    {

        if (useFTSensors_)
        {
          Matrix Rtemp(R_);

          R_=Matrix::Identity(getMeasurementSize(),getMeasurementSize());

          R_.block(0,0,6,6)=Rtemp.block(0,0,6,6);

          R_.block(6,6,functor_.getContactsNumber()*6,functor_.getContactsNumber()*6) =
            Matrix::Identity(functor_.getContactsNumber()*6,functor_.getContactsNumber()*6)*forceVariance_;
        }

        ekf_.setR(R_);
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getStateSize() const
    {
        return ekf_.getStateSize();
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getInputSize() const
    {
        return inputSize_;
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getMeasurementSize() const
    {
        return functor_.getMeasurementSize();
    }

    Matrix4 ModelBaseEKFFlexEstimatorIMU::getFlexibility()
    {

        Vector v (getFlexibilityVector());

        Vector6 v2;
        v2.head(3) = v.segment(kine::pos,3);
        v2.tail(3) = v.segment(kine::ori,3);

        return kine::vector6ToHomogeneousMatrix(v2);

    }

    timespec diff(const timespec & start, const timespec & end)
    {
            timespec temp;
            if ((end.tv_nsec-start.tv_nsec)<0) {
                    temp.tv_sec = end.tv_sec-start.tv_sec-1;
                    temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
            } else {
                    temp.tv_sec = end.tv_sec-start.tv_sec;
                    temp.tv_nsec = end.tv_nsec-start.tv_nsec;
            }
            return temp;
    }

    const Vector & ModelBaseEKFFlexEstimatorIMU::getFlexibilityVector()
    {
        timespec time1, time2, time3;

        if(on_==true & functor_.getContactsNumber() > 0)
        {
            if (ekf_.getMeasurementsNumber()>0)
            {
                clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
                clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
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
                    x_=ekf_.getEstimatedState(k_);

                    if (x_==x_)//detect NaN values
                    {
                        lastX_=x_;

                        ///regulate the part of orientation vector in the state vector
                        lastX_.segment(kine::ori,3)=kine::regulateOrientationVector(lastX_.segment(kine::ori,3));
                        ekf_.setState(lastX_,ekf_.getCurrentTime());
                    }
                    else //delete NaN values
                    {
                        ekf_.setState(lastX_,k_);

                        if(k_>1) //the first iteration give always nan without initialization
                        {
                            resetCovarianceMatrices();
                        }
                    }
                }
                clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);

                computeFlexibilityTime_=(double)diff(time2,time3).tv_nsec-(double)diff(time1,time2).tv_nsec;
            }
        }
        else
        {
            lastX_.setZero();
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

    void ModelBaseEKFFlexEstimatorIMU::setKfe(const Matrix3 & m)
    {
        functor_.setKfe(m);
    }

    void ModelBaseEKFFlexEstimatorIMU::setKfv(const Matrix3 & m)
    {
        functor_.setKfv(m);
    }

    void ModelBaseEKFFlexEstimatorIMU::setKte(const Matrix3 & m)
    {
        functor_.setKte(m);
    }

    void ModelBaseEKFFlexEstimatorIMU::setKtv(const Matrix3 & m)
    {
        functor_.setKtv(m);
    }

    double& ModelBaseEKFFlexEstimatorIMU::getComputeFlexibilityTime()
    {
        return computeFlexibilityTime_;
    }

    void ModelBaseEKFFlexEstimatorIMU::setWithForcesMeasurements(bool b)
    {
      if (useFTSensors_!= b)
      {
        functor_.setWithForceMeasurements(b);
        ekf_.setMeasureSize(functor_.getMeasurementSize());

        updateCovarianceMatrix_();
      }

      useFTSensors_=b;

    }

    void ModelBaseEKFFlexEstimatorIMU::setWithComBias(bool b)
    {
      if (withComBias_!= b)
      {
        withComBias_=b;
        functor_.setWithComBias(b);

        stateSize_=stateSizeBase_+(int)b*2;
        ekf_.setStateSize(stateSize_);

        Vector dx( Matrix::Constant(getStateSize(),1,dxFactor));//thanks Justin
        dx.segment(kine::ori,3).fill(1e-4) ;
        dx.segment(kine::angVel,3).fill(1e-4) ;
        dx_= dx;

        stateObservation::Vector x(lastX_); x=lastX_;
        lastX_.resize(stateSize_);
        if(b==true){
            lastX_.segment(0,x.size())=x;
        } else {
            lastX_=x.segment(0,lastX_.size());
        }

        ekf_.setState(lastX_,k_);
        resetCovarianceMatrices();
      }
    }

    void ModelBaseEKFFlexEstimatorIMU::setForceVariance(double d)
    {
        forceVariance_ = d;
        updateCovarianceMatrix_();
    }

}
}
