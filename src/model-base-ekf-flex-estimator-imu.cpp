#include <state-observation/flexibility-estimation/model-base-ekf-flex-estimator-imu.hpp>
#include <state-observation/tools/miscellaneous-algorithms.hpp>

const double dxFactor = 1.0e-8;
const int stateSize = 35;

namespace stateObservation
{
namespace flexibilityEstimation
{
    typedef IMUElasticLocalFrameDynamicalSystem::state state;

    ModelBaseEKFFlexEstimatorIMU::ModelBaseEKFFlexEstimatorIMU(double dt):
        EKFFlexibilityEstimatorBase
            (stateSize,measurementSizeBase_,inputSizeBase_,
                            Matrix::Constant(stateSize,1,dxFactor)),
        functor_(dt),
        stateSize_(stateSize),
        unmodeledForceVariance_(1e-6),
        forceVariance_(1e-4),
        absPosVariance_(1e-4),
        useFTSensors_(false),
        withAbsolutePos_(false),
        withComBias_(false),
        withUnmodeledForces_(false),
        limitOn_(true)
    {
        ekf_.setMeasureSize(functor_.getMeasurementSize());
        ekf_.setStateSize(stateSize_);
        ekf_.setInputSize(functor_.getInputSize());
        inputSize_=functor_.getInputSize();

        ModelBaseEKFFlexEstimatorIMU::resetCovarianceMatrices();

        Vector dx(Matrix::Constant(getStateSize(),1,dxFactor));

        dx.segment(state::ori,3).fill(1e-4) ;
        dx.segment(state::angVel,3).fill(1e-4) ;

        ModelBaseEKFFlexEstimatorIMU::useFiniteDifferencesJacobians(dx);

        Vector x0(ekf_.stateVectorZero());

        lastX_=x0;
        ekf_.setState(x0,0);

        ekf_.setStateCovariance(Q_);

        ekf_.setFunctor(& functor_);

        functor_.setFDstep(dx_);

        on_=true;

        Vector3 v1, v2;
        Matrix3 v;
        v1 << 100,
              100,
              1000;
        v2 << 10,
              10,
              100;
        v.setIdentity();
        v*=cst::gravityConstant*hrp2::m;
        limitForces_=v*v1;
        limitTorques_=v*v2;

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
            stateObservation::Matrix m; m.resize(6,6); m.setIdentity();
            Q_=ekf_.getQmatrixIdentity();

            Q_.block(state::pos,state::pos,3,3)=Matrix3::Identity()*1.e-8;
            Q_.block(state::ori,state::ori,3,3)=Matrix3::Identity()*1.e-8;
            Q_.block(state::linVel,state::linVel,3,3)=Matrix3::Identity()*1.e-8;
            Q_.block(state::angVel,state::angVel,3,3)=Matrix3::Identity()*1.e-8;

            Q_.block(state::fc,state::fc,3,3)=Matrix3::Identity()*1.e-4;
            Q_.block(state::fc+3,state::fc+3,3,3)=Matrix3::Identity()*1.e-4;
            Q_.block(state::fc+6,state::fc+6,3,3)=Matrix3::Identity()*1.e-4;
            Q_.block(state::fc+6+3,state::fc+6+3,3,3)=Matrix3::Identity()*1.e-4;

            if(withUnmodeledForces_)
                Q_.block(state::unmodeledForces,state::unmodeledForces,6,6)=m*1.e-2;
            else
                Q_.block(state::unmodeledForces,state::unmodeledForces,6,6).setZero();

            if(withComBias_)
              Q_.block(state::comBias,state::comBias,2,2)=Matrix::Identity(2,2)*2.5e-10;
            else
              Q_.block(state::comBias,state::comBias,2,2).setZero();

            if (withAbsolutePos_)
              Q_.block(state::drift,state::drift,3,3)=Matrix::Identity(3,3)*1e-8;
            else
              Q_.block(state::drift,state::drift,3,3).setZero();

            ekf_.setQ(Q_);
            resetStateCovarianceMatrix();
    }

    void ModelBaseEKFFlexEstimatorIMU::resetStateCovarianceMatrix()
    {
            P_=Q_;
            ekf_.setStateCovariance(P_);
    }

    void ModelBaseEKFFlexEstimatorIMU::setContactsNumber(unsigned i)
    {
        functor_.setContactsNumber(i);

        inputSize_ = functor_.getInputSize();
        ekf_.setInputSize(inputSize_);

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
                    must be 23x1 vector, 6x1 vector or 4x4 matrix");

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

            x_s.segment(state::pos,3)=x0.head<3> ();

            x_s.segment(state::ori,3)=x0.tail<3> ();

            ekf_.setState(x_s,k_);

            ekf_.setQ(Q_);
        }
    }

    void ModelBaseEKFFlexEstimatorIMU::setComBiasGuess(const stateObservation::Vector & x)
    {
        lastX_.segment(state::comBias,2)=x.segment(0,2);
        setFlexibilityGuess(lastX_);
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

    Vector ModelBaseEKFFlexEstimatorIMU::getMomenta()
    {
        if(on_==true)
        {
            return functor_.getMomenta(getFlexibilityVector(),getInput());
        }
        else
        {
            Vector6 v; v.setZero();
            return v;
        }
    }

    Vector ModelBaseEKFFlexEstimatorIMU::getForcesAndMoments()
    {
        const Vector & v (getFlexibilityVector());

        Vector v2; v2.resize(functor_.getContactsNumber()*6);
        v2 << v.segment(state::fc,functor_.getContactsNumber()*6);

        return v2;
    }

    void ModelBaseEKFFlexEstimatorIMU::updateCovarianceMatrix_()
    {
        int currIndex = 6;
        R_.conservativeResize(getMeasurementSize(),getMeasurementSize());

        if(useFTSensors_)
        {
          R_.block(currIndex,0,functor_.getContactsNumber()*6,currIndex).setZero();
          R_.block(0,currIndex,currIndex,functor_.getContactsNumber()*6).setZero();
          R_.block(currIndex,currIndex,functor_.getContactsNumber()*6,functor_.getContactsNumber()*6) =
            Matrix::Identity(functor_.getContactsNumber()*6,functor_.getContactsNumber()*6)*forceVariance_;

          currIndex += functor_.getContactsNumber()*6;
        }

        if(withAbsolutePos_)
        {
          R_.block(currIndex,0,6,currIndex).setZero();
          R_.block(0,currIndex,currIndex,6).setZero();
          R_.block(currIndex,currIndex,6,6) =   Matrix::Identity(6,6)*absPosVariance_;

          currIndex += 6;
        }

        ekf_.setR(R_);
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getStateSize() const
    {
        return ekf_.getStateSize();
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getInputSize() const
    {
        return ekf_.getInputSize();
    }

    unsigned ModelBaseEKFFlexEstimatorIMU::getMeasurementSize() const
    {
        return functor_.getMeasurementSize();
    }

    Matrix4 ModelBaseEKFFlexEstimatorIMU::getFlexibility()
    {

        const Vector & v (getFlexibilityVector());

        Vector6 v2;v2 << v.segment(state::pos,3), v.segment(state::ori,3);
        //v2.head<3>() = v.segment(kine::pos,3);
        //v2.tail<3>() = v.segment(kine::ori,3);

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

      if(on_==true)
      {
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
        preIterationCallback_();
        if (ekf_.getMeasurementsNumber()>0)
        {
          k_=ekf_.getMeasurementTime();

          unsigned i;
          for (i=ekf_.getCurrentTime()+1; i<=k_; ++i)
          {
            if (finiteDifferencesJacobians_)
            {
              ekf_.updatePredictedMeasurement();///triggers also ekf_.updatePrediction();

              //ekf_.setA(ekf_.getAMatrixFD(dx_));
              //ekf_.setC(ekf_.getCMatrixFD(dx_));
              ekf_.setA(functor_.stateDynamicsJacobian());
              ekf_.setC(functor_.measureDynamicsJacobian());
            }
            ekf_.getEstimatedState(i);
          }
          x_=ekf_.getEstimatedState(k_);

          if (! x_.hasNaN())//detect NaN values
          {
            lastX_=x_;

            ///regulate the part of orientation vector in the state vector
            lastX_.segment(state::ori,3)=kine::regulateOrientationVector(lastX_.segment(state::ori,3));
            if(limitOn_)
            {
              for(int i=0; i<3; i++)
              {
                // Saturation for bounded forces and torques
                lastX_[state::fc+6+i]=std::min(lastX_[state::fc+6+i],limitTorques_[i]);
                lastX_[state::fc+i]=std::min(lastX_[state::fc+i],limitForces_[i]);
                lastX_[state::fc+6+i]=std::max(lastX_[state::fc+6+i],-limitTorques_[i]);
                lastX_[state::fc+i]=std::max(lastX_[state::fc+i],-limitForces_[i]);
              }
            }
            ekf_.setState(lastX_,ekf_.getCurrentTime());
          }
          else //delete NaN values
          {
            ekf_.setState(lastX_,k_);

            if(k_>1) //the first iteration give always nan when not
              //initialized
            {
              resetCovarianceMatrices();
              resetStateCovarianceMatrix();
            }
          }
        }
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time3);

        computeFlexibilityTime_=(double)diff(time2,time3).tv_nsec-(double)diff(time1,time2).tv_nsec;

      }
      else
      {
        lastX_.setZero();
        ekf_.setState(lastX_,ekf_.getCurrentTime());

        if (ekf_.getMeasurementsNumber()>0)
        {
          ekf_.clearMeasurements();
          ekf_.clearInputs();
          resetStateCovarianceMatrix();
        }
      }

      return lastX_;
    }

    void ModelBaseEKFFlexEstimatorIMU::preIterationCallback_()
    {
      if (callback_.updateUnmodForceStateCov)
      {
        P_=ekf_.getStateCovariance();
        P_.diagonal().segment<6>(state::unmodeledForces).setConstant(unmodeledForceVariance_);
        ekf_.setStateCovariance(P_);
        callback_.updateUnmodForceStateCov=false;
      }
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

    Matrix ModelBaseEKFFlexEstimatorIMU::getKfe() const
    {
        return functor_.getKfe();
    }

    Matrix ModelBaseEKFFlexEstimatorIMU::getKfv() const
    {
        return functor_.getKfv();
    }

    Matrix ModelBaseEKFFlexEstimatorIMU::getKte() const
    {
        return functor_.getKte();
    }

    Matrix ModelBaseEKFFlexEstimatorIMU::getKtv() const
    {
        return functor_.getKtv();
    }

    double& ModelBaseEKFFlexEstimatorIMU::getComputeFlexibilityTime()
    {
        return computeFlexibilityTime_;
    }

    void ModelBaseEKFFlexEstimatorIMU::setWithForcesMeasurements(bool b)
    {
      if (useFTSensors_!= b)
      {
        useFTSensors_=b;
        functor_.setWithForceMeasurements(b);
        ekf_.setMeasureSize(functor_.getMeasurementSize());

        updateCovarianceMatrix_();
      }
    }

    void ModelBaseEKFFlexEstimatorIMU::setWithAbsolutePos(bool b)
    {
      if (withAbsolutePos_!= b)
      {
        functor_.setWithAbsolutePosition(b);
        ekf_.setMeasureSize(functor_.getMeasurementSize());
        withAbsolutePos_=b;
        updateCovarianceMatrix_();
      }
    }

    void ModelBaseEKFFlexEstimatorIMU::setWithUnmodeledForces(bool b)
    {
      if (withUnmodeledForces_!= b)
      {
        functor_.setWithUnmodeledForces(b);
        ekf_.setMeasureSize(functor_.getMeasurementSize());
        ekf_.setInputSize(functor_.getInputSize());
        withUnmodeledForces_=b;
        updateCovarianceMatrix_();
      }
    }

    bool ModelBaseEKFFlexEstimatorIMU::getWithForcesMeasurements()
    {
        return useFTSensors_;
    }


    void ModelBaseEKFFlexEstimatorIMU::setWithComBias(bool b)
    {
      if (withComBias_!= b)
      {
        withComBias_=b;
        functor_.setWithComBias(b);
      }
    }

    void ModelBaseEKFFlexEstimatorIMU::setUnmodeledForceVariance(double d)
    {
        unmodeledForceVariance_ = d;
        callback_.updateUnmodForceStateCov=true;
        updateCovarianceMatrix_();
    }

    void ModelBaseEKFFlexEstimatorIMU::setForceVariance(double d)
    {
        forceVariance_ = d;
        updateCovarianceMatrix_();
    }

    void ModelBaseEKFFlexEstimatorIMU::setAbsolutePosVariance(double d)
    {
        absPosVariance_ = d;
        updateCovarianceMatrix_();
    }

    void ModelBaseEKFFlexEstimatorIMU::setRobotMass(double m)
    {
        functor_.setRobotMass(m);
    }

}
}
