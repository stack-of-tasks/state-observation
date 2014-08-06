#include <state-observation/observer/extended-kalman-filter.hpp>

namespace stateObservation
{
    void ExtendedKalmanFilter::setFunctor(DynamicalSystemFunctorBase* f)
    {
        f_=f;
        //f_->reset();

    }

    void ExtendedKalmanFilter::clearFunctor()
    {
        f_=0x0;
    }

    void ExtendedKalmanFilter::setDirectInputOutputFeedthrough(bool b)
    {
        if (p_>0)
        {
            directInputOutputFeedthrough_=b;
        }
    }

    void ExtendedKalmanFilter::setDirectInputStateFeedthrough(bool b)
    {
        if (p_>0)
        {
            directInputStateProcessFeedthrough_=b;
        }
    }

    ObserverBase::StateVector ExtendedKalmanFilter::prediction_(unsigned k)
    {


        if (!this->xbar_.isSet() || this->xbar_.getTime()!=k)
        {
            ObserverBase::InputVector u;


            if ((p_>0) && (directInputStateProcessFeedthrough_))
            {

                BOOST_ASSERT(this->u_.size()>0 && this->u_.checkIndex(k-1) &&
                                        "ERROR: The input vector is not set");
                u=this->u_[k-1];
            }
            else
            {
                u = inputVectorZero();
            }

            //std::cout << "u" << u << std::endl;

            BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
           // std::cout << "calcul stateDynmic" << std::endl;
            xbar_.set(f_->stateDynamics(
                          this->x_(),
                          u,
                          this->x_.getTime()),
                      k);
        }

        //std::cout << "Je suis passé dans prediction_ et xbar_=" << xbar_().transpose() << std::endl;

        return xbar_();
    }

    ObserverBase::StateVector ExtendedKalmanFilter::getPrediction()
    {
      //  std::cout << "getPrediction -> prediction_" << std::endl;
        return prediction_(x_.getTime()+1);
    }

    ObserverBase::MeasureVector ExtendedKalmanFilter::simulateSensor_(const ObserverBase::StateVector& x, unsigned k)
    {
        BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
        ObserverBase::InputVector u;
        unsigned i;
        if (p_>0)
        {
            if (directInputOutputFeedthrough_)
            {
                BOOST_ASSERT(u_.checkIndex(k) &&
                "ERROR: The input feedthrough of the measurements is not set \
(the measurement at time k needs the input at time k which was not given) \
if you don't need the input in the computation of measurement, you \
must set directInputOutputFeedthrough to 'false' in the constructor");
            }

            if (u_.checkIndex(k))
            {
                u=u_[k];
            }
            else
            {
                u=inputVectorZero();
            }
        }

        return f_->measureDynamics(x,u,k);
    }

    KalmanFilterBase::Amatrix// ExtendedKalmanFilter<n,m,p>::Amatrix does not work
    ExtendedKalmanFilter::getAMatrixFD(const ObserverBase::StateVector
                                       &dx)
    {
        unsigned k=this->x_.getTime();
        Amatrix a(getAmatrixZero());
       // std::cout << "getAMatrix -> prediction_" << std::endl;
        StateVector fx=prediction_(k+1);
        StateVector x=this->x_();
        StateVector xp;

        InputVector u;

        if (p_>0)
            if (directInputStateProcessFeedthrough_)
                u=this->u_[k];
            else
                u=inputVectorZero();


        for (unsigned i=0;i<n_;++i)
        {
            unsigned it=(i-1)%n_;
            x[it]=this->x_()(it,0);
            x[i]=this->x_()(i,0);
            x[i]+=dx[i];
          //  std::cout << "getAMatrixFD -> stateDynamics " << std::endl;
            xp=f_->stateDynamics(x,u,k);

            xp-=fx;
            xp/=dx[i];

            a.block(0,i,n_,1)=xp;
        }

        return a;
    }

    KalmanFilterBase::Cmatrix
    ExtendedKalmanFilter::getCMatrixFD(const ObserverBase::StateVector
                                       &dx)
    {
        unsigned k=this->x_.getTime();

        Cmatrix c(getCmatrixZero());

      //  std::cout << "getCMatrix -> prediction_" << std::endl;
        StateVector xbar=prediction_(k+1);
        StateVector xbarInit = xbar;

        MeasureVector y=simulateSensor_( xbar, k+1);

        MeasureVector yp;

        for (unsigned i=0;i<n_;++i)
        {
            xbar[(i-1)%n_]=xbarInit[(i-1)%n_];
            xbar[i]=xbarInit[i];
            xbar[i]+= dx[i];

            yp=simulateSensor_(xbar, k+1);
            yp-=y;
            yp/=dx[i];

            //std::cout << "yi3" << yp.transpose() <<std::endl;

            c.block(0,i,m_,1)=yp;

        }

        return c;
    }

    void ExtendedKalmanFilter::reset()
    {
        KalmanFilterBase::reset();
        if (f_!=0x0)
            f_->reset();
    }

    DynamicalSystemFunctorBase* ExtendedKalmanFilter::functor() const
    {
        return f_;
    }

}
