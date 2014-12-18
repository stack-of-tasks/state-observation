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
            if ((p_>0) && (directInputStateProcessFeedthrough_))
            {

                BOOST_ASSERT(this->u_.size()>0 && this->u_.checkIndex(k-1) &&
                                        "ERROR: The input vector is not set");
                opt.u_=this->u_[k-1];
            }
            else
            {
                opt.u_ = inputVectorZero();
            }

            BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
            xbar_.set(f_->stateDynamics(
                      this->x_(),
                      opt.u_,
                      this->x_.getTime()),
                      k);
        }


        return xbar_();
    }

    ObserverBase::StateVector ExtendedKalmanFilter::getPrediction()
    {
        return prediction_(x_.getTime()+1);
    }

    ObserverBase::MeasureVector ExtendedKalmanFilter::simulateSensor_(const ObserverBase::StateVector& x, unsigned k)
    {
        BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");

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
                opt.u_=u_[k];
            }
            else
            {
                opt.u_=inputVectorZero();
            }
        }

        return f_->measureDynamics(x,opt.u_,k);
    }

    KalmanFilterBase::Amatrix// ExtendedKalmanFilter<n,m,p>::Amatrix does not work
    ExtendedKalmanFilter::getAMatrixFD(const ObserverBase::StateVector
                                       &dx)
    {
        unsigned k=this->x_.getTime();
        opt.a_=getAmatrixZero();
        opt.xbar_=prediction_(k+1);
        opt.x_=this->x_();

        if (p_>0)
        {
            if (directInputStateProcessFeedthrough_)
                opt.u_=this->u_[k];
            else
                opt.u_=inputVectorZero();
        }

        for (unsigned i=0;i<n_;++i)
        {
            unsigned it=(i-1)%n_;
            opt.x_[it]=this->x_()(it,0);
            opt.x_[i]+=dx[i];
            opt.xp_=f_->stateDynamics(opt.x_,opt.u_,k);

            opt.xp_-=opt.xbar_;
            opt.xp_/=dx[i];

            opt.a_.block(0,i,n_,1)=opt.xp_;
        }

        return opt.a_;
    }

    KalmanFilterBase::Cmatrix
    ExtendedKalmanFilter::getCMatrixFD(const ObserverBase::StateVector
                                       &dx)
    {
        unsigned k=this->x_.getTime();

        opt.c_.resize(m_,n_);

        opt.xbar_=prediction_(k+1);
        opt.xp_ = opt.xbar_;

        opt.y_=simulateSensor_( opt.xbar_, k+1);

        for (unsigned i=0;i<n_;++i)
        {
            unsigned it=(i-1)%n_;
            opt.xp_[it]=opt.xbar_[it];
            opt.xp_[i]+= dx[i];

            opt.yp_=simulateSensor_(opt.xp_, k+1);
            opt.yp_-=opt.y_;
            opt.yp_/=dx[i];

            opt.c_.block(0,i,m_,1)=opt.yp_;

        }

        return opt.c_;
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
