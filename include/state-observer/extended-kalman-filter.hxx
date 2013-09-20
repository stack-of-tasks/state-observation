template <unsigned n,unsigned m, unsigned p>
void ExtendedKalmanFilter<n,m,p>::setFunctor(typename ExtendedKalmanFilter<n,m,p>::DynamicsFunctorBase* f)
{
    f_=f;
    //f_->reset();

}

template <unsigned n,unsigned m, unsigned p>
void ExtendedKalmanFilter<n,m,p>::clearFunctor()
{
    f_=0x0;
}

template <unsigned n,unsigned m, unsigned p>
void ExtendedKalmanFilter<n,m,p>::setDirectInputOutputFeedthrough(bool b)
{
    directInputOutputFeedthrough_=b;
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector ExtendedKalmanFilter<n,m,p>::prediction_(unsigned k)
{
    if (!this->xbar_.isSet() || this->xbar_.getTime()!=k)
    {
        BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
        xbar_.set(f_->stateDynamics(
                      this->x_(),
                      this->u_[0](),
                      this->x_.getTime()),
                  k);
    }
    return xbar_();
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector ExtendedKalmanFilter<n,m,p>::getPrediction(unsigned k)
{
    BOOST_ASSERT(k==this->x_.getTime()+1 && "ERROR: The prediction can only be calculated for next sample (k+1)");
    BOOST_ASSERT(this->u_.size()>0 && this->u_[0].getTime()== k-1 && "ERROR: The input vector is not set");
    return prediction_(k);
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::MeasureVector
ExtendedKalmanFilter<n,m,p>::simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k)
{
    BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
    typename ObserverBase<n,m,p>::InputVector u= ObserverBase<n,m,p>::InputVector::Zero();
    unsigned i;
    for (i=0; i<this->u_.size()&&this->u_[i].getTime()<k;++i)
    {
    }
    if (directInputOutputFeedthrough_)
    {
        BOOST_ASSERT(i!=this->u_.size() && this->u_[i].getTime()==k &&
            "The input feedthrough of the measurements is not set\
                (if you don't need the input in the computation measurement, you need\
                to set directInputOutputFeedthrough to false in the constructor)");
    }
    if (i<this->u_.size())
        u=this->u_[i]();

    return f_->measureDynamics(x,u,k);
}

template <unsigned n,unsigned m, unsigned p>
typename KalmanFilterBase<n,m,p>::Amatrix// ExtendedKalmanFilter<n,m,p>::Amatrix does not work
ExtendedKalmanFilter<n,m,p>::getAMatrixFD(const typename ObserverBase<n,m,p>::StateVector
        &dx)
{
    unsigned k=this->x_.getTime();
    typename KalmanFilterBase<n,m,p>::Amatrix a;
    typename ObserverBase<n,m,p>::StateVector fx=prediction_(k+1);
    typename ObserverBase<n,m,p>::StateVector x=this->x_();
    typename ObserverBase<n,m,p>::StateVector xp;
    for (unsigned i=0;i<n;++i)
    {
        unsigned it=(i-1)%n;
        x[it]=this->x_()[it];
        x[i]=this->x_()[i]+dx[i];
        xp=(f_->stateDynamics(x,this->u_[0](),k)-fx)/dx[i];

        for (unsigned j=0;j<n;++j)
        {
            a(j,i)=xp[j];
        }
    }

    return a;
}

template <unsigned n,unsigned m, unsigned p>
typename KalmanFilterBase<n,m,p>::Cmatrix//typename ExtendedKalmanFilter<n,m,p>::Cmatrix does not work
ExtendedKalmanFilter<n,m,p>::getCMatrixFD(const typename ObserverBase<n,m,p>::StateVector
        &dx)
{
    unsigned k=this->x_.getTime();
    typename KalmanFilterBase<n,m,p>::Cmatrix c;
    typename ObserverBase<n,m,p>::MeasureVector y=simulateSensor_(this->x_(), k);
    typename ObserverBase<n,m,p>::StateVector x=this->x_();
    typename ObserverBase<n,m,p>::MeasureVector yp;

    for (unsigned i=0;i<n;++i)
    {
        x[(i-1)%n]=this->x_()[(i-1)%n];
        x[i]=this->x_()[i]+dx[i];
        yp=(simulateSensor_(x, k)-y)/dx[i];

        for (unsigned j=0;j<m;++j)
        {
            c(j,i)=yp[j];
        }
    }
    return c;
}

template <unsigned n,unsigned m, unsigned p>
void ExtendedKalmanFilter<n,m,p>::reset()
{
    KalmanFilterBase<n,m,p>::reset();
    if (f_!=0x0)
        f_->reset();
}



template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::StateVector ExtendedKalmanFilter<n,m,0>::prediction_(unsigned k)
{
    if (!(this->xbar_.isSet()) || (this->xbar_.getTime()!=k))
    {
        BOOST_ASSERT(f_!=0x0 && "ERROR: The Kalman filter functor is not set");
        this->xbar_.set(f_->stateDynamics(this->x_(),this->x_.getTime()),k);
    }
    return xbar_();
}

template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::StateVector ExtendedKalmanFilter<n,m,0>::getPrediction(unsigned k)
{
    BOOST_ASSERT(k==this->x_.getTime()+1 && "ERROR: The prediction can only be calculated for next sample (k+1)");
    return prediction_(k);
}

template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::MeasureVector
ExtendedKalmanFilter<n,m,0>::simulateSensor_(const typename ObserverBase<n,m,0>::StateVector& x, unsigned k)
{
    BOOST_ASSERT(f_!=0x0 && "ERROR: The Kalman filter functor is not set");
    return f_->measureDynamics(x,k);
}

template <unsigned n,unsigned m>
void ExtendedKalmanFilter<n,m,0>::setFunctor(typename ExtendedKalmanFilter<n,m,0>::DynamicsFunctorBase* p)
{
    f_=p;
}

template <unsigned n,unsigned m>
void ExtendedKalmanFilter<n,m,0>::clearFunctor()
{
    f_=0x0;
}



template <unsigned n,unsigned m>
typename KalmanFilterBase<n,m,0>::Amatrix  //typename ExtendedKalmanFilter<n,m,0>::Amatrix does not work
ExtendedKalmanFilter<n,m,0>::getAMatrixFD(const typename ObserverBase<n,m,0>::StateVector &dx)
{
    unsigned k=this->x_.getTime();
    typename KalmanFilterBase<n,m,0>::Amatrix a;
    typename ObserverBase<n,m,0>::stateVector fx=prediction_(k+1);
    typename ObserverBase<n,m,0>::StateVector x=this->x_();
    typename ObserverBase<n,m,0>::StateVector xp;
    for (unsigned i=0;i<n;++i)
    {
        x[i-1%n]=this->x_()[i-1%n];
        x[i]=this->x_[i]+dx[i];
        xp=(f_->stateDynamics(x,k)-fx)/dx[i];

        for (unsigned j=0;j<n;++j)
        {
            a(j,i)=xp[j];
        }
    }
    return a;
}

template <unsigned n,unsigned m>
typename KalmanFilterBase<n,m,0>::Cmatrix //typename ExtendedKalmanFilter<n,m,0>::Cmatrix does not work
ExtendedKalmanFilter<n,m,0>::getCMatrixFD(const typename ObserverBase<n,m,0>::StateVector &dx)
{
    unsigned k=this->x_.getTime();
    typename KalmanFilterBase<n,m,0>::Cmatrix c;
    typename ObserverBase<n,m,0>::MeasureVector y=simulateSensor_(this->x_(), k);
    typename ObserverBase<n,m,0>::StateVector x=this->x_;
    typename ObserverBase<n,m,0>::MeasureVector yp;

    for (unsigned i=0;i<n;++i)
    {
        x[i-1%n]=this->x_()[i-1%n];
        x[i]=this->x_[i]+dx[i];
        yp=(simulateSensor_(x, k)-y)/dx[i];

        for (unsigned j=0;j<m;++j)
        {
            c(j,i)=yp[j];
        }
    }
    return c;
}

template <unsigned n,unsigned m>
void ExtendedKalmanFilter<n,m,0>::reset()
{
    KalmanFilterBase<n,m,0>::reset();
    if (f_!=0x0)
        f_->reset();
}
