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
    if (p>0)
    {
        directInputOutputFeedthrough_=b;
    }
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector ExtendedKalmanFilter<n,m,p>::prediction_(unsigned k)
{
    typename ObserverBase<n,m,p>::InputVector u;

    if (p>0)
        u=this->u_[0]();


    if (!this->xbar_.isSet() || this->xbar_.getTime()!=k)
    {
        BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
        xbar_.set(f_->stateDynamics(
                      this->x_(),
                      u,
                      this->x_.getTime()),
                  k);
    }
    return xbar_();
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector ExtendedKalmanFilter<n,m,p>::getPrediction(unsigned k)
{
    BOOST_ASSERT(k==this->x_.getTime()+1 && "ERROR: The prediction can only be calculated for next sample (k+1)");
    if (p>0)
    {
        BOOST_ASSERT(this->u_.size()>0 && this->u_[0].getTime()== k-1 && "ERROR: The input vector is not set");
    }
    return prediction_(k);
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::MeasureVector
ExtendedKalmanFilter<n,m,p>::simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k)
{
    BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
    typename ObserverBase<n,m,p>::InputVector u (ObserverBase<n,m,p>::InputVector::Zero());
    unsigned i;
    if (p>0)
    {
        for (i=0; i<this->u_.size()&&this->u_[i].getTime()<k;++i)
        {
        }
        if (directInputOutputFeedthrough_)
        {
            BOOST_ASSERT(i!=this->u_.size() && this->u_[i].getTime()==k &&
                         "ERROR: The input feedthrough of the measurements is not set \
                         (the measurement at time k needs the input at time k which was not given) \
                         if you don't need the input in the computation of measurement, you \
                         must set directInputOutputFeedthrough to 'false' in the constructor");
        }

        if (i<this->u_.size())
            u=this->u_[i]();
    }

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

    typename ObserverBase<n,m,p>::InputVector u;

    if (p>0)
        u=this->u_[0]();

      //  cout << "input2: " << u << endl;

    for (unsigned i=0;i<n;++i)
    {
        unsigned it=(i-1)%n;
        x[it]=this->x_()[it];
        x[i]=this->x_()[i]+dx[i];
        xp=(f_->stateDynamics(x,u,k)-fx)/dx[i];

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


