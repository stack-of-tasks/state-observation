void ExtendedKalmanFilter::setFunctor(ExtendedKalmanFilter::DynamicsFunctorBase* f)
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

ObserverBase::StateVector ExtendedKalmanFilter::prediction_(unsigned k)
{
    ObserverBase::InputVector u;

    if (p_>0)
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

ObserverBase::StateVector ExtendedKalmanFilter::getPrediction(unsigned k)
{
    BOOST_ASSERT(k==this->x_.getTime()+1 && "ERROR: The prediction can only be calculated for next sample (k+1)");
    if (p_>0)
    {
        BOOST_ASSERT(this->u_.size()>0 && this->u_[0].getTime()== k-1 && "ERROR: The input vector is not set");
    }
    return prediction_(k);
}

ObserverBase::MeasureVector ExtendedKalmanFilter::simulateSensor_(const ObserverBase::StateVector& x, unsigned k)
{
    BOOST_ASSERT (f_!=0x0 && "ERROR: The Kalman filter functor is not set");
    ObserverBase::InputVector u (inputVectorZero());
    unsigned i;
    if (p_>0)
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

KalmanFilterBase::Amatrix// ExtendedKalmanFilter<n,m,p>::Amatrix does not work
ExtendedKalmanFilter::getAMatrixFD(const ObserverBase::StateVector
        &dx)
{
    unsigned k=this->x_.getTime();
    Amatrix a(getAmatrixZero());
    StateVector fx=prediction_(k+1);
    StateVector x=this->x_();
    StateVector xp;

    InputVector u;

    if (p_>0)
        u=this->u_[0]();

    for (unsigned i=0;i<n_;++i)
    {
        unsigned it=(i-1)%n_;
        x(it,0)=this->x_()(it,0);
        x(i,0)=this->x_()(i,0)+dx(i,0);
        xp=(f_->stateDynamics(x,u,k)-fx)/dx[i];

        for (unsigned j=0;j<n_;++j)
        {
            a(j,i)=xp[j];
        }
    }

    return a;
}

KalmanFilterBase::Cmatrix//typename ExtendedKalmanFilter<n,m,p>::Cmatrix does not work
ExtendedKalmanFilter::getCMatrixFD(const ObserverBase::StateVector
        &dx)
{
    unsigned k=this->x_.getTime();
    Cmatrix c(getCmatrixZero());
    MeasureVector y=simulateSensor_(this->x_(), k);
    StateVector x=this->x_();
    MeasureVector yp;

    for (unsigned i=0;i<n_;++i)
    {
        x[(i-1)%n_]=this->x_()((i-1)%n_,0);
        x[i]=this->x_()(i,0)+dx[i];
        yp=(simulateSensor_(x, k)-y)/dx[i];

        for (unsigned j=0;j<m_;++j)
        {
            c(j,i)=yp[j];
        }
    }
    return c;
}

void ExtendedKalmanFilter::reset()
{
    KalmanFilterBase::reset();
    if (f_!=0x0)
        f_->reset();
}


