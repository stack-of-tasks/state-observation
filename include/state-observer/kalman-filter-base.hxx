template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::setA(const typename KalmanFilterBase<n,m,p>::Amatrix& A)
{
    a_.set(A,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::clearA()
{
    a_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::setC( const typename KalmanFilterBase<n,m,p>::Cmatrix& C)
{
    c_.set(C,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::clearC()
{
    c_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::setR( const typename  KalmanFilterBase<n,m,p>::Rmatrix& R)
{
    r_.set(R,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::clearR()
{
    r_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::setQ( const typename   KalmanFilterBase<n,m,p>::Qmatrix& Q)
{
    q_.set(Q,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::clearQ()
{
    q_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::setStateCovariance(const typename   KalmanFilterBase<n,m,p>::Pmatrix& P)
{
    p_.set(P,this->x_.getTime());
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::clearStateCovariance()
{
    p_.reset();
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector KalmanFilterBase<n,m,p>::oneStepEstimation_()
{
    unsigned k=this->x_.getTime();
    BOOST_ASSERT(this->y_.size()> 0 && this->y_[0].getTime()==k+1 && "ERROR: The measurement vector is not set");
    BOOST_ASSERT(this->u_.size()> 0 && this->u_[0].getTime()==k && "ERROR: The input vector is not set");
    BOOST_ASSERT(a_.isSet() && "ERROR: The Matrix A is not initialized" );
    BOOST_ASSERT(c_.isSet() && "ERROR: The Matrix C is not initialized");
    BOOST_ASSERT(q_.isSet() && "ERROR: The Matrix Q is not initialized");
    BOOST_ASSERT(r_.isSet() && "ERROR: The Matrix R is not initialized");
    BOOST_ASSERT(p_.isSet() && "ERROR: The Matrix P is not initialized");

    Amatrix a=a_();
    Cmatrix c=c_();
    Pmatrix px=p_();

    //prediction
    typename ObserverBase<n,m,p>::StateVector xbar=prediction_(k+1);
    Pmatrix pbar=a*px*a.transpose()+q_();

    typename ObserverBase<n,m,p>::InputVector u;

    //innovation
    typename ObserverBase<n,m,p>::MeasureVector ino= this->y_[0]() - simulateSensor_(xbar,k+1);
    Rmatrix inoCov = c * pbar * c.transpose() + r_();

    //gain
    Kmatrix kGain = (pbar * c.transpose()) * inoCov.inverse();

    //update
    typename ObserverBase<n,m,p>::StateVector xhat=xbar+kGain*ino;

    this->x_.set(xhat,k+1);
    p_.set((Pmatrix::Identity()-kGain*c)*pbar,k+1);

    return xhat;
}

template <unsigned n,unsigned m, unsigned p>
typename KalmanFilterBase<n,m,p>::Pmatrix KalmanFilterBase<n,m,p>::getStateCovariance(unsigned k)
{
    this->getEstimateState(k);
    return p_();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilterBase<n,m,p>::reset()
{
    ZeroDelayObserver<n,m,p>::reset();

    a_.reset();
    c_.reset();
    a_.reset();
    c_.reset();
    q_.reset();
    r_.reset();
    p_.reset();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::setA( const typename KalmanFilterBase<n,m,0>::Amatrix& A)
{
    a_.set(A,0);
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::clearA()
{
    a_.reset();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::setC( const typename  KalmanFilterBase<n,m,0>::Cmatrix& C)
{
    c_.set(C,0);
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::clearC()
{
    c_.reset();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::setR(const typename  KalmanFilterBase<n,m,0>::Rmatrix& R)
{
    r_.set(R,0);
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::clearR()
{
    r_.reset();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::setQ(const typename   KalmanFilterBase<n,m,0>::Qmatrix& Q)
{
    q_.set(Q,0);
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::clearQ()
{
    q_.reset();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::setStateCovariance(const Pmatrix& P)
{
    p_.set(P,this->x_.getTime());
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::clearStateCovariance()
{
    p_.reset();
}


template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::StateVector KalmanFilterBase<n,m,0>::oneStepEstimation_()
{
    unsigned k=this->x_.getTime();
    BOOST_ASSERT(this->y_.size()!=0 && this->y_[0].getTime()==k+1 && "ERROR: The measurement vector is not set");
    BOOST_ASSERT(a_.isSet() && "ERROR: The Matrix A is not initialized" );
    BOOST_ASSERT(c_.isSet() && "ERROR: The Matrix C is not initialized");
    BOOST_ASSERT(q_.isSet() && "ERROR: The Matrix Q is not initialized");
    BOOST_ASSERT(r_.isSet() && "ERROR: The Matrix R is not initialized");
    BOOST_ASSERT(p_.isSet() && "ERROR: The Matrix P is not initialized");

    Amatrix a=a_();
    Cmatrix c=c_();
    Pmatrix px=p_();

    //prediction
    typename ObserverBase<n,m,0>::StateVector xbar=prediction_(k+1);
    Pmatrix pbar=a*px*a.transpose()+q_();

    //innovation
    typename ObserverBase<n,m,0>::MeasureVector ino= this->y_[0]() - simulateSensor_(xbar,k+1);
    Rmatrix inoCov = c * pbar * c.transpose() + r_();

    //gain
    Kmatrix kGain = pbar * c.transpose() * inoCov.inverse();

    //update
    typename ObserverBase<n,m,0>::StateVector xhat=xbar+kGain*ino;

    this->x_.set(xhat,k+1);
    p_.set((Pmatrix::Identity()-kGain*c)*pbar,k+1);

    return xhat;
}


template <unsigned n,unsigned m>
typename KalmanFilterBase<n,m,0>::Pmatrix KalmanFilterBase<n,m,0>::getStateCovariance(unsigned k)
{
    this->getEstimateState(k);
    return p_();
}

template <unsigned n,unsigned m>
void KalmanFilterBase<n,m,0>::reset()
{
    ZeroDelayObserver<n,m,0>::reset();

    a_.reset();
    c_.reset();
    a_.reset();
    c_.reset();
    q_.reset();
    r_.reset();
    p_.reset();
}
