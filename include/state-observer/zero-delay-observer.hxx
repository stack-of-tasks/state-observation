template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::setState
(const typename ObserverBase<n,m,p>::StateVector& x_k,unsigned k)
{
    x_.set(x_k,k);
    while (y_.size()>0 && y_[0].getTime()<=k)
    {
        y_.pop_front();
    }

    if (p>0)
        while (u_.size()>0 && u_[0].getTime()<k)
        {
            u_.pop_front();
        }
}

template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::clearState()
{
    x_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::setMeasurement
(const typename ObserverBase<n,m,p>::MeasureVector& y_k,unsigned k)
{
    if (y_.size()>0)
        BOOST_ASSERT (y_[y_.size()-1].getTime()==k-1 && "ERROR: The time is set incorrectly for the measurements (order or gap)");
    else
        BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k-1) && "ERROR: The time is set incorrectly for the measurements (must be [current_time+1])");

    typename ObserverBase<n,m,p>::Measure a(y_k,k);

    y_.push_back(a);
}

template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::clearMeasurements()
{
    y_.clear();
}

template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::setInput
(const typename ObserverBase<n,m,p>::InputVector& u_k,unsigned k)
{
    if (p>0)
    {
        if (u_.size()>0)
            BOOST_ASSERT (u_[u_.size()-1].getTime()==k-1 && "ERROR: The time is set incorrectly for the inputs (order or gap)");
        else
            BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k) && "ERROR: The time is set incorrectly for the inputs (must be [current_time])");

        typename ObserverBase<n,m,p>::Input a(u_k,k);

        u_.push_back(a);
    }
}

template <unsigned n,unsigned m, unsigned p>
void ZeroDelayObserver<n,m,p>::clearInputs()
{
    if (p>0)
        u_.clear();
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector
ZeroDelayObserver<n,m,p>::getEstimateState(unsigned k)
{
    unsigned k0=x_.getTime();

    BOOST_ASSERT(k0<=k && "ERROR: The observer cannot estimate previous states");

    for (unsigned i=k0;i<k;++i)
    {
        oneStepEstimation_();
        y_.pop_front();
        if (p>0)
            u_.pop_front();
    }

    return x_();
}


template <unsigned n,unsigned m, unsigned p>
unsigned ZeroDelayObserver<n,m,p>::getCurrentTime()const
{
    return x_.getTime();
}
