#include <state-observer/zero-delay-observer.hpp>

namespace stateObserver
{

    void ZeroDelayObserver::setState (const ObserverBase::StateVector& x_k,unsigned k)
    {
        BOOST_ASSERT(checkStateVector(x_k) && "The size of the state vector is incorrect");
        x_.set(x_k,k);
        while (y_.size()>0 && y_[0].getTime()<=k)
        {
            y_.pop_front();
        }

        if (p_>0)
            while (u_.size()>0 && u_[0].getTime()<k)
            {
                u_.pop_front();
            }
    }

    void ZeroDelayObserver::clearStates()
    {
        x_.reset();
    }

    void ZeroDelayObserver::setMeasurement (const ObserverBase::MeasureVector& y_k,unsigned k)
    {
        BOOST_ASSERT(checkMeasureVector(y_k) && "The size of the measure vector is incorrect");
        if (y_.size()>0)
            BOOST_ASSERT (y_[y_.size()-1].getTime()==k-1 && "ERROR: The time is set incorrectly for the measurements (order or gap)");
        else
            BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k-1) && "ERROR: The time is set incorrectly for the measurements (must be [current_time+1])");

        DiscreteTimeMatrix a(y_k,k);

        y_.push_back(a);
    }

    void ZeroDelayObserver::clearMeasurements()
    {
        y_.clear();
    }

    void ZeroDelayObserver::setInput (const ObserverBase::InputVector& u_k,unsigned k)
    {
        if (p_>0)
        {
            BOOST_ASSERT(checkInputVector(u_k) && "The size of the input vector is incorrect");

            if (u_.size()>0)
                BOOST_ASSERT (u_[u_.size()-1].getTime()==k-1 && "ERROR: The time is set incorrectly for the inputs (order or gap)");
            else
                BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k) && "ERROR: The time is set incorrectly for the inputs (must be [current_time])");

            DiscreteTimeMatrix a(u_k,k);

            u_.push_back(a);
        }
    }

    void ZeroDelayObserver::clearInputs()
    {
        if (p_>0)
            u_.clear();
    }

    ObserverBase::StateVector
    ZeroDelayObserver::getEstimateState(unsigned k)
    {
        unsigned k0=x_.getTime();

        BOOST_ASSERT(k0<=k && "ERROR: The observer cannot estimate previous states");

        for (unsigned i=k0;i<k;++i)
        {
            oneStepEstimation_();
            y_.pop_front();
            if (p_>0)
                u_.pop_front();
        }

        return x_();
    }


    unsigned ZeroDelayObserver::getCurrentTime()const
    {
        return x_.getTime();
    }

    void ZeroDelayObserver::setStateSize(unsigned n)
    {
        if (n!=n_)
        {
            ObserverBase::setStateSize(n);
            clearStates();
        }
    }

    void ZeroDelayObserver::setMeasureSize(unsigned m)
    {
        if (m!=m_)
        {
            ObserverBase::setMeasureSize(m);
            clearMeasurements();
        }
    }

    void ZeroDelayObserver::setInputSize(unsigned p)
    {
        if (p!=p_)
        {
            ObserverBase::setInputSize(p);
            clearInputs();
        }
    }
}
