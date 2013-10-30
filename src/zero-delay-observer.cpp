#include <state-observation/observer/zero-delay-observer.hpp>

namespace stateObservation
{

    void ZeroDelayObserver::setState
                            (const ObserverBase::StateVector& x_k,unsigned k)
    {
        BOOST_ASSERT(checkStateVector(x_k)
                            && "The size of the state vector is incorrect");
        x_.set(x_k,k);
        while (y_.size()>0 && y_.getFirstTime()<=k)
        {
            y_.popFront();
        }

        if (p_>0)
            while (u_.size()>0 && u_.getFirstTime()<k)
            {
                u_.popFront();
            }
    }

    void ZeroDelayObserver::clearStates()
    {
        x_.reset();
    }

    void ZeroDelayObserver::setMeasurement
                    (const ObserverBase::MeasureVector& y_k,unsigned k)
    {
        BOOST_ASSERT(checkMeasureVector(y_k)
                && "The size of the measure vector is incorrect");
        if (y_.size()>0)
            BOOST_ASSERT (y_.getLastTime()==k-1
                && "ERROR: The time is set incorrectly for \
                                the measurements (order or gap)");
        else
            BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k-1)
                && "ERROR: The time is set incorrectly for the measurements \
                                (must be [current_time+1])");

        y_.pushBack(y_k,k);
    }

    void ZeroDelayObserver::clearMeasurements()
    {
        y_.reset();
    }

    void ZeroDelayObserver::setInput
                    (const ObserverBase::InputVector& u_k,unsigned k)
    {
        if (p_>0)
        {
            BOOST_ASSERT(checkInputVector(u_k)
                        && "The size of the input vector is incorrect");

            if (u_.size()>0)
                BOOST_ASSERT (u_.getLastTime()==k-1
                        && "ERROR: The time is set incorrectly \
                                for the inputs (order or gap)");
            else
            {
                BOOST_ASSERT ( (!x_.isSet() || x_.getTime()==k || x_.getTime()==k-1)
                    && "ERROR: The time is set incorrectly for the \
                          inputs (must be [current_time] or [current_time+1])");
            }

            u_.pushBack(u_k,k);
        }
    }

    void ZeroDelayObserver::clearInputs()
    {
        if (p_>0)
            u_.reset();
    }

    ObserverBase::StateVector
    ZeroDelayObserver::getEstimateState(unsigned k)
    {
        unsigned k0=x_.getTime();

        BOOST_ASSERT(k0<=k
                && "ERROR: The observer cannot estimate previous states");

        for (unsigned i=k0;i<k;++i)
        {
            oneStepEstimation_();
            y_.popFront();
            if (p_>0)
                u_.popFront();
        }

        return x_();
    }


    unsigned ZeroDelayObserver::getCurrentTime()const
    {
        return x_.getTime();
    }

    unsigned ZeroDelayObserver::getInputTime()const
    {
        if (u_.size()>0)
        {
            return u_.getLastTime();
        }
        else
        {
            return 0;
        }
    }

    unsigned ZeroDelayObserver::getMeasurementTime()const
    {
        if (y_.size()>0)
        {
            return y_.getLastTime();
        }
        else
        {
            return 0;
        }
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
