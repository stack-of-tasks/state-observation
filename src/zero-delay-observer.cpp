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

        y_.setValue(y_k,k);

//        cout << "tata (" << y_k.rows() << "," << y_k.cols() << ")" << endl;
//        cout << "vecteur (" <<  getMeasurement(k).rows() << "," << getMeasurement(k).cols() << ")" << endl;

    }

    void ZeroDelayObserver::clearMeasurements()
    {
        y_.reset();
    }

    void ZeroDelayObserver::setInput
                    (const ObserverBase::InputVector& u_k,unsigned k)
    {

    // if(k%100==0)
       // cout << "\n\n\n\n\n ======> k: " << k << endl;

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

            u_.setValue(u_k,k);
        }

//        if (u_.size()>0 && y_.size()>0)
//       {
//            cout << "time setinput toto " << this->u_.getFirstTime() << ' ' << u_.getLastTime() << " " << k << endl;
//            cout << "input time: " << u_.getLastTime() << "\nmeasurement time: " << y_.getLastTime() << endl;
//       }
//       else
//       {
//            if(u_.size()<=0)
//                cout << "empty input" << endl;
//
//            if(y_.size()<=0)
//                cout << "empty  measurement" << endl;
//       }

    }

    void ZeroDelayObserver::clearInputs()
    {
        if (p_>0)
            u_.reset();
    }

    ObserverBase::StateVector
    ZeroDelayObserver::getEstimatedState(unsigned k)
    {
        unsigned k0=x_.getTime();

        BOOST_ASSERT(k0<=k
                && "ERROR: The observer cannot estimate previous states");

        for (unsigned i=k0;i<k;++i)
        {
            oneStepEstimation_();
            if (y_.getFirstTime()<k)
                y_.popFront();

            if (p_>0)
                if (u_.getFirstTime()<k)
                    u_.popFront();
        }

       // std::cout << "x_ getEstsimatedState " << x_() << std::endl;

        return x_();
    }


    unsigned ZeroDelayObserver::getCurrentTime()const
    {
        return x_.getTime();
    }

    Vector ZeroDelayObserver::getInput(unsigned k) const
    {
        return u_[k];
    }

    unsigned ZeroDelayObserver::getInputsNumber()const
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

    Vector ZeroDelayObserver::getMeasurement(unsigned k) const
    {
        return y_[k];
    }

    unsigned ZeroDelayObserver::getMeasurementTime()const
    {
        BOOST_ASSERT(y_.size()>0
                && "ERROR: There is no measurements registered (past measurements are erased)");
        return y_.getLastTime();
    }

    unsigned ZeroDelayObserver::getMeasurementsNumber()const
    {

//        cout << "y_.size(): " << y_.size() << endl;
//        Vector measure(getMeasurement(y_.getLastTime()));
//        cout << "y_: (" << measure.rows() << "," << measure.cols() << ") " << endl;

        return y_.size();
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
