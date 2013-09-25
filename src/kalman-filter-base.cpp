#include <state-observer/kalman-filter-base.hpp>

namespace stateObserver
{

    void KalmanFilterBase::setA(const Amatrix& A)
    {
        BOOST_ASSERT(checkAmatrix(A)&& "ERROR: The A matrix dimensions are wrong");
        a_.set(A,0);
    }

    void KalmanFilterBase::clearA()
    {
        a_.reset();
    }

    void KalmanFilterBase::setC( const Cmatrix& C)
    {
        BOOST_ASSERT(checkCmatrix(C)&& "ERROR: The C matrix dimensions are wrong");
        c_.set(C,0);
    }

    void KalmanFilterBase::clearC()
    {
        c_.reset();
    }


    void KalmanFilterBase::setR( const Rmatrix& R)
    {
        BOOST_ASSERT(checkRmatrix(R)&& "ERROR: The R matrix dimensions are wrong");
        r_.set(R,0);
    }

    void KalmanFilterBase::clearR()
    {
        r_.reset();
    }

    void KalmanFilterBase::setQ( const Qmatrix& Q)
    {
        BOOST_ASSERT(checkQmatrix(Q)&& "ERROR: The Q matrix dimensions are wrong");
        q_.set(Q,0);
    }

    void KalmanFilterBase::clearQ()
    {
        q_.reset();
    }

    void KalmanFilterBase::setStateCovariance(const Pmatrix& P)
    {
        BOOST_ASSERT(checkPmatrix(P)&& "ERROR: The P matrix dimensions are wrong");
        pr_.set(P,this->x_.getTime());
    }


    void KalmanFilterBase::clearStateCovariance()
    {
        pr_.reset();
    }

    ObserverBase::StateVector KalmanFilterBase::oneStepEstimation_()
    {
        unsigned k=this->x_.getTime();
        BOOST_ASSERT(this->y_.size()> 0 && this->y_[0].getTime()==k+1 && "ERROR: The measurement vector is not set");
        if (p_>0)
            BOOST_ASSERT(this->u_.size()> 0 && this->u_[0].getTime()==k && "ERROR: The input vector is not set");

        BOOST_ASSERT(a_.isSet() && "ERROR: The Matrix A is not initialized" );
        BOOST_ASSERT(c_.isSet() && "ERROR: The Matrix C is not initialized");
        BOOST_ASSERT(q_.isSet() && "ERROR: The Matrix Q is not initialized");
        BOOST_ASSERT(r_.isSet() && "ERROR: The Matrix R is not initialized");
        BOOST_ASSERT(pr_.isSet() && "ERROR: The Matrix P is not initialized");

        Amatrix a=a_();
        Cmatrix c=c_();
        Pmatrix px=pr_();

        //prediction
        StateVector xbar=prediction_(k+1);
        Pmatrix pbar=a*px*a.transpose()+q_();

        //innovation
        MeasureVector ino= this->y_[0]() - simulateSensor_(xbar,k+1);
        Rmatrix inoCov = c * pbar * c.transpose() + r_();

        //gain
        Kmatrix kGain = (pbar * c.transpose()) * inoCov.inverse();

        //update
        StateVector xhat=xbar+kGain*ino;

        this->x_.set(xhat,k+1);
        pr_.set((getPmatrixIdentity()-kGain*c)*pbar,k+1);

        return xhat;
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getStateCovariance(unsigned k)
    {
        this->getEstimateState(k);
        return pr_();
    }

    void KalmanFilterBase::reset()
    {
        ZeroDelayObserver::reset();

        a_.reset();
        c_.reset();
        a_.reset();
        c_.reset();
        q_.reset();
        r_.reset();
        pr_.reset();
    }


    KalmanFilterBase::Amatrix KalmanFilterBase::getAmatrixConstant(double c) const
    {
        return Amatrix::Constant(n_,n_,c);
    }

    KalmanFilterBase::Amatrix KalmanFilterBase::getAmatrixRandom() const
    {
        return Amatrix::Random(n_,n_);
    }

    KalmanFilterBase::Amatrix KalmanFilterBase::getAmatrixZero() const
    {
        return Amatrix::Zero(n_,n_);
    }

    KalmanFilterBase::Amatrix KalmanFilterBase::getAmatrixIdentity() const
    {
        return Amatrix::Identity(n_,n_);
    }

    bool KalmanFilterBase::checkAmatrix(const Amatrix & a) const
    {
        return (a.rows()==n_ && a.cols()==n_);
    }

    KalmanFilterBase::Cmatrix KalmanFilterBase::getCmatrixConstant(double c) const
    {
        return Cmatrix::Constant(m_,n_,c);
    }

    KalmanFilterBase::Cmatrix KalmanFilterBase::getCmatrixRandom() const
    {
        return Cmatrix::Random(m_,n_);
    }

    KalmanFilterBase::Cmatrix KalmanFilterBase::getCmatrixZero() const
    {
        return Cmatrix::Zero(m_,n_);
    }

    bool KalmanFilterBase::checkCmatrix(const Cmatrix & a) const
    {
        return (a.rows()==m_ && a.cols()==n_);
    }

    KalmanFilterBase::Qmatrix KalmanFilterBase::getQmatrixConstant(double c) const
    {
        return Qmatrix::Constant(n_,n_,c);
    }

    KalmanFilterBase::Qmatrix KalmanFilterBase::getQmatrixRandom() const
    {
        return Qmatrix::Random(n_,n_);
    }

    KalmanFilterBase::Qmatrix KalmanFilterBase::getQmatrixZero() const
    {
        return Qmatrix::Zero(n_,n_);
    }

    KalmanFilterBase::Qmatrix KalmanFilterBase::getQmatrixIdentity() const
    {
        return Qmatrix::Identity(n_,n_);
    }

    bool KalmanFilterBase::checkQmatrix(const Qmatrix & a) const
    {
        return (a.rows()==n_ && a.cols()==n_);
    }

    KalmanFilterBase::Rmatrix KalmanFilterBase::getRmatrixConstant(double c) const
    {
        return Cmatrix::Constant(m_,m_,c);
    }

    KalmanFilterBase::Rmatrix KalmanFilterBase::getRmatrixRandom() const
    {
        return Cmatrix::Random(m_,m_);
    }

    KalmanFilterBase::Rmatrix KalmanFilterBase::getRmatrixZero() const
    {
        return Rmatrix::Zero(m_,m_);
    }

    KalmanFilterBase::Rmatrix KalmanFilterBase::getRmatrixIdentity() const
    {
        return Rmatrix::Identity(m_,m_);
    }

    bool KalmanFilterBase::checkRmatrix(const Rmatrix & a) const
    {
        return (a.rows()==m_ && a.cols()==m_);
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getPmatrixConstant(double c) const
    {
        return Pmatrix::Constant(n_,n_,c);
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getPmatrixRandom() const
    {
        return Pmatrix::Random(n_,n_);
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getPmatrixZero() const
    {
        return Pmatrix::Zero(n_,n_);
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getPmatrixIdentity() const
    {
        return Pmatrix::Identity(n_,n_);
    }

    bool KalmanFilterBase::checkPmatrix(const Pmatrix & a) const
    {
        return (a.rows()==n_ && a.cols()==n_);
    }

    void KalmanFilterBase::setStateSize(unsigned n)
    {
        if (n!=n_)
        {
            ZeroDelayObserver::setStateSize(n);
            a_.reset();
            c_.reset();
            q_.reset();
            pr_.reset();
        }
    }

    void KalmanFilterBase::setMeasureSize(unsigned m)
    {
        if (m!=m_)
        {
            ZeroDelayObserver::setMeasureSize(m);
            c_.reset();
            r_.reset();
        }
    }
}
