#include <state-observation/observer/kalman-filter-base.hpp>

namespace stateObservation
{

    void KalmanFilterBase::setA(const Amatrix& A)
    {
        BOOST_ASSERT(checkAmatrix(A)&& "ERROR: The A matrix dimensions are wrong");
        a_=A;
    }

    void KalmanFilterBase::clearA()
    {
        a_.resize(0,0);
    }

    void KalmanFilterBase::setC( const Cmatrix& C)
    {
        BOOST_ASSERT(checkCmatrix(C)&& "ERROR: The C matrix dimensions are wrong");
        c_=C;
    }

    void KalmanFilterBase::clearC()
    {
        c_.resize(0,0);
    }


    void KalmanFilterBase::setR( const Rmatrix& R)
    {
        BOOST_ASSERT(checkRmatrix(R)&& "ERROR: The dimensions of the measurement noise covariance matrix R are wrong");
        r_=R;
    }

    void KalmanFilterBase::clearR()
    {
        r_.resize(0,0);
    }

    void KalmanFilterBase::setQ( const Qmatrix& Q)
    {
        BOOST_ASSERT(checkQmatrix(Q)&& "ERROR: The dimensions of the process noise covariance matrix Q are wrong");
        q_=Q;
    }

    void KalmanFilterBase::clearQ()
    {
        q_.resize(0,0);
    }

    void KalmanFilterBase::setStateCovariance(const Pmatrix& P)
    {
        BOOST_ASSERT(checkPmatrix(P)&& "ERROR: The P matrix dimensions are wrong");
        pr_=P;
    }


    void KalmanFilterBase::clearStateCovariance()
    {
        pr_.resize(0,0);
    }

    ObserverBase::StateVector KalmanFilterBase::oneStepEstimation_()
    {

        unsigned k=this->x_.getTime();
        BOOST_ASSERT(this->y_.size()> 0 && this->y_.checkIndex(k+1) && "ERROR: The measurement vector is not set");

        BOOST_ASSERT(checkAmatrix(a_) && "ERROR: The Matrix A is not initialized" );
        BOOST_ASSERT(checkCmatrix(c_) && "ERROR: The Matrix C is not initialized");
        BOOST_ASSERT(checkQmatrix(q_) && "ERROR: The Matrix Q is not initialized");
        BOOST_ASSERT(checkRmatrix(r_) && "ERROR: The Matrix R is not initialized");
        BOOST_ASSERT(checkPmatrix(pr_) && "ERROR: The Matrix P is not initialized");

        //prediction
        StateVector xbar=prediction_(k+1);
        Pmatrix pbar=a_*pr_*a_.transpose()+q_;

        //innovation
        MeasureVector ino= this->y_[k+1] - simulateSensor_(xbar,k+1);
        Rmatrix inoCov = c_ * pbar * c_.transpose() + r_;

        //gain
        Kmatrix kGain = (pbar * c_.transpose()) * inoCov.inverse();

        //update
        StateVector xhat=xbar+kGain*ino;

        this->x_.set(xhat,k+1);
        pr_=(getPmatrixIdentity()-kGain*c_)*pbar;

        return xhat;
    }

    KalmanFilterBase::Pmatrix KalmanFilterBase::getStateCovariance(unsigned k)
    {
        this->getEstimatedState(k);
        return pr_;
    }

    void KalmanFilterBase::reset()
    {
        ZeroDelayObserver::reset();

        clearA();
        clearC();
        clearQ();
        clearR();
        clearStateCovariance();
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
            clearA();
            clearC();
            clearQ();
            clearStateCovariance();
        }
    }

    void KalmanFilterBase::setMeasureSize(unsigned m)
    {
        if (m!=m_)
        {
            ZeroDelayObserver::setMeasureSize(m);
            clearC();
            clearR();
        }
    }
}
