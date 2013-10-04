#include <state-observation/observer/kalman-filter.hpp>

namespace stateObservation
{

    void KalmanFilter::setB(const Bmatrix& B)
    {
        BOOST_ASSERT(checkBmatrix(B)&& "ERROR: The B matrix size is incorrect");
        b_.set(B);
    }

    void KalmanFilter::clearB()
    {
        b_.reset();
    }

    void KalmanFilter::setD(const Dmatrix& D)
    {
        BOOST_ASSERT(checkDmatrix(D)&& "ERROR: The D matrix size is incorrect");

        d_.set(D);
    }

    void KalmanFilter::clearD()
    {
        d_.reset();
    }


    ObserverBase::StateVector KalmanFilter::prediction_(unsigned k)
    {
        (void)k; //unused

        if (p_>0)
            return this->a_()*this->x_()+this->b_()*this->u_[k-1];
        else
            return this->a_()*this->x_();

    }

    ObserverBase::MeasureVector KalmanFilter::simulateSensor_(const StateVector& x, unsigned k)
    {

        InputVector u ( inputVectorZero());

        if (p_>0)
        {
            if (d_()!=getDmatrixZero())
            {
                BOOST_ASSERT(u_.checkIndex(k) &&
                             "ERROR: The input feedthrough of the measurements is not set \
                             (the measurement at time k needs the input at time k which was not given) \
                             if you don't need the input in the computation of measurement, you \
                             must set D matrix to zero");
                u=u_[k];
            }
            return c_()*x+d_()*u;
        }
        else
        {
            return c_()*x;
        }

    }

    void KalmanFilter::reset()
    {
        KalmanFilterBase::reset();

        b_.reset();
        d_.reset();
    }

    KalmanFilter::Bmatrix KalmanFilter::getBmatrixConstant(double c) const
    {
        return Bmatrix::Constant(n_,p_,c);
    }

    KalmanFilter::Bmatrix KalmanFilter::getBmatrixRandom() const
    {
        return Bmatrix::Random(n_,p_);
    }

    KalmanFilter::Bmatrix KalmanFilter::getBmatrixZero() const
    {
        return Bmatrix::Zero(n_,p_);
    }

    bool KalmanFilter::checkBmatrix(const Bmatrix & a) const
    {
        return (a.rows()==n_ && a.cols()==p_);
    }

    KalmanFilter::Dmatrix KalmanFilter::getDmatrixConstant(double c) const
    {
        return Dmatrix::Constant(m_,p_,c);
    }

    KalmanFilter::Dmatrix KalmanFilter::getDmatrixRandom() const
    {
        return Dmatrix::Random(m_,p_);
    }

    KalmanFilter::Dmatrix KalmanFilter::getDmatrixZero() const
    {
        return Dmatrix::Zero(m_,p_);
    }

    bool KalmanFilter::checkDmatrix(const Dmatrix & a) const
    {
        return (a.rows()==m_ && a.cols()==p_);
    }


    void KalmanFilter::setStateSize(unsigned n)
    {
        if (n!=n_)
        {
            KalmanFilterBase::setStateSize(n);
            b_.reset();
        }
    }

    void KalmanFilter::setMeasureSize(unsigned m)
    {
        if (m!=m_)
        {
            KalmanFilterBase::setMeasureSize(m);
            d_.reset();
        }
    }

    void KalmanFilter::setInputSize(unsigned p)
    {
        if (p!=p_)
        {
            KalmanFilterBase::setInputSize(p);
            b_.reset();
            d_.reset();
        }
    }

}
