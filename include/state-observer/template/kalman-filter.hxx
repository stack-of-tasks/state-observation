template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::setB(const typename  KalmanFilter<n,m,p>::Bmatrix& B)
{
    b_.set(B,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::clearB()
{
    b_.reset();
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::setD(const typename  KalmanFilter<n,m,p>::Dmatrix& D)
{
    d_.set(D,0);
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::clearD()
{
    d_.reset();
}


template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::StateVector KalmanFilter<n,m,p>::prediction_(unsigned k)
{
    (void)k; //unused

    if (p>0)
        return this->a_()*this->x_()+this->b_()*this->u_[0]();
    else
        return this->a_()*this->x_();

}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::MeasureVector
KalmanFilter<n,m,p>::simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k)
{

    typename ObserverBase<n,m,p>::InputVector u ( ObserverBase<n,m,p>::InputVector::Zero());

    if (p>0)
    {
        if (d_()!=Dmatrix::Zero())
        {
            unsigned i;
            for (i=0; i<this->u_.size()&&this->u_[i].getTime()<k;++i)
            {
            }

            BOOST_ASSERT(i!=this->u_.size() && this->u_[i].getTime()<=k &&
                         "ERROR: The input feedthrough of the measurements is not set \
                         (the measurement at time k needs the input at time k which was not given) \
                         if you don't need the input in the computation of measurement, you \
                         must set D matrix to zero");
            u=this->u_[i]();
        }
        return this->c_()*x+this->d_()*u;
    }
    else
    {
        return this->c_()*x;
    }

}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::reset()
{
    KalmanFilterBase<n,m,p>::reset();

    b_.reset();
    d_.reset();
}

