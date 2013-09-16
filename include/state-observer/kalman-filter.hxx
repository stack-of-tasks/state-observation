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
    return this->a_()*this->x_()+this->b_()*this->u_[0]();
}

template <unsigned n,unsigned m, unsigned p>
typename ObserverBase<n,m,p>::MeasureVector
KalmanFilter<n,m,p>::simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k)
{

    typename ObserverBase<n,m,p>::InputVector u= ObserverBase<n,m,p>::InputVector::Zero();


    if (d_()!=Dmatrix::Zero())
    {
        unsigned i;
        for (i=0; i<this->u_.size()&&this->u_[i].getTime()<k;++i)
        {
        }

        if (i==this->u_.size()||this->u_[i].getTime()>k)
            throw TimeException("The input feedthrough of the measurements is not set");
        else
            u=this->u_[i]();
    }
    return this->c_()*x+this->d_()*u;
}

template <unsigned n,unsigned m, unsigned p>
void KalmanFilter<n,m,p>::reset()
{
    KalmanFilterBase<n,m,p>::reset();

    b_.reset();
    d_.reset();
}

template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::StateVector KalmanFilter<n,m,0>::prediction_(unsigned k)
{
    return this->a_()*this->x_();
}

template <unsigned n,unsigned m>
typename ObserverBase<n,m,0>::MeasureVector
KalmanFilter<n,m,0>::simulateSensor_(const typename ObserverBase<n,m,0>::StateVector& x, unsigned k)
{
    return this->c_()*x;
}


