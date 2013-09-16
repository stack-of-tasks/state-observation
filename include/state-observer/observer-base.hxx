template<unsigned r,unsigned c>
DiscreteTimeMatrix<r, c>::DiscreteTimeMatrix(const MatrixT& v,unsigned k):
        isSet_(true),
        k_(k),
        v_(v)
{
}

template<unsigned r,unsigned c>
DiscreteTimeMatrix<r, c>::DiscreteTimeMatrix():
        isSet_(false),
        k_(0)
{
}

template<unsigned r,unsigned c>
void DiscreteTimeMatrix<r, c>::set(const MatrixT& v,unsigned k)
{
    k_=k;
    v_=v;
    isSet_=true;
}

template<unsigned r,unsigned c>
void DiscreteTimeMatrix<r, c>::reset()
{
    k_=0;
    isSet_=false;
}

template<unsigned r,unsigned c>
typename DiscreteTimeMatrix<r, c>::MatrixT DiscreteTimeMatrix<r, c>::operator()()const
{
    if (isSet_)
        return v_;
    else
        throw InitializationException("Vector not initialized");
}

template<unsigned r,unsigned c>
const unsigned & DiscreteTimeMatrix<r, c>::getTime()const
{
    if (isSet_)
        return k_;
    else
        throw InitializationException("Vector not initialized");
}


template<unsigned r,unsigned c>
const bool & DiscreteTimeMatrix<r, c>::isSet()const
{
    return isSet_;
}

template <unsigned n,unsigned m, unsigned p>
void ObserverBase<n,m,p>::reset()
{
    clearState();
    clearMeasurements();
    clearInputs();
}
