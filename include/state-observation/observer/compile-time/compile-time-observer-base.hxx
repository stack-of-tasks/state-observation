template<unsigned r,unsigned c>
IndexedMatrix<r, c>::IndexedMatrix(const MatrixT& v,unsigned k):
        isSet_(true),
        k_(k),
        v_(v)
{
}

template<unsigned r,unsigned c>
IndexedMatrix<r, c>::IndexedMatrix():
        isSet_(false),
        k_(0)
{
}

template<unsigned r,unsigned c>
void IndexedMatrix<r, c>::set(const MatrixT& v,unsigned k)
{
    k_=k;
    v_=v;
    isSet_=true;
}

template<unsigned r,unsigned c>
void IndexedMatrix<r, c>::reset()
{
    k_=0;
    isSet_=false;
}

template<unsigned r,unsigned c>
typename IndexedMatrix<r, c>::MatrixT IndexedMatrix<r, c>::operator()()const
{
    check_();
    return v_;
}

template<unsigned r,unsigned c>
const unsigned & IndexedMatrix<r, c>::getTime()const
{
    check_();
    return k_;
}

template<unsigned r,unsigned c>
const bool & IndexedMatrix<r, c>::isSet()const
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

template<unsigned r,unsigned c>
void IndexedMatrix<r, c>::check_()const
{
    BOOST_ASSERT(isSet_ && "Error : Matrix not initialized");
}
