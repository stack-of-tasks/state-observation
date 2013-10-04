///Set the value of the matrix and the time sample
void CheckedMatrix::set(const Matrix& v)
{
    v_=v;
    isSet_=true;
}

///Get the matrix value
Matrix CheckedMatrix::operator()()const
{
    check_();
    return v_;
}

///Says whether the matrix is initialized or not
const bool & CheckedMatrix::isSet()const
{
    return isSet_;
}

///Switch off the initalization flag, the value is no longer accessible
void CheckedMatrix::reset()
{
    isSet_=false;
}


///Checks whether the matrix is set or not (assert)
///does nothing in release mode
void CheckedMatrix::check_()const
{
    BOOST_ASSERT(isSet_ && "Error: Matrix not initialized");
}

///Set the value of the matrix and the time sample
inline void DiscreteTimeMatrix::set(const Matrix& v,unsigned k)
{
    k_=k;
    v_=v;
    isSet_=true;
}

///Get the matrix value
Matrix DiscreteTimeMatrix::operator()()const
{
    check_();
    return v_;
}

///Get the time index
const unsigned & DiscreteTimeMatrix::getTime()const
{
    check_();
    return k_;
}

///Says whether the matrix is initialized or not
const bool & DiscreteTimeMatrix::isSet()const
{
    return isSet_;
}

///Switch off the initalization flag, the value is no longer accessible
void DiscreteTimeMatrix::reset()
{
    k_=0;
    isSet_=false;
}

///Checks whether the matrix is set or not (assert)
///does nothing in release mode
void DiscreteTimeMatrix::check_()const
{
    BOOST_ASSERT(isSet_ && "Error: Matrix not initialized");
}

///Set the value of the matrix and the time sample
void DiscreteTimeArray::pushBack(const Matrix& v,unsigned k)
{
    checkNext_(k);
    if (v_.size()==0)
        k_=k;

    v_.push_back(v);
}

void DiscreteTimeArray::popFront()
{
    check_();
    v_.pop_front();
    ++k_;
}

///Get the matrix value
Matrix DiscreteTimeArray::operator[](unsigned time)const
{
    check_(time);
    return v_[time - k_];
}

///Get the matrix value
Matrix & DiscreteTimeArray::operator[](unsigned time)
{
    check_(time);
    return v_[time - k_];
}

///Get the time index
unsigned DiscreteTimeArray::getLastTime()const
{
    check_();
    return k_+v_.size()-1;
}

///Get the time index
unsigned DiscreteTimeArray::getFirstTime()const
{
    check_();
    return k_;
}

unsigned DiscreteTimeArray::size() const
{
    return v_.size();
}

///Switch off the initalization flag, the value is no longer accessible
void DiscreteTimeArray::reset()
{
    k_=0;
    v_.clear();
}


bool DiscreteTimeArray::checkIndex(unsigned time) const
{
    return (v_.size()>0 && k_<=time && k_+v_.size() > time);
}


///Checks whether the matrix is set or not (assert)
///does nothing in release mode
void DiscreteTimeArray::check_(unsigned time)const
{
    BOOST_ASSERT(checkIndex(time) && "Error: Time out of range");
}


///Checks whether the matrix is set or not (assert)
///does nothing in release mode
void DiscreteTimeArray::check_()const
{
    BOOST_ASSERT(v_.size() && "Error: Matrix array not initialized");
}

void DiscreteTimeArray::checkNext_(unsigned time)const
{
    BOOST_ASSERT( (v_.size()==0 || k_+v_.size() == time )&& "Error: Time instant must be consecutive");
}

