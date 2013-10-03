#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    DiscreteTimeMatrix::DiscreteTimeMatrix(const Matrix& v,unsigned k):
            isSet_(true),
            k_(k),
            v_(v)
    {
    }

    DiscreteTimeMatrix::DiscreteTimeMatrix():
            isSet_(false),
            k_(0)
    {
    }


    ///Default constructor
    CheckedMatrix::CheckedMatrix():
            isSet_(false)
    {
    }

    ///A constructor with a given matrix value and a time index
    CheckedMatrix::CheckedMatrix(const Matrix& v):
            isSet_(false)
    {
    }



    ///Default constructor
    DiscreteTimeArray::DiscreteTimeArray():
        k_(0)
    {
    }

    std::vector<Matrix> DiscreteTimeArray::getArray() const
    {
        std::vector<Matrix> v;

        for (unsigned i=0;i<v_.size();++i)
        {
            v.push_back(v_[i]);
        }

        return v;
    }

    void DiscreteTimeArray::truncate(unsigned time)
    {
        for (unsigned i=v_.size(); (i > time - k_) && (i > 0) ;--i)
        {
            v_.pop_back();
        }
    }

}
