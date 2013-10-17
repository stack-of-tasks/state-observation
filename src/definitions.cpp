#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
    DiscreteTimeMatrix::DiscreteTimeMatrix(const Matrix& v,unsigned k):
            k_(k),
            v_(v)
    {
    }

    DiscreteTimeMatrix::DiscreteTimeMatrix():
            k_(0)
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
        if (v_.size()>0)
        {
            if (time > getFirstTime())
            {
                for (unsigned i=getLastTime(); i>=time ;--i)
                {
                    v_.pop_back();
                }
            }
            else
            {
                v_.clear();
            }
        }
    }

}
