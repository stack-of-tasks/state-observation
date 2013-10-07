#include <state-observation/dynamical-system/dynamical-system-functor-base.hpp>

namespace stateObservation
{

    DynamicalSystemFunctorBase::DynamicalSystemFunctorBase()
    {
        //ctor
    }

    DynamicalSystemFunctorBase::~DynamicalSystemFunctorBase()
    {
        //dtor
    }

    bool DynamicalSystemFunctorBase::checkStateVector(const Vector & v)
    {
        return (v.rows()==getStateSize() && v.cols()==1);
    }

    bool DynamicalSystemFunctorBase::checkInputvector(const Vector & v)
    {
        return (v.rows()==getInputSize() && v.cols()==1);
    }
}
