#ifndef EXTENDEDKALMANFILTERHPP
#define EXTENDEDKALMANFILTERHPP

#include <state-observer/kalman-filter-base.hpp>

namespace observation
{

    //n : size of the state vector
    //m : size of the measurements vector
    //p : size of the input vector
    template <unsigned n,unsigned m, unsigned p=0>
    class ExtendedKalmanFilter: public KalmanFilterBase<n,m,p>
    {
    public:
        class DynamicsFunctorBase
        {
        public:
            virtual typename ObserverBase<n,m,p>::StateVector stateDynamics
            (const typename ObserverBase<n,m,p>::StateVector& x,
             const typename ObserverBase<n,m,p>::InputVector& u,
             unsigned k)=0;

            virtual typename ObserverBase<n,m,p>::MeasureVector measureDynamics
            (const typename ObserverBase<n,m,p>::StateVector& x,
             const typename ObserverBase<n,m,p>::InputVector& u,
             unsigned k)=0;

            virtual void reset(){}

        };

        ExtendedKalmanFilter(bool directInputOutputFeedthrough=0):
                directInputOutputFeedthrough_(directInputOutputFeedthrough),f_(0x0)
        {
        }


        void setFunctor(DynamicsFunctorBase* f);

        void clearFunctor();

        void setDirectInputOutputFeedthrough(bool b=true);

        virtual typename ObserverBase<n,m,p>::StateVector getPrediction(unsigned k);

        //use the finite differences method to get the A and C matrices using
        //finite difference method (the forward difference method)
        //dx is the step vector
        //x is the linearization center
        virtual typename KalmanFilterBase<n,m,p>::Amatrix getAMatrixFD(const typename ObserverBase<n,m,p>::StateVector &dx);
        virtual typename KalmanFilterBase<n,m,p>::Cmatrix getCMatrixFD(const typename ObserverBase<n,m,p>::StateVector &dx);

        virtual void reset();

    protected:
        virtual typename ObserverBase<n,m,p>::StateVector prediction_(unsigned k);
        virtual typename ObserverBase<n,m,p>::MeasureVector simulateSensor_
        (const typename ObserverBase<n,m,p>::StateVector& x, unsigned k);

        typename ObserverBase<n,m,p>::State xbar_;

        bool directInputOutputFeedthrough_;

        DynamicsFunctorBase* f_;
    };

    template <unsigned n,unsigned m>
    class ExtendedKalmanFilter<n,m,0>: public KalmanFilterBase<n,m,0>
    {
    public:
        class DynamicsFunctorBase
        {
        public:
            virtual typename ObserverBase<n,m,0>::StateVector stateDynamics(const typename ObserverBase<n,m,0>::StateVector& x, unsigned k)=0;
            virtual typename ObserverBase<n,m,0>::MeasureVector measureDynamics(const typename ObserverBase<n,m,0>::StateVector& x, unsigned k)=0;
            virtual void reset();
        };

        ExtendedKalmanFilter():f_(0x0)
        {}

        void setFunctor(DynamicsFunctorBase* f);

        void clearFunctor();

        virtual typename ObserverBase<n,m,0>::StateVector getPrediction(unsigned k);

        //use the finite differences method to get the A and C matrices using
        //finite difference method (the forward difference method)
        //dx is the step vector
        //x is the linearization center
        virtual typename KalmanFilterBase<n,m,0>::Amatrix getAMatrixFD(const typename ObserverBase<n,m,0>::StateVector &dx);
        virtual typename KalmanFilterBase<n,m,0>::Cmatrix getCMatrixFD(const typename ObserverBase<n,m,0>::StateVector &dx);

        virtual void reset();

    protected:
        virtual typename ObserverBase<n,m,0>::StateVector prediction_(unsigned k);
        virtual typename ObserverBase<n,m,0>::MeasureVector simulateSensor_(const typename ObserverBase<n,m,0>::StateVector& x, unsigned k);

        const typename ObserverBase<n,m,0>::State xbar_;
        DynamicsFunctorBase* f_;
    };

#include <state-observer/extended-kalman-filter.hxx>

}

#endif //EXTENDEDKALMANFILTERHPP
