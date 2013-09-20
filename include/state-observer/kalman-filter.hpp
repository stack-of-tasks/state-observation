#ifndef KALMANFILTERHPP
#define KALMANFILTERHPP

#include <state-observer/kalman-filter-base.hpp>

namespace observation
{

    //n : size of the state vector
    //m : size of the measurements vector
    //p : size of the input vector
    template <unsigned n,unsigned m, unsigned p=0>
    class KalmanFilter: public KalmanFilterBase<n,m,p>
    {
    public:
        typedef Eigen::Matrix<double, n,p> Bmatrix;
        typedef Eigen::Matrix<double, m,p> Dmatrix;

        virtual void setB(const Bmatrix& B);
        virtual void clearB();

        virtual void setD(const Dmatrix& D);
        virtual void clearD();

        virtual void reset();

    protected:
        virtual typename ObserverBase<n,m,p>::StateVector prediction_(unsigned k);
        virtual typename ObserverBase<n,m,p>::MeasureVector simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k);

        typedef DiscreteTimeMatrix<n,p>  Bmat_;
        typedef DiscreteTimeMatrix<m,p>  Dmat_;

        Dmat_ d_;
        Bmat_ b_;
    };


#include <state-observer/kalman-filter.hxx>

}

#endif //KALMANFILTERHPP
