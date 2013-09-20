/**
 * \file      kalman-filter.hpp
 * \author    Mehdi Benallegue
 * \date       2012
 * \brief      Defines the class of a Linear Kalman filter
 *
 *              It implements the Kalman filter for linear systems (LTI-LTV).
 *              This is the class to instanciate when you want to use Kalman filtering
 *              for linear systems
 *
 *             x_{k+1}=A_k x_k+ B_k u_k + v_k
 *
 *             y_k=C_k x_k + D_k u_k + w_k
 *
 * \details
 *
 *
 */


#ifndef STATEOBSERVER_KALMANFILTERHPP
#define STATEOBSERVER_KALMANFILTERHPP

#include <state-observer/kalman-filter-base.hpp>

namespace observation
{

/**
     * \class  KalmanFilter
     * \brief
     *        The class of a Linear Kalman filter
     *
     *        It implements the Kalman filter for linear systems (LTI-LTV).
     *        This is the class to instanciate when you want to use Kalman filtering
     *        for linear systems. To use this class, one needs to provide the
     *        matrices that describe the dynamics of the state and the measurement.
     *
     *             x_{k+1}=A_k x_k+ B_k u_k + v_k
     *
     *             y_k=C_k x_k + D_k u_k + w_k
     *
     *         \li n : size of the state vector
     *         \li m : size of the measurements vector
     *         \li p : size of the input vector
     *
     * \details
     *
     */

    template <unsigned n,unsigned m, unsigned p=0>
    class KalmanFilter: public KalmanFilterBase<n,m,p>
    {
    public:
        /// The type of the matrix linking the input to the state
        typedef Eigen::Matrix<double, n,p> Bmatrix;

        /// The type of the matrix linking the input to the measurement
        typedef Eigen::Matrix<double, m,p> Dmatrix;

        /// Set the value of the input-state matrix
        virtual void setB(const Bmatrix& B);

        ///Clear the value of the input-state Matrix
        virtual void clearB();

        ///Set the value of the input-measurement matrix
        virtual void setD(const Dmatrix& D);

        ///Clear the value of the input-measurement matrix
        virtual void clearD();

        ///Reset all the observer
        virtual void reset();

    protected:
        /// The implementation of the (linear) prediction (state dynamics)
        virtual typename ObserverBase<n,m,p>::StateVector prediction_(unsigned k);

        /// The implementation of the (linear) measurement (state dynamics)
        virtual typename ObserverBase<n,m,p>::MeasureVector simulateSensor_(const typename ObserverBase<n,m,p>::StateVector& x, unsigned k);

        /// The type of the timed Input-State matrix
        typedef DiscreteTimeMatrix<n,p>  Bmat_;

        /// The type of the timed Input-Measurement matrix
        typedef DiscreteTimeMatrix<m,p>  Dmat_;

        /// The container of the Input-State matrix
        Dmat_ d_;

        /// The container of the Input-Measurement matrix
        Bmat_ b_;
    };


#include <state-observer/kalman-filter.hxx>

}

#endif //KALMANFILTERHPP
