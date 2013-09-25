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

namespace stateObserver
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

    class KalmanFilter: public KalmanFilterBase
    {
    public:

        KalmanFilter(){}

        KalmanFilter(unsigned n,unsigned m,unsigned p=0)
            :KalmanFilterBase(n,m,p){}

        /// The type of the matrix linking the input to the state
        typedef Eigen::MatrixXd Bmatrix;

        /// The type of the matrix linking the input to the measurement
        typedef Eigen::MatrixXd Dmatrix;

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

        /// The type of the jacobian df/dx
        Bmatrix getBmatrixConstant(double c) const;

        Bmatrix getBmatrixRandom() const;

        Bmatrix getBmatrixZero() const;

        bool checkBmatrix(const Bmatrix & ) const;


        /// The type of the jacobian dh/dx
        Dmatrix getDmatrixConstant(double c) const;

        Dmatrix getDmatrixRandom() const;

        Dmatrix getDmatrixZero() const;

        bool checkDmatrix(const Dmatrix &) const;


        virtual void setStateSize(unsigned n);

        virtual void setMeasureSize(unsigned m);

        virtual void setInputSize(unsigned p);

    protected:
        /// The implementation of the (linear) prediction (state dynamics)
        virtual StateVector prediction_(unsigned k);

        /// The implementation of the (linear) measurement (state dynamics)
        virtual MeasureVector simulateSensor_(const StateVector& x, unsigned k);

        /// The container of the Input-State matrix
        DiscreteTimeMatrix d_;

        /// The container of the Input-Measurement matrix
        DiscreteTimeMatrix b_;
    };


}

#endif //KALMANFILTERHPP
