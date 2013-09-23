#ifndef OBSERVATIONTOOLSHPP
#define OBSERVATIONTOOLSHPP

#include <Eigen/Core>
#include <boost/random.hpp>

namespace stateObserver
{
    namespace unitTesting
    {
        class Tools
        {
        public:

            //add White Gaussian Noise to a vector
            //having a given bias and standard deviation(std)
            template <unsigned n>
            static Eigen::Matrix<double,n,1> getWGNoise( const Eigen::Matrix<double,n,n>& std,
                    const Eigen::Matrix<double,n,1>& bias);

        protected:
            static boost::lagged_fibonacci1279 gen_;

        };

#       include "tools.hxx"
    }
}



#endif //OBSERVATIONTOOLSHPP
