#include <boost/random.hpp>

boost::lagged_fibonacci1279 Tools::gen_;

template <unsigned n>
Eigen::Matrix<double,n,1> Tools::getWGNoise(const Eigen::Matrix<double,n,n>& std,
	const Eigen::Matrix<double,n,1>& bias)
{
	boost::normal_distribution<> g(0, 1);
	Eigen::Matrix<double,n,1> ret;
	for (unsigned i=0;i<n;++i)
	{
		ret[i]=g(gen_);
	}
	ret=std*ret+bias;

	return ret;
}


