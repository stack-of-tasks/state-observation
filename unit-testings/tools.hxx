#include <boost/random.hpp>

boost::lagged_fibonacci1279 Tools::gen_;

Eigen::MatrixXd Tools::getWGNoise(const Eigen::MatrixXd& std,
	const Eigen::MatrixXd& bias,unsigned rows, unsigned cols)
{
	boost::normal_distribution<> g(0, 1);
	Eigen::MatrixXd ret= Eigen::MatrixXd::Zero(rows,cols);
	for (unsigned i=0;i<rows;++i)
	{
	    for (unsigned j=0; j<cols; ++j)
            ret(i,j)=g(gen_);
	}
	ret=std*ret+bias;

	return ret;
}


