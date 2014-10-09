#include <fstream>

#include <state-observation/tools/definitions.hpp>


namespace stateObservation
{
    IndexedMatrix::IndexedMatrix(const Matrix& v,unsigned k):
            k_(k),
            v_(v)
    {
    }

    IndexedMatrix::IndexedMatrix():
            k_(0),
            v_(Matrix::Zero(0,0))
    {
    }


    ///Default constructor
    IndexedMatrixArray::IndexedMatrixArray():
            k_(0)
    {
    }

    std::vector<Matrix> IndexedMatrixArray::getArray() const
    {
        std::vector<Matrix> v;

        for (unsigned i=0;i<v_.size();++i)
        {
            v.push_back(v_[i]);
        }

        return v;
    }

    void IndexedMatrixArray::truncate(unsigned time)
    {
        if (v_.size()>0)
        {
            if (time > getFirstIndex())
            {
                for (unsigned i=getLastIndex(); i>=time ;--i)
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

    void IndexedMatrixArray::getFromFile(char * filename , size_t rows, size_t cols)
    {
        reset();

        std::ifstream f;

        f.open(filename);

        Matrix m(Matrix::Zero(rows,cols));

        bool continuation=true;

        while (continuation)
        {
            unsigned k;

            f >> k;

            if (f.eof())
                continuation=false;
            else
            {
                for (size_t i = 0 ; i<rows; ++i)
                {
                    for (size_t j = 0 ; j<cols; ++j)
                    {
                        f >> m(i,j);
                    }
                }

                setValue(m,k);
            }
        }
    }


    void IndexedMatrixArray::writeInFile(char * filename)
    {
    	std::ofstream f;

        f.open(filename);

        for (size_t k=getFirstIndex();k<=getLastIndex();++k)
        {

            f << k;

            Matrix & m = operator[](k);

            for (size_t i = 0 ; i< m.rows(); ++i)
            {
            	for (size_t j = 0 ; j< m.cols(); ++j)
                {
            		f << " "<< m(i,j);
                }
            }

            f << std::endl;
        }
    }

}
