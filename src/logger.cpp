#include <fstream>
#include <state-observation/tools/logger.hpp>


namespace stateObservation
{
  namespace tools
  {
    Logger::Logger()
    {
      scalar_.resize(1,1);
    }

    Logger::~Logger()
    {
      //dtor
    }

     void Logger::setPath(const std::string& path)
    {
      path_=path;

    }

    void Logger::save()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        if (i->second.filename!=std::string(""))
        {
          i->second.array.writeInFile((path_+std::string("/")+i->second.filename).c_str());
        }
      }
    }

    void Logger::clear()
    {
      logs_.clear();
    }

    const IndexedMatrixArray & Logger::getRecord(const Matrix& matrix) const
    {
      Tmap::const_iterator i= logs_.find(& matrix);
      if (i==logs_.end())
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
      else
        return i->second.array;
    }

    IndexedMatrixArray & Logger::getRecord(const Matrix& matrix)
    {
      Tmap::iterator i= logs_.find(& matrix);
      if (i==logs_.end())
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
      else
        return i->second.array;
    }

    void Logger::push()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        update_(i);
      }
    }

    void Logger::update_(const Tmap::iterator & i)
    {
      if (*i->second.type == typeid(Matrix))
      {
        i->second.array.pushBack(*static_cast<const Matrix *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Vector))
      {
        i->second.array.pushBack(*static_cast<const Vector *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Matrix3))
      {
        i->second.array.pushBack(*static_cast<const Matrix3 *>(i->first));
        return;
      }
      if (*i->second.type == typeid(Vector3))
      {
        i->second.array.pushBack(*static_cast<const Vector3 *>(i->first));
        return;
      }
      else if (*i->second.type == typeid(int))
      {
        scalar_(0,0)=*static_cast<const int *>(i->first);
      }
      else if (*i->second.type == typeid(unsigned))
      {
        scalar_(0,0)=*static_cast<const unsigned *>(i->first);
      }
      else if (*i->second.type == typeid(long))
      {
        scalar_(0,0)=*static_cast<const long *>(i->first);
      }
      else if (*i->second.type == typeid(double))
      {
        scalar_(0,0)=*static_cast<const double *>(i->first);
      }
      else if (*i->second.type == typeid(float))
      {
        scalar_(0,0)=*static_cast<const float *>(i->first);
      }



      i->second.array.pushBack(scalar_);
    }

  }
}
