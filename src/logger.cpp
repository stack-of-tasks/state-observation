#include <fstream>
#include <stdexcept>
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

    void Logger::record(const Matrix& matrix, const std::string& filename)
    {
      logs_.insert(std::pair<const void *, log_s>(&matrix,log_s(filename,log_s::typematrix)));
    }

    void Logger::record(const double& scalar, const std::string& filename)
    {
      logs_.insert(std::pair<const void *, log_s>(&scalar,log_s(filename,log_s::typedouble)));
    }

    void Logger::record(const float& scalar, const std::string& filename)
    {
      logs_.insert(std::pair<const void *, log_s>(&scalar,log_s(filename,log_s::typefloat)));
    }

    void Logger::record(const int& integer, const std::string& filename)
    {
      logs_.insert(std::pair<const void *, log_s>(&integer,log_s(filename,log_s::typeint)));
    }

    void Logger::record(const long& integer, const std::string& filename)
    {
      logs_.insert(std::pair<const void *, log_s>(&integer,log_s(filename,log_s::typelong)));
    }



    void Logger::setPath(const std::string& path)
    {
      path_=path;

    }

    void Logger::insert(const Matrix & matrix)
    {
      std::map<const void *, log_s >::iterator i= logs_.find( & matrix );

      if (i!=logs_.end())
        i->second.array.pushBack(matrix);
      else
        throw std::invalid_argument("The logger cannot find the data, please use record function");
    }

    void Logger::save()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        if (i->second.filename!=std::string(""))
        {
          i->second.array.writeInFile((path_+i->second.filename).c_str());
        }
      }
    }

    void Logger::clear()
    {
      logs_.clear();
    }

    const IndexedMatrixArray & Logger::getRecord(const Matrix& matrix) const
    {
      return logs_.find(& matrix)->second.array;
    }

    IndexedMatrixArray & Logger::getRecord(const Matrix& matrix)
    {
      return logs_.find(& matrix)->second.array;
    }

    void Logger::update()
    {
      for (std::map<const void *, log_s >::iterator i=logs_.begin();i!=logs_.end();++i)
      {
        switch (i->second.type)
        {
        case log_s::typematrix:
          i->second.array.pushBack(*static_cast<const Matrix *>(i->first));
          return;
        case log_s::typeint:
            scalar_(1,1)=*static_cast<const int *>(i->first);
          break;
        case log_s::typedouble:
            scalar_(1,1)=*static_cast<const double *>(i->first);
          break;
        }
        //this line is not executed in case i->second.type=log_s::typematrix:
        i->second.array.pushBack(scalar_);


      }
    }

  }
}
