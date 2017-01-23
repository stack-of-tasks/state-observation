#ifndef STATEOBSERVATIONLOGGER_H
#define STATEOBSERVATIONLOGGER_H
#include <map>
#include <string>

#include <state-observation/tools/definitions.hpp>

namespace stateObservation
{
  namespace tools
  {
    class Logger
    {
    public:
      Logger();
      virtual ~Logger();

      ///Use this function to start the recoding of this variable.
      /// WARNING: Be sure that the recorded variable keeps the same memory address
       void record(const Matrix & matrix, const std::string & filename=std::string(""));
      //The use of two different overloads for this function serves two purposes
      // 1- have a different logger for doubles and floats (obvious)
      // 2- avoid implicit conversions from another floating type into a double or a float.
      //    by making the conversion ambiguous.
      void record(const double & d, const std::string & filename=std::string(""));
      void record(const float & d, const std::string & filename=std::string(""));

      void record(const int & i, const std::string & filename=std::string(""));
      void record(const long & i, const std::string & filename=std::string(""));

      ///set the Path for the log files, the filenames will be appended to this path
      void setPath(const std::string & path);

      ///update the log with a new value of
      void insert(const Matrix & matrix);
      void insert(const double & d);
      void insert(const float & d);
      void insert(const int & i);
      void insert(const long & i);

      void update();

      const IndexedMatrixArray & getRecord(const Matrix & matrix) const;

      IndexedMatrixArray & getRecord(const Matrix & matrix);

      void save();

      void clear();

    protected:
      struct log_s
      {

        IndexedMatrixArray array;
        std::string filename;
        enum DataType {typematrix,typedouble, typefloat, typeint, typelong} type;

        log_s(const std::string & newfilename, DataType t)
        {
          filename = newfilename;
          type =t;
        }
      };
      std::string path_;
      std::map<const void *, log_s > logs_;
      Matrix scalar_;


    private:

    };
  }
}

#endif // STATEOBSERVATIONLOGGER_H
