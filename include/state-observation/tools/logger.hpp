#ifndef STATEOBSERVATIONLOGGER_H
#define STATEOBSERVATIONLOGGER_H

#include <typeinfo>
#include <map>
#include <string>
#include <stdexcept>

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
      ///otherwise use updateAddress
      template <typename T>
        void record(const T * address, const std::string & filename=std::string(""));
      template <typename T>
        void record(const T & reference, const std::string & filename=std::string(""));

      ///updates the address of a recorded variable with a new address
      template <typename T>
        void updateAddress(const void * oldAddress, const T * newAddress);


      ///set the Path for the log files, the filenames will be appended to this path
      void setPath(const std::string & path);

      ///update the log with a new value of the reference
      template <typename T>
      void push(const T & reference);

      template <typename T>
      void push(const T * address);

      ///updates all the logs for all recorded variables
      void push();

      const IndexedMatrixArray & getRecord(const Matrix & matrix) const;

      IndexedMatrixArray & getRecord(const Matrix & matrix);

      void save();

      void clear();

    protected:
      struct log_s
      {

        IndexedMatrixArray array;
        std::string filename;
        const std::type_info * type;

        log_s(const std::string & newfilename):
          type (0x0)
        {
          filename = newfilename;

        }

        template <typename T>
        void setType()
        {
          type =&typeid(T);
        }
      };

      typedef std::map<const void *, log_s > Tmap;
      typedef std::pair<const void *, log_s > Tpair;

      void update_(const Tmap::iterator & i);

      std::string path_;
      std::map<const void *, log_s > logs_;
      Matrix scalar_;


    private:

    };
  }
}

#include <state-observation/tools/logger.hxx>

#endif // STATEOBSERVATIONLOGGER_H
