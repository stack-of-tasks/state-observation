
namespace stateObservation
{
  namespace tools
  {
    template <typename T>
    void Logger::record(const T * address, const std::string& filename)
    {
      log_s log(filename);
      log.setType<T>();
      if    (*log.type!=typeid(Matrix)
          && *log.type!=typeid(Matrix3)
          && *log.type!=typeid(Vector)
          && *log.type!=typeid(Vector3)
          && *log.type!=typeid(int)
          && *log.type!=typeid(long)
          && *log.type!=typeid(double)
          && *log.type!=typeid(float))
      {
        throw std::invalid_argument
                ((std::string("The type ") +
                    std::string("log.type->name")
                      +std::string(" is not supported")).c_str());
        BOOST_ASSERT(false && "Tentative to record an unsupported type");
      }

      logs_.insert(Tpair(address,log));
    }

    template <typename T>
    void Logger::record(const T & reference, const std::string& filename)
    {
      record<T>(&reference,filename);
    }

    template <typename T>
    void Logger::push(const T * address)
    {
      Tmap::iterator i= logs_.find( * address );

      if (i!=logs_.end())
        update_(i);
      else
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
    }

    template <typename T>
    void Logger::push(const T & reference)
    {
      push(*reference);
    }

    template <typename T>
    void Logger::updateAddress(const void * oldAddress, const T * newAddress)
    {
      Tmap::iterator i= logs_.find(oldAddress);
      if (i==logs_.end())
      {
        throw std::invalid_argument
                ("The logger cannot find the data, please use record function");
        BOOST_ASSERT(false && "Logger : the address to update is not found");
      }
      else
      {
        log_s log(i->second.filename);
        log.setType<T>();
        logs_.insert(Tpair(newAddress,log));
        logs_.erase(i);
      }
    }

  }
}
