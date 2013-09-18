#ifndef PERFORMANCE_DESCRIPTION_HPP
#define PERFORMANCE_DESCRIPTION_HPP

#include "ProcessingDevice.hpp"
#include "PThreadProcessingDevice.hpp"

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace zipHMM {
  typedef std::string DeviceTypeId;

  /////////////////////////// DeviceFactory ///////////////////////////

  class DeviceFactory {
    DeviceTypeId typeId;
    bool multipleDevices;

  public:
    DeviceFactory(const DeviceTypeId &typeId, bool multipleDevices)
      : typeId(typeId), multipleDevices(multipleDevices)
    {}

    virtual ~DeviceFactory()
    {}

    virtual ProcessingDevice *createDevice(unsigned id) const = 0;

    const DeviceTypeId &getTypeId() const { return typeId; }
    bool supportsMultipleDevices() const { return multipleDevices; }
  }; // end class


  /////////////////////////// PThreadProcessingDeviceFactory ///////////////////////////

  struct PThreadProcessingDeviceFactory : public DeviceFactory {
    
    PThreadProcessingDeviceFactory()
      : DeviceFactory("PThreadProcessingDevice", true)
    {}
    
    virtual PThreadProcessingDevice *createDevice(unsigned id) const {
      return new PThreadProcessingDevice(id);
    }
  }; // end class
  
  
  extern const PThreadProcessingDeviceFactory PTHREAD_PROCESSING_DEVICE_FACTORY;
  extern const DeviceFactory *ALL_DEVICE_FACTORIES[];
  extern const unsigned N_DEVICE_FACTORIES;

  const DeviceFactory &deviceIdToFactory(const DeviceTypeId &typeId);

  /////////////////////////// DeviceDescriptor ///////////////////////////

  class DeviceDescriptor {
    const DeviceFactory *factory;
    size_t nDevices;

  public:
    DeviceDescriptor()
      : factory(0), nDevices()
    {}

    DeviceDescriptor(const DeviceFactory &factory, size_t nDevices)
      : factory(&factory), nDevices(nDevices)
    {}

    const DeviceTypeId getTypeId() const { return factory->getTypeId(); }

    const std::string &getDeviceName() const { return factory->getTypeId(); }

    size_t getNDevices() const {
      return nDevices;
    }

    const DeviceFactory &getFactory() const { return *factory; }
        
    ProcessingDevice *createDevice(unsigned i) const { return factory->createDevice(i); }
  }; // end class


  extern const std::string DEFAULT_DEVICE_FILENAME;

  void readDescriptors(std::vector<DeviceDescriptor> &resultDescriptors, const std::string &filename);

  void readDescriptors(std::vector<DeviceDescriptor> &resultDescriptors, std::istream &in);

  void writeDescriptors(const std::vector<DeviceDescriptor> &descriptors, const std::string filename);

  void writeDescriptors(const std::vector<DeviceDescriptor> &descriptors, std::ostream &out);
}

#endif
