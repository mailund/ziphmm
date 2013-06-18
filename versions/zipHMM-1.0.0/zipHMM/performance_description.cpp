#include "performance_description.hpp"

#include "io_utils.hpp"

#include <fstream>

namespace zipHMM {
  
  const PThreadProcessingDeviceFactory PTHREAD_PROCESSING_DEVICE_FACTORY;
  const DeviceFactory *ALL_DEVICE_FACTORIES[] = { &PTHREAD_PROCESSING_DEVICE_FACTORY };
  const unsigned N_DEVICE_FACTORIES = 1;

  const std::string DEFAULT_DEVICE_FILENAME = get_user_dir() + "/.ziphmm.devices";

  void readDescriptors(std::vector<DeviceDescriptor> &resultDescriptors, const std::string &filename) {
    std::ifstream in(filename.c_str());
    
    if(!in) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }

    readDescriptors(resultDescriptors, in);

    in.close();
  }

  const DeviceFactory &deviceIdToFactory(const DeviceTypeId &typeId) {
    for(unsigned i = 0; i < N_DEVICE_FACTORIES; ++i) {
      if(ALL_DEVICE_FACTORIES[i]->getTypeId() == typeId)
	return *(ALL_DEVICE_FACTORIES[i]);
    }
    
    std::cerr << typeId << " is not a valid device name." << std::endl;
    exit(-1);
  }

  namespace {

    template<typename T>
    void assertEqual(const T& x1, const T& x2, const std::string &text = "Assertion failed:") {
      if(x1 != x2) {
	std::cerr << text << " Found " << x1 << " but expected " << x2 << std::endl;
	std::exit(-1);
      }
    }
    
  } // namespace


void readDescriptors(std::vector<DeviceDescriptor> &resultDescriptors, std::istream &in) {

    resultDescriptors.clear();

    read_token_or_exit(in, "nDeviceKinds");
    const unsigned nDeviceKinds = read_or_exit<unsigned>(in, "number of kinds of devices");

    for(unsigned iDevice = 0; iDevice < nDeviceKinds; ++iDevice) {
      read_token_or_exit(in, "device");
      const std::string deviceName = read_or_exit<std::string>(in, "device name");
      
      read_token_or_exit(in, "algorithm");
      const std::string algorithmName = read_or_exit<std::string>(in, "algorithm name");
      assertEqual(std::string("likelihood"), algorithmName, "is not a valid algorithm.");
	
      read_token_or_exit(in, "nDevices");
      size_t nDevices = read_or_exit<unsigned>(in, "number of devices");
    
      resultDescriptors.push_back(DeviceDescriptor(deviceIdToFactory(deviceName), nDevices));
    }
  }
  
  void writeDescriptors(const std::vector<DeviceDescriptor> &descriptors, const std::string filename) {
    std::ofstream out(filename.c_str());
    
    if(!out) {
      std::cerr << "Unable to open \"" << filename << "\"" << std::endl;
      exit(-1);
    }
    
    writeDescriptors(descriptors, out);

    out.close();
  }

  void writeDescriptors(const std::vector<DeviceDescriptor> &descriptors, std::ostream &out) {
    out << "nDeviceKinds " << descriptors.size() << std::endl;
    for(unsigned iDevice = 0; iDevice < descriptors.size(); ++iDevice) {
      const DeviceDescriptor &desc = descriptors[iDevice];
      out << "  device " << desc.getDeviceName() << std::endl;
      
      out << "    algorithm " << "likelihood" << std::endl;
      out << "      nDevices " << desc.getNDevices() << std::endl;
    }
  }
  
} // namespace zipHMM
