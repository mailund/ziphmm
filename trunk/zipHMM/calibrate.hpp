#ifndef CALIBRATE_HPP
#define CALIBRATE_HPP

#include "performance_description.hpp"

#include <vector>
#include <string>


namespace zipHMM {

  void calibrateDevice(DeviceDescriptor &result, const DeviceTypeId &typeId);
  void calibrateDevices(std::vector<DeviceDescriptor> &result);
  void calibrate(const std::string &filename);

}

#endif
