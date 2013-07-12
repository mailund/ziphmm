#ifndef CALIBRATE_HPP
#define CALIBRATE_HPP

#include "performance_description.hpp"

#include <vector>



namespace zipHMM {

    void calibrateDevice(DeviceDescriptor &result, const DeviceTypeId &typeId);
    void calibrateDevices(std::vector<DeviceDescriptor> &result);

}

#endif
