#include "performance_description.hpp"
#include "test_utils.hpp"

using namespace zipHMM;

int main(int argc, char ** argv) {
    const std::string inFilename = "test_data/test1.devices";
    std::vector<DeviceDescriptor> devices;
    readDescriptors(devices, inFilename);

    const std::string outFilename = "test_data/postTest1.devices";
    writeDescriptors(devices, outFilename);

    assertFilesEqual(inFilename, outFilename);
}
