/*
 * Copyright (c) 2013, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#include "tool.h"

#include <iostream>
#include <sstream>

#include <CL/cl.hpp>

namespace Helium {

  void checkError(cl_int ret, const std::string &msg)
  {
    if (ret != CL_SUCCESS) {
      std::string code;
      switch (ret) {
        case CL_DEVICE_NOT_FOUND:
          code = "CL_DEVICE_NOT_FOUND";
          break;
        case CL_DEVICE_NOT_AVAILABLE:
          code = "CL_DEVICE_NOT_AVAILABLE";
          break;
        case CL_COMPILER_NOT_AVAILABLE:
          code = "CL_COMPILER_NOT_AVAILABLE";
          break;
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
          code = "CL_MEM_OBJECT_ALLOCATION_FAILURE";
          break;
        case CL_OUT_OF_RESOURCES:
          code = "CL_OUT_OF_RESOURCES";
          break;
        case CL_OUT_OF_HOST_MEMORY:
          code = "CL_OUT_OF_HOST_MEMORY";
          break;
        case CL_PROFILING_INFO_NOT_AVAILABLE:
          code = "CL_PROFILING_INFO_NOT_AVAILABLE";
          break;
        case CL_MEM_COPY_OVERLAP:
          code = "CL_MEM_COPY_OVERLAP";
          break;
        case CL_IMAGE_FORMAT_MISMATCH:
          code = "CL_IMAGE_FORMAT_MISMATCH";
          break;
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
          code = "CL_IMAGE_FORMAT_NOT_SUPPORTED";
          break;
        case CL_BUILD_PROGRAM_FAILURE:
          code = "CL_BUILD_PROGRAM_FAILURE";
          break;
        case CL_MAP_FAILURE:
          code = "CL_MAP_FAILURE";
          break;
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
          code = "CL_MISALIGNED_SUB_BUFFER_OFFSET";
          break;
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
          code = "CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
          break;
        case CL_COMPILE_PROGRAM_FAILURE:
          code = "CL_COMPILER_PROGRAM_FAILURE";
          break;
        case CL_LINKER_NOT_AVAILABLE:
          code = "CL_LINKER_NOT_AVAILABLE";
          break;
        case CL_LINK_PROGRAM_FAILURE:
          code = "CL_LINK_PROGRAM_FAILURE";
          break;
        case CL_DEVICE_PARTITION_FAILED:
          code = "CL_DEVICE_PARTITION_FAILED";
          break;
        case CL_KERNEL_ARG_INFO_NOT_AVAILABLE:
          code = "CL_KERNEL_ARG_INFO_NOT_AVAILABLE";
          break;
        case CL_INVALID_VALUE:
          code = "CL_INVALID_VALUE";
          break;
        case CL_INVALID_DEVICE_TYPE:
          code = "CL_INVALID_DEVICE_TYPE";
          break;
        case CL_INVALID_PLATFORM:
          code = "CL_INVALID_PLATFORM";
          break;
        case CL_INVALID_DEVICE:
          code = "CL_INVALID_DEVICE";
          break;
        case CL_INVALID_CONTEXT:
          code = "CL_INVALID_CONTEXT";
          break;
        case CL_INVALID_QUEUE_PROPERTIES:
          code = "CL_INVALID_QUEUE_PROPERTIES";
          break;
        case CL_INVALID_COMMAND_QUEUE:
          code = "CL_INVALID_COMMAND_QUEUE";
          break;
        case CL_INVALID_HOST_PTR:
          code = "CL_INVALID_HOST_PTR";
          break;
        case CL_INVALID_MEM_OBJECT:
          code = "CL_INVALID_MEM_OBJECT";
          break;
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
          code = "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
          break;
        case CL_INVALID_IMAGE_SIZE:
          code = "CL_INVALID_IMAGE_SIZE";
          break;
        case CL_INVALID_SAMPLER:
          code = "CL_INVALID_SAMPLER";
          break;
        case CL_INVALID_BINARY:
          code = "CL_INVALID_BINARY";
          break;
        case CL_INVALID_BUILD_OPTIONS:
          code = "CL_INVALID_BUILD_OPTIONS";
          break;
        case CL_INVALID_PROGRAM:
          code = "CL_INVALID_PROGRAM";
          break;
        case CL_INVALID_PROGRAM_EXECUTABLE:
          code = "CL_INVALID_PROGRAM_EXECUTABLE";
          break;
        case CL_INVALID_KERNEL_NAME:
          code = "CL_INVALID_KERNEL_NAME";
          break;
        case CL_INVALID_KERNEL_DEFINITION:
          code = "CL_INVALID_KERNEL_DEFINITION";
          break;
        case CL_INVALID_KERNEL:
          code = "CL_INVALID_KERNEL";
          break;
        case CL_INVALID_ARG_INDEX:
          code = "CL_INVALID_ARG_INDEX";
          break;
        case CL_INVALID_ARG_VALUE:
          code = "CL_INVALID_ARG_VALUE";
          break;
        case CL_INVALID_ARG_SIZE:
          code = "CL_INVALID_ARG_SIZE";
          break;
        case CL_INVALID_KERNEL_ARGS:
          code = "CL_INVALID_KERNEL_ARGS";
          break;
        case CL_INVALID_WORK_DIMENSION:
          code = "CL_INVALID_WORK_DIMENSION";
          break;
        case CL_INVALID_WORK_GROUP_SIZE:
          code = "CL_INVALID_GROUP_SIZE";
          break;
        case CL_INVALID_WORK_ITEM_SIZE:
          code = "CL_INVALID_ITEM_SIZE";
          break;
        case CL_INVALID_GLOBAL_OFFSET:
          code = "CL_INVALID_GLOBAL_OFFSET";
          break;
        case CL_INVALID_EVENT_WAIT_LIST:
          code = "CL_INVALID_EVENT_WAIT_LIST";
          break;
        case CL_INVALID_EVENT:
          code = "CL_INVALID_EVENT";
          break;
        case CL_INVALID_OPERATION:
          code = "CL_INVALID_OPERATION";
          break;
        case CL_INVALID_GL_OBJECT:
          code = "CL_INVALID_GL_OBJECT";
          break;
        case CL_INVALID_BUFFER_SIZE:
          code = "CL_INVALID_BUFFER_SIZE";
          break;
        case CL_INVALID_MIP_LEVEL:
          code = "CL_INVALID_MIP_LEVEL";
          break;
        case CL_INVALID_GLOBAL_WORK_SIZE:
          code = "CL_INVALID_GLOBAL_WORK_SIZE";
          break;
        case CL_INVALID_PROPERTY:
          code = "CL_INVALID_PROPERTY";
          break;
        case CL_INVALID_IMAGE_DESCRIPTOR:
          code = "CL_INVALID_IMAGE_DESCRIPTO";
          break;
        case CL_INVALID_COMPILER_OPTIONS:
          code = "CL_INVALID_COMPILER_OPTIONS";
          break;
        case CL_INVALID_LINKER_OPTIONS:
          code = "CL_INVALID_LINKER_OPTIONS";
          break;
        case CL_INVALID_DEVICE_PARTITION_COUNT:
          code = "CL_INVALID_PARTITION_COUNT";
          break;
        default:
          code = "CL_???";
          break;
      }

      std::cerr << "Error (" << code << "): " << msg << std::endl;
      std::exit(-1);
    }
  }
  
  /**
   * Tool for listing OpenCL information
   */
  class OpenCLTool : public HeliumTool
  {
    public:
      /**
       * Perform tool action.
       */
      int run(int argc, char **argv)
      {
        // get a list of OpenCL platforms
        std::vector<cl::Platform> platforms;
        checkError(cl::Platform::get(&platforms), "Could not get platforms");

        if (platforms.empty()) {
          std::cout << "No OpenCL platforms found" << std::endl;
          return 0;
        }

        // foreach platform
        for (std::size_t i = 0; i < platforms.size(); ++i) {
          cl::Platform &platform = platforms[i];
          std::cout << "Platform " << i + 1 << ":" << std::endl;

          // get platform info
          std::string name, profile, vendor, version, extensions;
          checkError(platform.getInfo(CL_PLATFORM_NAME, &name), "Could not get platform name");
          checkError(platform.getInfo(CL_PLATFORM_PROFILE, &profile), "Could not get platform profile");
          checkError(platform.getInfo(CL_PLATFORM_VENDOR, &vendor), "Could not get platform vendor");
          checkError(platform.getInfo(CL_PLATFORM_VERSION, &version), "Could not get platform version");
          checkError(platform.getInfo(CL_PLATFORM_EXTENSIONS, &extensions), "Could not get platform extensions");

          std::cout << "    Name: " << name << std::endl;
          std::cout << "    Profile: " << profile << std::endl;
          std::cout << "    Vendor: " << vendor << std::endl;
          std::cout << "    Version: " << version << std::endl;
          std::cout << "    Extensions: " << extensions << std::endl;
          std::cout << std::endl;

          // get a list of devices for this platform
          std::vector<cl::Device> devices;
          checkError(platform.getDevices(CL_DEVICE_TYPE_ALL, &devices), "Could not get platform devices");

          // foreach device
          for (std::size_t j = 0; j < devices.size(); ++j) {
            cl::Device &device = devices[j];

            std::cout << "    Device " << j + 1 << ":" << std::endl;

            // get device info
            std::string name, profile, vendor, version, extensions, opencl_c_version, driver_version, built_in_kernels; 
            checkError(device.getInfo(CL_DEVICE_NAME, &name), "Could not get device name");
            checkError(device.getInfo(CL_DEVICE_PROFILE, &profile), "Could not get device profile");
            checkError(device.getInfo(CL_DEVICE_VENDOR, &vendor), "Could not get device vendor");
            checkError(device.getInfo(CL_DEVICE_VERSION, &version), "Could not get device version");
            checkError(device.getInfo(CL_DEVICE_EXTENSIONS, &extensions), "Could not get device extensions");
            checkError(device.getInfo(CL_DEVICE_OPENCL_C_VERSION, &opencl_c_version), "Could not get device opencl_c_version");
            checkError(device.getInfo(CL_DRIVER_VERSION, &driver_version), "Could not get device driver_version");
            //checkError(device.getInfo(CL_DEVICE_BUILT_IN_KERNELS, &built_in_kernels), "Could not get device built_in_kernels");
            std::vector<size_t> max_work_item_sizes;
            checkError(device.getInfo(CL_DEVICE_MAX_WORK_ITEM_SIZES, &max_work_item_sizes), "Could not get device max work item sizes");
            //std::vector<cl_device_partition_property> partition_properties, partition_type;
            //checkError(device.getInfo(CL_DEVICE_PARTITION_PROPERTIES, &partition_properties), "Could not get device partition properties");
            //checkError(device.getInfo(CL_DEVICE_PARTITION_TYPE, &partition_type), "Could not get device partition type");

            std::cout << "        Name: " << name << std::endl;
            std::cout << "        Profile: " << profile << std::endl;
            std::cout << "        Vendor: " << vendor << std::endl;
            std::cout << "        Version: " << version << std::endl;
            std::cout << "        Extensions: " << extensions << std::endl;
            std::cout << "        OpenCL C Version: " << opencl_c_version << std::endl;
            std::cout << "        Driver Version: " << driver_version << std::endl;
            std::cout << "        Max Work Item Sizes: ";
            for (std::size_t k = 0; k < max_work_item_sizes.size(); ++k)
              std::cout << max_work_item_sizes[k] << " ";
            std::cout << std::endl;

            std::cout << std::endl;
          }

        }


        return 0;
      }

  };

  class OpenCLToolFactory : public HeliumToolFactory
  {
    public:
      HELIUM_TOOL("opencl", "Get OpenCL information", 0, OpenCLTool);

      /**
       * Get usage information.
       */
      std::string usage(const std::string &command) const
      {
        std::stringstream ss;
        ss << "Usage: " << command << std::endl;
        ss << std::endl;
        return ss.str();
      }
  };

  OpenCLToolFactory theOpenCLToolFactory;

}
