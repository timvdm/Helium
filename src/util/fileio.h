#ifndef HELIUM_UTIL_FILEIO_H
#define HELIUM_UTIL_FILEIO_H

#include <istream>
#include <ostream>

namespace Helium {

  template<typename T>
  void write64(std::ostream &os, T value)
  {
    unsigned long copy = value;
    os.write(reinterpret_cast<const char*>(&copy), sizeof(unsigned long));
  }

  template<typename T>
  void write32(std::ostream &os, T value)
  {
    unsigned int copy = value;
    os.write(reinterpret_cast<const char*>(&copy), sizeof(unsigned int));
  }

  template<typename T>
  void write16(std::ostream &os, T value)
  {
    unsigned short copy = value;
    os.write(reinterpret_cast<const char*>(&copy), sizeof(unsigned short));
  }

  template<typename T>
  void write8(std::ostream &os, T value)
  {
    unsigned char copy = value;
    os.write(reinterpret_cast<const char*>(&copy), sizeof(unsigned char));
  }


  template<typename T>
  void read64(std::istream &is, T &value)
  {
    StaticAssert<sizeof(T) == sizeof(long)>();
    is.read(reinterpret_cast<char*>(&value), sizeof(unsigned long));
  }

  template<typename T>
  void read32(std::istream &is, T &value)
  {
    StaticAssert<sizeof(T) == sizeof(int)>();
    is.read(reinterpret_cast<char*>(&value), sizeof(unsigned int));
  }

  template<typename T>
  void read16(std::istream &is, T &value)
  {
    StaticAssert<sizeof(T) == sizeof(short)>();
    is.read(reinterpret_cast<char*>(&value), sizeof(unsigned short));
  }

  template<typename T>
  void read8(std::istream &is, T &value)
  {
    StaticAssert<sizeof(T) == sizeof(char)>();
    is.read(reinterpret_cast<char*>(&value), sizeof(unsigned char));
  }

}

#endif
