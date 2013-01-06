#ifndef HELIUM_UTIL_STRING_H
#define HELIUM_UTIL_STRING_H

#include <string>
#include <sstream>
#include <vector>

namespace Helium {

  /**
   * Strings
   */
  template<typename T1>
  inline std::string make_string(const T1 &t1)
  {
    std::stringstream ss;
    ss << t1;
    return ss.str();
  }

  template<typename T1, typename T2>
  inline std::string make_string(const T1 &t1, const T2 &t2)
  {
    std::stringstream ss;
    ss << t1 << t2;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9;
    return ss.str();
  }

  template<typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  inline std::string make_string(const T1 &t1, const T2 &t2, const T3 &t3, const T4 &t4, const T5 &t5, const T6 &t6, const T7 &t7, const T8 &t8, const T9 &t9, const T10 &t10)
  {
    std::stringstream ss;
    ss << t1 << t2 << t3 << t4 << t5 << t6 << t7 << t8 << t9 << t10;
    return ss.str();
  }
  
  inline void replace_first(std::string &str, const std::string &what, const std::string &with = std::string())
  {
    std::size_t pos = str.find(what);
    if (pos != std::string::npos)
      str.replace(pos, what.size(), with);
  }

  inline std::string replace_first(const std::string &str, const std::string &what, const std::string &with = std::string())
  {
    std::string tmp(str);
    replace_first(tmp, what, with);
    return tmp;
  }

  inline void replace_all(std::string &str, const std::string &what, const std::string &with = std::string())
  {
    std::size_t last_pos = 0, pos;
    while ((pos = str.find(what, last_pos)) != std::string::npos) {
      str.replace(pos, what.size(), with);
      last_pos = pos + 1;
    }
  }

  inline std::string replace_all(const std::string &str, const std::string &what, const std::string &with = std::string())
  {
    std::string tmp(str);
    replace_all(tmp, what, with);
    return tmp;
  }

  inline void strip(std::string &str)
  {
    while (str.size() && str[0] == ' ')
      str = str.substr(1);
    while (str.size() && str[str.size() - 1] == ' ')
      str.resize(str.size() - 1);
  }

  inline std::string stripped(const std::string &str)
  {
    std::string result(str);
    strip(result);
    return result;
  }

  template<typename T>
  inline T string2number(const std::string &str)
  {
    std::stringstream ss(str);
    T number;
    ss >> number;
    return number;
  }

  /**
   * @param repeat The delimiter may be repeated. Only works for single
   *        character delimiter.
   */
  inline std::vector<std::string> tokenize(const std::string &str, const std::string &delimiter, bool repeat = false)
  {
    std::vector<std::string> tokens;
    std::size_t currpos = 0, nextpos = 0;
    //std::cout << "tokenize: \"" << str << "\"" << std::endl;

    while ((nextpos = str.find(delimiter, currpos)) != std::string::npos) {
      if (repeat)
        while (nextpos < str.size() && str[nextpos] == delimiter[0])
          ++nextpos;
      if (nextpos == str.size())
        return tokens;
      tokens.push_back(str.substr(currpos, nextpos - currpos - 1));
      //std::cout << "token: \"" << tokens.back() << "\"" << std::endl;
      currpos = nextpos;
    }
    tokens.push_back(str.substr(currpos, str.length() - currpos));
    //std::cout << "token: \"" << tokens.back() << "\"" << std::endl;

    return tokens;
  }

}

#endif
