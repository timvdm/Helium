/**
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
#include "../src/util.h"

#include <map>

/**
 * Utility class to help processing command line arguments.
 */
class ParseArgs
{
  /**
   * Helper structure to hold information about optional arguments.
   */
  struct OptionalArg
  {
    OptionalArg(const std::string &option_) : option(option_)
    {
    }

    OptionalArg(const std::string &option_, const std::vector<std::string> &args_)
        : option(option_), args(args_)
    {
    }

    std::string option;
    std::vector<std::string> args;
  };

  public:
    /**
     * Helper function to create a list of arguments.
     */
    static std::vector<std::string> Args(
        const std::string &str1 = std::string(), const std::string &str2 = std::string(),
        const std::string &str3 = std::string(), const std::string &str4 = std::string(),
        const std::string &str5 = std::string(), const std::string &str6 = std::string(),
        const std::string &str7 = std::string(), const std::string &str8 = std::string(),
        const std::string &str9 = std::string(), const std::string &str10 = std::string(),
        const std::string &str11 = std::string(), const std::string &str12 = std::string(),
        const std::string &str13 = std::string(), const std::string &str14 = std::string())
    {
      std::vector<std::string> result;
      if (!str1.empty()) result.push_back(str1);
      if (!str2.empty()) result.push_back(str2);
      if (!str3.empty()) result.push_back(str3);
      if (!str4.empty()) result.push_back(str4);
      if (!str5.empty()) result.push_back(str5);
      if (!str6.empty()) result.push_back(str6);
      if (!str7.empty()) result.push_back(str7);
      if (!str8.empty()) result.push_back(str8);
      if (!str9.empty()) result.push_back(str9);
      if (!str10.empty()) result.push_back(str10);
      if (!str11.empty()) result.push_back(str11);
      if (!str12.empty()) result.push_back(str12);
      if (!str13.empty()) result.push_back(str13);
      if (!str14.empty()) result.push_back(str14);
      return result;
    }

    /**
     * Constructor. Use the Args() helper function to specify the @p optional
     * and @p required arguments. The optional arguments may also have parameters.
     * For example:
     *
     * @code
     * ParseArgs::Args("-foo(bar1,bar2)")
     * @endcode
     *
     * This is an optional argument -foo that has two arguments. Could be used
     * like this:
     *
     * @code
     * ./tool -foo bar 42
     * @endcode
     *
     * The actual names of the optional parameters are not used and are only there
     * for human readablility.
     *
     * @param argc The argument count (from main).
     * @param argv The argument array (from main).
     * @param optional Optional arguments including parameters.
     * @param required Required arguments.
     */
    ParseArgs(int argc, char **argv, const std::vector<std::string> &optional, const std::vector<std::string> &required) : m_valid(false)
    {
      if (static_cast<std::size_t>(argc) < required.size())
        return;

      if (static_cast<std::size_t>(argc) < required.size() + 1) {
        std::cerr << argv[0] << " requires at least " << required.size() << " argument(s)." << std::endl;
        m_valid = false;
        return;
      }

      std::vector<OptionalArg> optionalArgs;
      for (std::size_t i = 0; i < optional.size(); ++i) {
        std::size_t open = optional[i].find("(");
        std::size_t close = optional[i].find(")");
        if (open != std::string::npos && close != std::string::npos)
          optionalArgs.push_back(OptionalArg(optional[i].substr(0, open), Helium::tokenize(optional[i].substr(open + 1, close - open - 1), ",")));
        else
          optionalArgs.push_back(optional[i]);
      }

      for (std::size_t i = 1; i < argc - required.size(); ++i) {
        for (std::size_t j = 0; j < optionalArgs.size(); ++j) {
          std::string opt = optionalArgs[j].option;
          bool val = false;
          if (opt.rfind("$") != std::string::npos) {
            opt = opt.substr(0, opt.size() - 1);
            val = true;
          }

          if (std::string(argv[i]).substr(0, opt.size()) == opt) {
            if (val)
              m_args[opt] = std::string(argv[i]).substr(opt.size());
            else if (optionalArgs[j].args.size()) {
              if (i + optionalArgs[j].args.size() == static_cast<std::size_t>(argc)) {
                std::cerr << "Option " << opt << " requires " << optionalArgs[j].args.size() << "argument(s)." << std::endl;
                m_valid = false;
                return;
              }
              for (std::size_t k = 0; k < optionalArgs[j].args.size(); ++k)
                m_argargs[opt].push_back(argv[i + k + 1]);
              i += optionalArgs[j].args.size();
            } else
              m_args[opt] = "true";
          }
        }
      }

      int req = argc - required.size();
      for (std::size_t j = 0; j < required.size(); ++j)
        m_args[required[j]] = argv[req + j];

      m_valid = true;
    }

    /**
     * Were the arguments valid (i.e. enough arguments, recognised, ...)
     */
    bool IsValid() const
    {
      return m_valid;
    }

    /**
     * Get the number of arguments.
     */
    int GetArgCount() const
    {
      return m_args.size() + m_argargs.size();
    }

    /**
     * Check if the @p key argument was set (e.g. "-foo").
     */
    bool IsArg(const std::string &key)
    {
      if (m_args.find(key) != m_args.end())
        return true;
      return m_argargs.find(key) != m_argargs.end();
    }

    std::string GetArgString(const std::string &key)
    {
      if (m_args.find(key) != m_args.end())
        return m_args[key];
      return "";
    }

    int GetArgInt(const std::string &key)
    {
      if (m_args.find(key) == m_args.end())
        return -1;
      std::stringstream ss(m_args[key]);
      int result;
      ss >> result;
      return result;
    }

    double GetArgDouble(const std::string &key)
    {
      if (m_args.find(key) == m_args.end())
        return -1;
      std::stringstream ss(m_args[key]);
      double result;
      ss >> result;
      return result;
    }

    std::string GetArgString(const std::string &key, int arg)
    {
      if (m_argargs.find(key) != m_argargs.end())
        return m_argargs[key][arg];
      return "";
    }

    int GetArgInt(const std::string &key, int arg)
    {
      if (m_argargs.find(key) == m_argargs.end())
        return -1;
      std::stringstream ss(m_argargs[key][arg]);
      int result;
      ss >> result;
      return result;
    }

    double GetArgDouble(const std::string &key, int arg)
    {
      if (m_argargs.find(key) == m_argargs.end())
        return -1;
      std::stringstream ss(m_argargs[key][arg]);
      double result;
      ss >> result;
      return result;
    }

  private:
    std::map<std::string, std::string> m_args;
    std::map<std::string, std::vector<std::string> > m_argargs;
    bool m_valid;
};
