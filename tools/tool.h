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
#ifndef HELIUM_HELIUMTOOL_H
#define HELIUM_HELIUMTOOL_H

#include <Helium/config.h>

#include <string>
#include <vector>

namespace Helium {

  class HeliumTool;

  class HeliumToolFactory
  {
    public:
      /**
       * Get tool name (e.g. "index", "query", ...)
       */
      virtual std::string name() const = 0;
      /**
       * Get description.
       */
      virtual std::string description() const = 0;
      /**
       * Get usage information.
       *
       * @param command The command (e.g. "helium index", "helium query", ...)
       */
      virtual std::string usage(const std::string &command) const = 0;
      /**
       * Get the number of required arguments.
       */
      virtual int requiredArgs() const = 0;
      /**
       * Get an instance of the tool.
       */
      virtual HeliumTool* instance() const = 0;
  };

#define HELIUM_TOOL(tool, desc, args, klass) \
  klass##Factory() { HeliumTool::registerTool(this); } \
  std::string name() const { return tool; } \
  std::string description() const { return desc; } \
  int requiredArgs() const { return args; } \
  HeliumTool* instance() const { return new klass; }

  /**
   * Base class for command line tools.
   */
  class HeliumTool
  {
    public:
      /**
       * Constructor.
       */
      virtual ~HeliumTool()
      {
      }

      /**
       * Get a list of all plugins.
       */
      static std::vector<HeliumToolFactory*>& toolFactories()
      {
        static std::vector<HeliumToolFactory*> *list = 0;
        if (!list)
          list = new std::vector<HeliumToolFactory*>();
        return *list;
      }

      /**
       * Register a tool facotry.
       */
      static void registerTool(HeliumToolFactory *factory)
      {
        toolFactories().push_back(factory);
      }

      /**
       * Perform tool action.
       */
      virtual int run(int argc, char **argv) = 0;
  };

}

#endif
