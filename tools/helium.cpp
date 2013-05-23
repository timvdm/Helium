#include <Helium/util.h>
#include "tool.h"

using namespace Helium;

int main(int argc, char**argv)
{
  std::vector<HeliumToolFactory*> factories = HeliumTool::toolFactories();

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <tool> [options]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Tools:" << std::endl;
    for (std::size_t i = 0; i < factories.size(); ++i)
      std::cerr << "    " << factories[i]->name() << "\t" << factories[i]->description() << std::endl;
    std::cerr << std::endl;
    return 0;
  }

  HeliumToolFactory *factory = 0;
  for (std::size_t i = 0; i < factories.size(); ++i)
    if (factories[i]->name() == argv[1]) {
      factory = factories[i];
      break;
    }

  if (!factory) {
    std::cerr << argv[1] << " is not a recognised tool" << std::endl;
    return 0;
  }

  HeliumTool *tool = factory->instance();

  if (argc < 2 + factory->requiredArgs()) {
    std::cerr << factory->usage(make_string(argv[0], " ", argv[1]));
    return 0;
  }

  int ret = tool->run(argc, argv);

  delete tool;

  return ret;
}
