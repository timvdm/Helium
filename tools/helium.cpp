#include "header.h"
#include "index.h"
#include "fold.h"
#include "similarity.h"

#include "../src/util.h"

using namespace Helium;

int main(int argc, char**argv)
{
  std::vector<HeliumTool*> tools;
  tools.push_back(new HeaderTool);
  tools.push_back(new IndexTool);
  tools.push_back(new FoldTool);
  tools.push_back(new SimilarityTool);

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <tool> [options]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Tools:" << std::endl;
    for (std::size_t i = 0; i < tools.size(); ++i)
      std::cerr << "    " << tools[i]->name() << "\t" << tools[i]->description() << std::endl;
    std::cerr << std::endl;
    return 0;
  }

  HeliumTool *tool = 0;
  for (std::size_t i = 0; i < tools.size(); ++i)
    if (tools[i]->name() == argv[1]) {
      tool = tools[i];
      break;
    }

  if (!tool) {
    std::cerr << argv[1] << " is not a recognised tool" << std::endl;
    return 0;
  }

  if (argc < 2 + tool->requiredArgs()) {
    std::cerr << tool->usage(make_string(argv[0], " ", argv[1]));
    return 0;
  }

  int ret = tool->run(argc, argv);

  for (std::size_t i = 0; i < tools.size(); ++i)
    delete tools[i];  

  return ret;
}
