#ifndef HELIUM_KEKULIZE_H
#define HELIUM_KEKULIZE_H

#include <Helium/molecule.h>

#include <iostream>

#define KEKULIZE_DEBUG 1

namespace Helium {

  namespace impl {






  }

  template<typename MoleculeType>
  void kekulize(const MoleculeType &mol)
  {
    std::vector<bool> assigned(num_atoms(mol));




  }

}

#endif
