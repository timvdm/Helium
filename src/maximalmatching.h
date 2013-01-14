#ifndef HELIUM_INVARIANTS_H
#define HELIUM_INVARIANTS_H

#include "tie.h"
#include "molecule.h"
#include "util.h"

#include <Eigen/Dense>

#include <vector>
#include <iostream>

namespace Helium {

  unsigned int uniform(unsigned int min, unsigned int max)
  {
    return (static_cast<double>(rand()) / (static_cast<unsigned long>(RAND_MAX) + 1)) * (max - min + 1) + min;
  }

  template<typename MoleculeType>
  std::vector<bool> maximal_matching(MoleculeType *mol)
  {
    srand(time(NULL));
    typedef typename molecule_traits<MoleculeType>::mol_bond_iter mol_bond_iter;

    unsigned int n = num_atoms(mol);
    // compute adjacency matrix D
    Eigen::MatrixXi D(n, n);
    D.fill(0);

    std::vector<unsigned int> w_ij;

    mol_bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol);
    for (; bond != end_bonds; ++bond) {
      unsigned int source = get_index(mol, get_source(mol, *bond));
      unsigned int target = get_index(mol, get_target(mol, *bond));
      D(source, target) = 1;
      D(target, source) = 1;
      // choose w_ij from [1, 2m] independantly and uniformly
      //w_ij.push_back(1 + rand() % (2 * num_bonds(mol)));
      //w_ij.push_back(2 * (1 + get_index(mol, *bond)));
      w_ij.push_back(uniform(1, 2 * num_bonds(mol)));
    }
    std::cout << "w_ij: " << w_ij << std::endl;
    std::cout << "D:" << std::endl << D << std::endl;

    unsigned int indeterminant = 1;

    // B: replace 1's with 2^w_ij 
    Eigen::MatrixXi &B(D);
    for (int i = 0; i < B.rows(); ++i)
      for (int j = i + 1; j < B.cols(); ++j) {
        if (!B(i, j))
          continue;
        unsigned int wij = w_ij[get_index(mol, get_bond(mol, get_atom(mol, i), get_atom(mol, j)))];
        unsigned int power = std::pow(2, wij);
        B(i, j) = power;
        B(j, i) = -(get_index(mol, get_bond(mol, get_atom(mol, i), get_atom(mol, j))) + 1);
      }
    std::cout << "B:" << std::endl << B << std::endl;
    
    int det = B.determinant();
    std::cout << "|B| = " << det << std::endl;

    std::vector<bool> result;
    return result;
  }


}

#endif
