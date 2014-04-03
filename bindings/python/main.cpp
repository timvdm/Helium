#include <boost/python.hpp>

using namespace boost::python;

void export_bitvec();
void export_common();
void export_element();
void export_error();
void export_molecule();
void export_smarts();
void export_smiles();
void export_smirks();
void export_substructure();
void export_rings();

void export_canonical();
void export_components();
void export_dfs();
void export_enumerate_paths();
void export_enumerate_subgraphs();
void export_extended_connectivities();
void export_invariants();

void export_file();
void export_moleculefile();
void export_fps();
void export_fingerprintfiles();

void export_fingerprints();
void export_similarity();

BOOST_PYTHON_MODULE(helium) {

  export_bitvec();
  export_common();
  export_element();
  export_error();
  export_molecule();
  export_smiles();
  export_smirks();
  export_smarts();
  //export_substructure();
  export_rings();

  // algorithms
  export_canonical();
  export_components();
  export_dfs();
  export_enumerate_paths();
  export_enumerate_subgraphs();
  export_extended_connectivities();
  export_invariants();

  //  fileio
  export_file();
  export_moleculefile();
  export_fps();
  export_fingerprintfiles();

  // fingerprints
  export_fingerprints();
  export_similarity();

}
