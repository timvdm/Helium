#include <Helium/stereo.h>

#include "test.h"

using namespace Helium;

void test_tetrahedral_class()
{
  Stereo::Ref refs[4] = {0, 3, 2, 1};
  auto storage = StereoStorage::create(Stereo::Tetrahedral, 4, refs, refs + 4);

  COMPARE(Stereo::TH1, tetrahedral_class(storage, 0, 3, 2, 1));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 0, 2, 1, 3));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 0, 1, 3, 2));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 0, 1, 2, 3));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 0, 2, 3, 1));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 0, 3, 1, 2));

  COMPARE(Stereo::TH1, tetrahedral_class(storage, 1, 2, 3, 0));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 1, 3, 0, 2));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 1, 0, 2, 3));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 1, 0, 3, 2));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 1, 3, 2, 0));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 1, 2, 0, 3));

  COMPARE(Stereo::TH1, tetrahedral_class(storage, 2, 1, 0, 3));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 2, 3, 1, 0));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 2, 0, 3, 1));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 2, 3, 0, 1));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 2, 0, 1, 3));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 2, 1, 3, 0));

  COMPARE(Stereo::TH1, tetrahedral_class(storage, 3, 1, 2, 0));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 3, 2, 0, 1));
  COMPARE(Stereo::TH1, tetrahedral_class(storage, 3, 0, 1, 2));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 3, 0, 2, 1));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 3, 2, 1, 0));
  COMPARE(Stereo::TH2, tetrahedral_class(storage, 3, 1, 0, 2));
}

int main()
{
  test_tetrahedral_class();
}
