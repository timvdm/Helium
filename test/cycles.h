#include <Helium/smiles.h>
#include <Helium/fileio/molecules.h>

#include "test.h"

#include <numeric>

using namespace Helium;

class TestCycle
{
  public:
    TestCycle(Size numEdges = 0,
        unsigned int e1 = -1,
        unsigned int e2 = -1,
        unsigned int e3 = -1,
        unsigned int e4 = -1,
        unsigned int e5 = -1,
        unsigned int e6 = -1,
        unsigned int e7 = -1,
        unsigned int e8 = -1,
        unsigned int e9 = -1,
        unsigned int e10 = -1,
        unsigned int e11 = -1,
        unsigned int e12 = -1,
        unsigned int e13 = -1,
        unsigned int e14 = -1,
        unsigned int e15 = -1,
        unsigned int e16 = -1,
        unsigned int e17 = -1,
        unsigned int e18 = -1,
        unsigned int e19 = -1,
        unsigned int e20 = -1)
    {
      m_cycle.resize(numEdges);
      if (e1 != -1) m_cycle[e1] = true;
      if (e2 != -1) m_cycle[e2] = true;
      if (e3 != -1) m_cycle[e3] = true;
      if (e4 != -1) m_cycle[e4] = true;
      if (e5 != -1) m_cycle[e5] = true;
      if (e6 != -1) m_cycle[e6] = true;
      if (e7 != -1) m_cycle[e7] = true;
      if (e8 != -1) m_cycle[e8] = true;
      if (e9 != -1) m_cycle[e9] = true;
      if (e10 != -1) m_cycle[e10] = true;
      if (e11 != -1) m_cycle[e11] = true;
      if (e12 != -1) m_cycle[e12] = true;
      if (e13 != -1) m_cycle[e13] = true;
      if (e14 != -1) m_cycle[e14] = true;
      if (e15 != -1) m_cycle[e15] = true;
      if (e16 != -1) m_cycle[e16] = true;
      if (e17 != -1) m_cycle[e17] = true;
      if (e18 != -1) m_cycle[e18] = true;
      if (e19 != -1) m_cycle[e19] = true;
      if (e20 != -1) m_cycle[e20] = true;
    }

    std::size_t size() const
    {
      return std::accumulate(m_cycle.begin(), m_cycle.end(), 0);
    }

    const std::vector<bool>& cycle() const
    {
      return m_cycle;
    }

    std::vector<bool>& cycle()
    {
      return m_cycle;
    }

    bool operator==(const TestCycle &other) const
    {
      return m_cycle == other.m_cycle;
    }

    static std::vector<TestCycle> list(const TestCycle &cycle1 = TestCycle(),
                                       const TestCycle &cycle2 = TestCycle(),
                                       const TestCycle &cycle3 = TestCycle(),
                                       const TestCycle &cycle4 = TestCycle(),
                                       const TestCycle &cycle5 = TestCycle(),
                                       const TestCycle &cycle6 = TestCycle(),
                                       const TestCycle &cycle7 = TestCycle(),
                                       const TestCycle &cycle8 = TestCycle(),
                                       const TestCycle &cycle9 = TestCycle(),
                                       const TestCycle &cycle10 = TestCycle(),
                                       const TestCycle &cycle11 = TestCycle(),
                                       const TestCycle &cycle12 = TestCycle(),
                                       const TestCycle &cycle13 = TestCycle(),
                                       const TestCycle &cycle14 = TestCycle(),
                                       const TestCycle &cycle15 = TestCycle())
    {
      std::vector<TestCycle> cycles;
      if (cycle1.size()) cycles.push_back(cycle1);
      if (cycle2.size()) cycles.push_back(cycle2);
      if (cycle3.size()) cycles.push_back(cycle3);
      if (cycle4.size()) cycles.push_back(cycle4);
      if (cycle5.size()) cycles.push_back(cycle5);
      if (cycle6.size()) cycles.push_back(cycle6);
      if (cycle7.size()) cycles.push_back(cycle7);
      if (cycle8.size()) cycles.push_back(cycle8);
      if (cycle9.size()) cycles.push_back(cycle9);
      if (cycle10.size()) cycles.push_back(cycle10);
      if (cycle11.size()) cycles.push_back(cycle11);
      if (cycle12.size()) cycles.push_back(cycle12);
      if (cycle13.size()) cycles.push_back(cycle13);
      if (cycle14.size()) cycles.push_back(cycle14);
      if (cycle15.size()) cycles.push_back(cycle15);
      return cycles;
    }

  private:
    std::vector<bool> m_cycle; // incident vector (i.e. bond indices)
};


template<typename CyclePerceptionAlgorithm>
void test_cycle_perception(const CyclePerceptionAlgorithm &algorithm, const std::string &name,
    const HeMol &mol, const std::vector<TestCycle> &correct)
{
  std::cout << "Testing: " << name << std::endl;

  std::vector<TestCycle> cycles = algorithm(mol);

  for (std::size_t i = 0; i < correct.size(); ++i)
    ASSERT(correct[i].cycle().size() == num_bonds(mol));
  for (std::size_t i = 0; i < cycles.size(); ++i)
    ASSERT(cycles[i].cycle().size() == num_bonds(mol));

  COMPARE(correct.size(), cycles.size());

  for (std::size_t i = 0; i < correct.size(); ++i) {
    bool found = false;
    for (std::size_t j = 0; j < cycles.size(); ++j) {
      if (correct[i] == cycles[j]) {
        found = true;
        break;
      }
    }

    ASSERT(found);
  }
}

template<typename CyclePerceptionAlgorithm>
void test_cycle_perception(const CyclePerceptionAlgorithm &algorithm, const std::string &name,
    const std::string &smiles, const std::vector<TestCycle> &correct)
{
  HeMol mol;
  try {
    parse_smiles(smiles, mol);
  } catch (const Smiley::Exception &e) {
    std::cout << e.what() << std::endl;
    return;
  }

  test_cycle_perception(algorithm, name, mol, correct);
}

template<typename CyclePerceptionAlgorithm>
void test_cycle_perception(const CyclePerceptionAlgorithm &algorithm)
{
  //
  // Structures from Figueras paper page 989
  //
  test_cycle_perception(algorithm, "Figueras 1", "C1CCCCC1",
      TestCycle::list(TestCycle(6, 0, 1, 2, 3, 4, 5)));
  test_cycle_perception(algorithm, "Figueras 2", "C1CC2CCC3C4CCCC4CCC3C2CC1",
      TestCycle::list(TestCycle(20, 0, 1, 16, 17, 18, 19),
                      TestCycle(20, 2, 3, 4, 14, 15, 16),
                      TestCycle(20, 5, 10, 11, 12, 13, 14),
                      TestCycle(20, 6, 7, 8, 9, 10)));
  test_cycle_perception(algorithm, "Figueras 3", "C1CC2CCC3CCC4CCC5CCC6CCC1C7C6C5C4C3C27",
      TestCycle::list(TestCycle(30, 0, 1, 17, 18, 28, 29),
                      TestCycle(30, 2, 3, 4, 26, 27, 28),
                      TestCycle(30, 5, 6, 7, 24, 25, 26),
                      TestCycle(30, 8, 9, 10, 22, 23, 24),
                      TestCycle(30, 11, 12, 13, 20, 21, 22),
                      TestCycle(30, 14, 15, 16, 18, 19, 20),
                      TestCycle(30, 19, 21, 23, 25, 27, 29)));
  test_cycle_perception(algorithm, "Figueras 4", "C12CCC(C2)CC1",
      TestCycle::list(TestCycle(8, 0, 1, 2, 3, 4),
                      TestCycle(8, 3, 4, 5, 6, 7)));
  test_cycle_perception(algorithm, "Figueras 5", "C123CCC(CC2)(CC3)CC1",
      TestCycle::list(TestCycle(12, 0, 1, 2, 3, 4, 5),
                      TestCycle(12, 0, 1, 2, 6, 7, 8),
                      TestCycle(12, 0, 1, 2, 9, 10, 11),
                      TestCycle(12, 3, 4, 5, 6, 7, 8),
                      TestCycle(12, 3, 4, 5, 9, 10, 11),
                      TestCycle(12, 6, 7, 8, 9, 10, 11)));
  test_cycle_perception(algorithm, "Figueras 6", "C1CCCCC1CCC2CCCCC2",
      TestCycle::list(TestCycle(15, 0, 1, 2, 3, 4, 5),
                      TestCycle(15, 9, 10, 11, 12, 13, 14)));
  test_cycle_perception(algorithm, "Figueras 7", "C1C(C2)C3C(C4)C5C(CC6)CCC6C5C4C3C2C1",
      TestCycle::list(TestCycle(23, 0, 1, 20, 21, 22),
                      TestCycle(23, 1, 2, 18, 19, 20),
                      TestCycle(23, 3, 4, 16, 17, 18),
                      TestCycle(23, 4, 5, 14, 15, 16),
                      TestCycle(23, 6, 7, 8, 12, 13, 14),
                      TestCycle(23, 6, 9, 10, 11, 13, 14),
                      TestCycle(23, 7, 8, 9, 10, 11, 12)));
  test_cycle_perception(algorithm, "Figueras 8", "C12CCC(C1)CCC2",
      TestCycle::list(TestCycle(9, 0, 1, 2, 3, 4),
                      TestCycle(9, 3, 4, 5, 6, 7, 8)));
  test_cycle_perception(algorithm, "Figueras 9", "C12CCC(CC1)CCCC2",
      TestCycle::list(TestCycle(11, 0, 1, 2, 3, 4, 5),
                      TestCycle(11, 3, 4, 5, 6, 7, 8, 9, 10),
                      TestCycle(11, 0, 1, 2, 6, 7, 8, 9, 10)));
  test_cycle_perception(algorithm, "Figueras 10", "C12CCC(C3CC1)C3CC2",
      TestCycle::list(TestCycle(12, 0, 1, 2, 7, 9, 10, 11),
                      TestCycle(12, 3, 7, 8),
                      TestCycle(12, 4, 5, 6, 8, 9, 10, 11),
                      TestCycle(12, 0, 1, 2, 3, 4, 5, 6)));
  test_cycle_perception(algorithm, "Figueras 11", "C1C2C3CCCC3C(C2)CC1",
      TestCycle::list(TestCycle(13, 0, 8, 9, 10, 11, 12),
                      TestCycle(13, 1, 6, 7, 8, 9),
                      TestCycle(13, 2, 3, 4, 5, 6)));
  test_cycle_perception(algorithm, "Figueras 12", "C1CCC2C3CCCC(CC4)C3CC2C41",
      TestCycle::list(TestCycle(18, 0, 1, 2, 14, 15, 17),
                      TestCycle(18, 3, 11, 12, 13, 14),
                      TestCycle(18, 4, 5, 6, 7, 10, 11),
                      TestCycle(18, 8, 9, 10, 12, 13, 15, 16)));
  test_cycle_perception(algorithm, "Figueras 13", "C1C23CC34CC421",
      TestCycle::list(TestCycle(9, 0, 7, 8),
                      TestCycle(9, 1, 2, 3),
                      TestCycle(9, 3, 6, 7),
                      TestCycle(9, 4, 5, 6)));
  test_cycle_perception(algorithm, "Figueras 14", "C17C2C3C4C5C1C6C5C4C3C2C67",
      TestCycle::list(TestCycle(18, 0, 1, 2, 3, 4, 5),
                      TestCycle(18, 7, 9, 11, 13, 15, 16),
                      TestCycle(18, 0, 14, 15, 17),
                      TestCycle(18, 1, 12, 13, 14),
                      TestCycle(18, 2, 10, 11, 12),
                      TestCycle(18, 3, 8, 9, 10),
                      TestCycle(18, 4, 6, 7, 8),
                      TestCycle(18, 5, 6, 16, 17)));
  test_cycle_perception(algorithm, "Figueras 15", "C12C3C4C1C5C4C3C25",
      TestCycle::list(TestCycle(12, 0, 8, 9, 10),
                      TestCycle(12, 1, 6, 7, 8),
                      TestCycle(12, 2, 4, 5, 6),
                      TestCycle(12, 3, 4, 10, 11),
                      TestCycle(12, 0, 1, 2, 3),
                      TestCycle(12, 5, 7, 9, 11)));
  test_cycle_perception(algorithm, "Figueras 16", "C12C3C4C5C6C1C27C65C437",
      TestCycle::list(TestCycle(15, 0, 7, 13, 14),
                      TestCycle(15, 1, 12, 13),
                      TestCycle(15, 2, 10, 11, 12),
                      TestCycle(15, 3, 9, 10),
                      TestCycle(15, 4, 6, 8, 9),
                      TestCycle(15, 5, 6, 7),
                      TestCycle(15, 8, 11, 14)));
  test_cycle_perception(algorithm, "Figueras 17", "C12C3C4C5C6C1C2C(C65)C43",
      TestCycle::list(TestCycle(15, 0, 7, 8, 12, 14),
                      TestCycle(15, 2, 9, 11, 12, 13),
                      TestCycle(15, 4, 6, 8, 9, 10),
                      TestCycle(15, 1, 13, 14),
                      TestCycle(15, 3, 10, 11),
                      TestCycle(15, 5, 6, 7)));
  test_cycle_perception(algorithm, "Figueras 18", "C12C3C4C1C56C47CCCC37C25CCC6",
      TestCycle::list(TestCycle(20, 0, 11, 13, 14),
                      TestCycle(20, 1, 6, 11, 12),
                      TestCycle(20, 2, 4, 5, 6),
                      TestCycle(20, 3, 4, 14, 15),
                      TestCycle(20, 0, 1, 2, 3),
                      TestCycle(20, 5, 12, 13, 15),
                      TestCycle(20, 7, 8, 9, 10, 12),
                      TestCycle(20, 15, 16, 17, 18, 19)));
  test_cycle_perception(algorithm, "Figueras 19", "C1CC2CC3CCCC(C4)C3C5C2C6C1CC7C8C6C9C5C4CCC9CC8CCC7",
      TestCycle::list(TestCycle(38, 0, 1, 13, 14, 15, 16),
                      TestCycle(38, 2, 3, 10, 11, 12, 13),
                      TestCycle(38, 4, 5, 6, 7, 9, 10),
                      TestCycle(38, 8, 9, 11, 24, 25, 26),
                      TestCycle(38, 12, 14, 21, 22, 23, 24),
                      TestCycle(38, 15, 17, 18, 19, 20, 21),
                      TestCycle(38, 19, 33, 34, 35, 36, 37),
                      TestCycle(38, 20, 22, 30, 31, 32, 33),
                      TestCycle(38, 23, 25, 27, 28, 29, 30)));
  test_cycle_perception(algorithm, "Figueras 20", "C1C2C3C4CC5C46C37C28C1CC8C7C6C5",
      TestCycle::list(TestCycle(22, 0, 10, 11, 12),
                      TestCycle(22, 1, 8, 9, 10),
                      TestCycle(22, 2, 6, 7, 8),
                      TestCycle(22, 3, 4, 5, 6),
                      TestCycle(22, 5, 19, 20, 21),
                      TestCycle(22, 7, 17, 18, 19),
                      TestCycle(22, 9, 15, 16, 17),
                      TestCycle(22, 11, 13, 14, 15)));

  test_cycle_perception(algorithm, "Vismara, page 7", "C1C(CC2)CCC2CCCC(CC3)CCC31",
      TestCycle::list(TestCycle(18, 1, 2, 3, 4, 5, 6),
                      TestCycle(18, 11, 12, 13, 14, 15, 16),
                      TestCycle(18, 0, 1, 2, 6, 7, 8, 9, 10, 11, 12, 16, 17),
                      TestCycle(18, 0, 1, 2, 6, 7, 8, 9, 10, 13, 14, 15, 17),
                      TestCycle(18, 0, 3, 4, 5, 7, 8, 9, 10, 11, 12, 16, 17),
                      TestCycle(18, 0, 3, 4, 5, 7, 8, 9, 10, 13, 14, 15, 17)));

  test_cycle_perception(algorithm, "Horton, Figure 1: Counterexample", "C123C45C67C18C69C741C52C1983",
      TestCycle::list(TestCycle(17, 0, 10, 11),
                      TestCycle(17, 8, 9, 10),
                      TestCycle(17, 1, 7, 8),
                      TestCycle(17, 5, 6, 7),
                      TestCycle(17, 2, 4, 5),
                      TestCycle(17, 4, 14, 15),
                      TestCycle(17, 3, 15, 16),
                      TestCycle(17, 11, 12, 16),
                      TestCycle(17, 9, 12, 13),
                      TestCycle(17, 6, 13, 14)));

  test_cycle_perception(algorithm, "Berger, Figure 2 (a)", "C1C23CC24CC45CC531",
      TestCycle::list(TestCycle(12, 0, 10, 11),
                      TestCycle(12, 1, 2, 3),
                      TestCycle(12, 4, 5, 6),
                      TestCycle(12, 7, 8, 9),
                      TestCycle(12, 3, 6, 9, 10)));
  test_cycle_perception(algorithm, "Berger, Figure 2 (b)", "C1C23C45CC67C18C69C74C52C389",
      TestCycle::list(TestCycle(18, 0, 5, 15, 16),
                      TestCycle(18, 1, 12, 13),
                      TestCycle(18, 2, 3, 9, 10),
                      TestCycle(18, 4, 6, 7),
                      TestCycle(18, 7, 8, 9),
                      TestCycle(18, 10, 11, 12),
                      TestCycle(18, 13, 14, 15),
                      TestCycle(18, 6, 16, 17),
                      TestCycle(18, 8, 11, 14, 17)));
  test_cycle_perception(algorithm, "Berger, Figure 4 (a)", "C1C23C4C5C16C7C2C8C7C6C5C4C38",
      TestCycle::list(TestCycle(20, 0, 1, 2, 3, 4),
                      TestCycle(20, 0, 4, 5, 6, 7),
                      TestCycle(20, 1, 16, 17, 18),
                      TestCycle(20, 2, 14, 15, 16),
                      TestCycle(20, 3, 12, 13, 14),
                      TestCycle(20, 5, 10, 11, 12),
                      TestCycle(20, 6, 8, 9, 10),
                      TestCycle(20, 7, 8, 18, 19)));
  /* SLOW
  test_cycle_perception(algorithm, "Berger, Figure 4 (b)", "C1C234C56C78C19%10C%11%12C2%13C%11%14C%129C%107C85C63C4%13%14",
      TestCycle::list(TestCycle(26, 0, 1, 2, 3, 4),
                      TestCycle(26, 0, 4, 5, 6, 7),
                      TestCycle(26, 1, 20, 21),
                      TestCycle(26, 21, 22, 23),
                      TestCycle(26, 2, 17, 18),
                      TestCycle(26, 18, 19, 20),
                      TestCycle(26, 3, 14, 15),
                      TestCycle(26, 15, 16, 17),
                      TestCycle(26, 5, 11, 12),
                      TestCycle(26, 12, 13, 14),
                      TestCycle(26, 6, 8, 9),
                      TestCycle(26, 9, 10, 11),
                      TestCycle(26, 7, 23, 24),
                      TestCycle(26, 8, 24, 25)));
  */
  test_cycle_perception(algorithm, "Berger, Figure 5 (a)", "C12C34C1C25C6C5C67C3C74",
      TestCycle::list(TestCycle(15, 0, 1, 2),
                      TestCycle(15, 2, 3, 4),
                      TestCycle(15, 5, 6, 7),
                      TestCycle(15, 6, 8, 9),
                      TestCycle(15, 10, 12, 13),
                      TestCycle(15, 11, 12, 14),
                      TestCycle(15, 0, 4, 5, 9, 13, 14),
                      TestCycle(15, 0, 4, 5, 9, 10, 11),
                      TestCycle(15, 0, 4, 7, 8, 13, 14),
                      TestCycle(15, 0, 4, 7, 8, 10, 11),
                      TestCycle(15, 1, 3, 5, 9, 13, 14),
                      TestCycle(15, 1, 3, 5, 9, 10, 11),
                      TestCycle(15, 1, 3, 7, 8, 13, 14),
                      TestCycle(15, 1, 3, 7, 8, 10, 11)));
  test_cycle_perception(algorithm, "Berger, Figure 5 (b)", "C123C4C1CC45C6C5C67C2C73",
      TestCycle::list(TestCycle(16, 0, 1, 2),
                      TestCycle(16, 1, 3, 4, 5),
                      TestCycle(16, 6, 7, 8),
                      TestCycle(16, 7, 9, 10),
                      TestCycle(16, 11, 13, 14),
                      TestCycle(16, 12, 13, 15),
                      TestCycle(16, 0, 5, 8, 9, 11, 12),
                      TestCycle(16, 0, 5, 8, 9, 14, 15),
                      TestCycle(16, 0, 5, 6, 10, 11, 12),
                      TestCycle(16, 0, 5, 6, 10, 14, 15)));
  test_cycle_perception(algorithm, "Berger, Figure 7 (a)", "C12CCC3C1CCC32",
      TestCycle::list(TestCycle(10, 3, 4, 8, 9),
                      TestCycle(10, 0, 1, 2, 3, 4),
                      TestCycle(10, 0, 1, 2, 8, 9),
                      TestCycle(10, 3, 5, 6, 7, 8),
                      TestCycle(10, 4, 5, 6, 7, 9)));
  test_cycle_perception(algorithm, "Berger, Figure 7 (b)", "C12CC3CC4CC(C1)C56C47C38C25C9CC8CC7CC6C9",
      TestCycle::list(TestCycle(28, 0, 1, 12, 13, 14),
                      TestCycle(28, 2, 3, 10, 11, 12),
                      TestCycle(28, 4, 5, 8, 9, 10),
                      TestCycle(28, 6, 7, 8, 14, 15),
                      TestCycle(28, 9, 22, 23, 24, 25),
                      TestCycle(28, 11, 19, 20, 21, 22),
                      TestCycle(28, 13, 16, 17, 18, 19),
                      TestCycle(28, 15, 16, 25, 26, 27),
                      TestCycle(28, 9, 11, 13, 15)));
  test_cycle_perception(algorithm, "Berger, Figure 8 (a)", "C12C3C4CC5CC6CC1C7C6C5C4CC3CC2C7",
      TestCycle::list(TestCycle(24, 0, 18, 19, 20, 21),
                      TestCycle(24, 1, 15, 16, 17, 18),
                      TestCycle(24, 2, 3, 13, 14, 15),
                      TestCycle(24, 4, 5, 11, 12, 13),
                      TestCycle(24, 6, 7, 9, 10, 11),
                      TestCycle(24, 8, 9, 21, 22, 23),
                      TestCycle(24, 0, 1, 8, 9, 10, 12, 14, 15)));
  test_cycle_perception(algorithm, "Berger, Figure 8 (b)", "C12CCCC3C4C1C56CC47CC38CC2(C5)C876",
      TestCycle::list(/*TestCycle(23, 0, 1, 2, 3, 4, 5, 6),*/
                      TestCycle(23, 4, 10, 11, 12, 13),
                      TestCycle(23, 4, 10, 13, 20, 21),
                      TestCycle(23, 5, 7, 8, 9, 10),
                      TestCycle(23, 5, 7, 10, 21, 22),
                      TestCycle(23, 6, 7, 16, 17, 18),
                      TestCycle(23, 6, 7, 16, 19, 22),
                      TestCycle(23, 8, 9, 21, 22),
                      TestCycle(23, 11, 12, 20, 21),
                      TestCycle(23, 14, 15, 19, 20),
                      TestCycle(23, 17, 18, 19, 22)));

  //test_cycle_perception(algorithm, "nanotube", "[C@@H]12C[C@H]3C[C@H]4C[C@H]5C[C@H]6C[C@H]7C[C@H]8C[C@H]9C[C@H]%10C[C@H]%11C[C@H]%12C[C@H]%13C[C@H]%14C[C@@H](C1)[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@@H]3[C@H]2[C@H]1[C@H]1[C@H]2[C@@H]3[C@H]3[C@@H]4[C@H]4[C@@H]5[C@H]5[C@@H]6[C@H]6[C@@H]7[C@H]7[C@@H]8[C@H]8[C@@H]9[C@H]9[C@@H]%10[C@H]%10[C@@H]%11[C@H]%11[C@@H]%12[C@H]%12[C@@H]%13[C@H]%13[C@@H]%14[C@H]%14[C@@H]%15[C@@H]1[C@@H]1[C@@H]%15[C@H]%14[C@@H]%14[C@H]%13[C@@H]%13[C@H]%12[C@@H]%12[C@H]%11[C@@H]%11[C@H]%10[C@@H]%10[C@H]9[C@@H]9[C@H]8[C@@H]8[C@H]7[C@@H]7[C@H]6[C@@H]6[C@H]5[C@@H]5[C@H]4[C@@H]4[C@H]3[C@H]([C@H]2C1)[C@H](C)C[C@@H]4C[C@@H]5C[C@@H]6C[C@@H]7C[C@@H]8C[C@@H]9C[C@@H]%10C[C@@H]%11C[C@@H]%12C[C@@H]%13C[C@@H]%14CC%15", TestCycle::list());

  /*
  test_cycle_perception(algorithm, "Figueras 16", "",
      TestCycle::list(TestCycle(1, ),
                      TestCycle(1, ),
                      TestCycle(1, )));
  */


}



