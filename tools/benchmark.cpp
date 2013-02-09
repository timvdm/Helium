#include "../src/molecule.h"
#include "../src/fileio.h"
#include "../src/bitvec.h"

#include "../src/enumeratesubgraphs.h"
#include "../src/substructure.h"
#include "../src/extendedconnectivities.h"
#include "../src/canonical.h"

#include <numeric>

#include <boost/timer/timer.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "args.h"

const boost::timer::nanosecond_type one_milisecond(1000000L);

using namespace Helium;

std::string get_filename(const std::string &filename)
{
  return DATADIR + std::string("/") + filename;
}

using namespace boost::gregorian;
using namespace boost::posix_time;

class BenchmarkRecorder
{
  public:

    BenchmarkRecorder() : m_first(true)
    {
      std::stringstream filename;

      date d(day_clock::local_day());
      filename << d.month() << "_" << d.day() << "_" << d.year() << "_";

      ptime t(microsec_clock::local_time());
      time_duration now = t.time_of_day();
      filename << now.hours() << "_" << now.minutes() << "_" << now.seconds();

      filename << ".benchmark";

      m_ofs.open(filename.str().c_str());

      m_ofs << "{ \"time\": \"" << d.month() << " " << d.day() << " " << d.year() << ", ";
      if (now.hours() < 10)
        m_ofs << 0;
      m_ofs << now.hours() << ":";
      if (now.minutes() < 10)
        m_ofs << 0;
      m_ofs << now.minutes() << ":";
      if (now.seconds() < 10)
        m_ofs << 0;
      m_ofs << now.seconds() << "\", " << "\"benchmarks\": [";
    }

    ~BenchmarkRecorder()
    {
      m_ofs << "]}";
    }

    void record(const std::string &label, int ms)
    {
      if (m_first)
        m_first = false;
      else
        m_ofs << ", ";

      m_ofs << "{ \"benchmark\": \"" << label << "\", \"ms\": " << ms << "}";
    }

  private:
    std::ofstream m_ofs;
    bool m_first;
};

BenchmarkRecorder recorder;


/**
 *
 * Benchmark reading *.hem files.
 *
 */

void benchmark_read(const std::string &filename)
{
  boost::timer::cpu_timer timer;

  // open molecule file
  MoleculeFile file(get_filename(filename));
  HeMol mol;

  // read molecules
  while (file.read_molecule(mol));

  boost::timer::cpu_times elapsed = timer.elapsed();
  int ms = (elapsed.system + elapsed.user) / one_milisecond;
  std::cout << "read_file (" << file.current() << " molecules): " << ms << " ms" << std::endl;
  recorder.record("read_file", ms);
}

/**
 *
 * Benchmark subgraph enumeration.
 *
 */

template<typename MoleculeType>
struct EnumerateSubgraphsCallback
{
  EnumerateSubgraphsCallback() : count(0) {}
  void operator()(const Subgraph &subgraph) { ++count; }
  int count;
};

void benchmark_enumerate_subgraphs(const std::string &filename, int size = 7)
{
  boost::timer::cpu_timer timer;

  // open molecule file
  MoleculeFile file(get_filename(filename));
  HeMol mol;

  // read molecules
  while (file.read_molecule(mol)) {
    EnumerateSubgraphsCallback<HeMol> callback;
    enumerate_subgraphs(mol, callback, size);
  }

  boost::timer::cpu_times elapsed = timer.elapsed();
  int ms = (elapsed.system + elapsed.user) / one_milisecond;
  std::cout << "subgraph_enumeration (" << file.current() << " molecules): " << ms << " ms" << std::endl;
  recorder.record("subgraph_enumeration", ms);
}

void benchmark_canonicalize(const std::string &filename)
{
  boost::timer::cpu_timer timer;

  // open molecule file
  MoleculeFile file(get_filename(filename));
  HeMol mol;

  // read molecules
  while (file.read_molecule(mol)) {
    canonicalize(mol, extended_connectivities(mol));
  }

  boost::timer::cpu_times elapsed = timer.elapsed();
  int ms = (elapsed.system + elapsed.user) / one_milisecond;
  std::cout << "canonicalize (" << file.current() << " molecules): " << ms << " ms" << std::endl;
  recorder.record("canonicalize", ms);
}







int main(int argc, char**argv)
{
  benchmark_read("1M.hem");
  benchmark_enumerate_subgraphs("10K.hem");
  benchmark_canonicalize("10K.hem");
}
