#ifndef HELIUM_RESTFULSERVICE_H
#define HELIUM_RESTFULSERVICE_H

#include <Helium/fileio/fingerprints.h>
#include <Helium/fileio/moleculefile.h>
#include <Helium/fingerprints/similarity.h>
#include <Helium/bitvec.h>

#include <string>

struct mg_connection;

namespace Helium {

class RESTfulService
{
public:
  RESTfulService();
  ~RESTfulService();
  void init(const std::string &similarityIndex,
      const std::string &subStructureIndex, const std::string &moleculeFile);
  int handleRequest(struct mg_connection *conn);
  void logInfo(const std::string &msg);
  void logError(const std::string &msg);

private:
  InMemoryRowMajorFingerprintStorage m_similarityStorage;
  MemoryMappedMoleculeFile m_moleculeFile;
  SimilaritySearchIndex<Helium::InMemoryRowMajorFingerprintStorage> *m_similarityIndex;

  InMemoryColumnMajorFingerprintStorage m_subStructureStorage;

  std::string createResponse(const std::string &status, const std::string &type,
      const std::string &content);
  Word* computeFingerprint(const std::string &settings,
      const std::string &smiles);
  Word* computeFingerprint(const std::string &settings,
      HeMol &query);

  std::string similaritySearch(const std::string &query, bool pretty,
      unsigned int limit);
  std::string subStructureSearch(const std::string &query, bool pretty,
      unsigned int limit);
  void screen(Word *fingerprint, Word *result);
};

}

#endif
