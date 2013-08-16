#include "restfulservice.h"
#include "mongoose/mongoose.h"

#include <Helium/smiles.h>
#include <Helium/algorithms/isomorphism.h>

#include <signal.h>

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iterator>

#ifdef _WIN32
#define WINCDECL __cdecl
#else
#define WINCDECL
#endif

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/lexical_cast.hpp>

using namespace std;

// Have to define global instance as callback is C
static Helium::RESTfulService service;

// This function will be called by mongoose on every new request.
static int begin_request_handler(struct mg_connection *conn)
{
  return service.handleRequest(conn);
}

namespace Helium {

static string PREFIX = "/helium/";
static string SIMILARITY = "similarity";
static string SUBSTRUCTURE = "substructure";
static double TMIN =  0.7 - 10e-5;

RESTfulService::RESTfulService()
  : m_similarityIndex(NULL)
{

}

RESTfulService::~RESTfulService()
{
  delete m_similarityIndex;
}

string RESTfulService::createResponse(const string &status, const string &type,
    const string &content)
{
  ostringstream response;
  response << "HTTP/1.1 ";
  response << status << "\r\n";
  response << "Content-Type: ";
  response << type << "\r\n";
  response << "Content-Length: ";
  response << content.size() << "\r\n";
  response << "\r\n";
  response << content;

  return response.str();
}

int RESTfulService::handleRequest(struct mg_connection *conn)
{
  const struct mg_request_info *request_info = mg_get_request_info(conn);

  string method = request_info->request_method;
  string uri = request_info->uri;
  string queryString = request_info->query_string;

  // Is this a supported method
  if (method != "GET") {
    string response = createResponse("405 Method Not Allowed", "text/plain",
        "Unsupported method");
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }

  // Now parse the URI
  if (uri.find(PREFIX) == string::npos) {
    string response = createResponse("404 Not Found", "text/plain",
        "Not Found");
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }

  // Get the query var
  char buf[1024];

  int length = mg_get_var(request_info->query_string, queryString.length(),
      "q", buf, 1024);
  if (length <= 0) {
    string response = createResponse("400 Bad Request", "text/plain",
        "A query value must be provided");
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }

  string query(buf, length);

  // Get the pretty var
  length = mg_get_var(request_info->query_string, queryString.length(),
      "pretty", buf, 1024);

  bool pretty = false;
  if (length == 0) {
    pretty = true;
  }
  else if (length > 0) {
    pretty = string(buf, length) == "true";
  }

  // Substructure
  if (uri.find(SUBSTRUCTURE, PREFIX.length()) == PREFIX.length()) {

    // Are the necessary files loaded
    if (m_subStructureStorage.numFingerprints() == 0
        || m_moleculeFile.numMolecules() == 0) {
      string response = createResponse("503 Bad  Service Unavailable", "text/plain",
           "Service unavailable while indexes are loading");
      mg_write(conn, response.c_str(), response.length());
      return 1;
    }

    string result = subStructureSearch(query, pretty);
    string response = createResponse("200 OK", "application/json", result);
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }
  // Similarity
  else if (uri.find(SIMILARITY, PREFIX.length()) == PREFIX.length()) {

    // Is the index loaded?
    if (m_similarityIndex == NULL) {
      string response = createResponse("503 Bad  Service Unavailable", "text/plain",
           "Service unavailable while indexes are loading");
      mg_write(conn, response.c_str(), response.length());
      return 1;
    }

    string result = similaritySearch(query, pretty);
    string response = createResponse("200 OK", "application/json", result);
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }
  else {
    string response = createResponse("404 Not Found", "text/plain",
          "Supported searches are \"similarity\" and \"substructure\"");
    mg_write(conn, response.c_str(), response.length());
    return 1;
  }

}

Word* RESTfulService::computeFingerprint(const std::string &settings,
    const std::string &smiles)
{
  HeMol mol;
  parse_smiles(smiles, mol);

  return computeFingerprint(settings, mol);
}

Word* RESTfulService::computeFingerprint(const std::string &settings,
    HeMol &mol)
{
  Json::Reader reader;
  Json::Value data;

  if (!reader.parse(settings, data)) {
    std::cerr << reader.getFormattedErrorMessages() << std::endl;
    return 0;
  }

  int bits = data["num_bits"].asInt();
  int words = bitvec_num_words_for_bits(bits);
  int k = data["fingerprint"]["k"].asInt();
  int prime = data["fingerprint"]["prime"].asInt();
  std::string type = data["fingerprint"]["type"].asString();

  Word *fingerprint = new Word[words];

  if (type == "Helium::paths_fingerprint") {
    path_fingerprint(mol, fingerprint, k, words, prime);
    return fingerprint;
  }

  if (type == "Helium::trees_fingerprint") {
    tree_fingerprint(mol, fingerprint, k, words, prime);
    return fingerprint;
  }

  if (type == "Helium::subgraph_fingerprint") {
    subgraph_fingerprint(mol, fingerprint, k, words, prime);
    return fingerprint;
  }

  logError("Fingerprint type \"" + type + "\" not recognised");

  delete [] fingerprint;
  return 0;
}

string RESTfulService::similaritySearch(const string &query, bool pretty)
{
  Word *fingerPrint = computeFingerprint(m_similarityStorage.header(), query);

  vector<pair<unsigned int, double> > result
    = m_similarityIndex->search(fingerPrint, TMIN);

  Json::Value data;
    for (std::size_t j = 0; j < result.size(); ++j) {
      data["hits"][Json::ArrayIndex(j)] = Json::Value(Json::objectValue);
      Json::Value &obj = data["hits"][Json::ArrayIndex(j)];
      obj["index"] = result[j].first;
      obj["tanimoto"] = result[j].second;
    }

    string content;

    if (pretty) {
      Json::StyledWriter writer;
      content = writer.write(data);
    } else {
      Json::FastWriter writer;
      content = writer.write(data);
    }

  delete[] fingerPrint;

  return content;
}


void RESTfulService::screen(Word *fingerprint, Word *result)
{
  bool first = true;
  for (int i = 0; i < m_subStructureStorage.numBits(); ++i) { // foreach bit
    // skip this bit if it is not set
    if (!bitvec_get(i, fingerprint))
      continue;

    if (first) {
      // if this is the first bit, just set result
      for (unsigned int j = 0; j < bitvec_num_words_for_bits(m_subStructureStorage.numFingerprints()); ++j)
        result[j] = m_subStructureStorage.bit(i)[j];
      first = false;
    } else {
      // do set intersection
      for (unsigned int j = 0; j < bitvec_num_words_for_bits(m_subStructureStorage.numFingerprints()); ++j)
        result[j] &= m_subStructureStorage.bit(i)[j];
    }
  }
}

string RESTfulService::subStructureSearch(const string &smiles, bool pretty)
{
  HeMol query;
  parse_smiles(smiles, query);

  // compute query fingerprint
  Word *queryFingerprint = computeFingerprint(m_subStructureStorage.header(),
                                              query);

  // perform search
  Word *candidates = new Word[bitvec_num_words_for_bits(
                                m_subStructureStorage.numFingerprints())];

  screen(queryFingerprint, candidates);

  HeMol mol;
  std::vector<Index> result;
  int can =0;
  for (unsigned int i = 0; i < m_moleculeFile.numMolecules(); ++i) {
    if (bitvec_get(i, candidates)) {
      ++can;
      m_moleculeFile.read_molecule(i, mol);
      if (isomorphism_search<DefaultAtomMatcher, DefaultBondMatcher>(mol, query))
        result.push_back(i);
    }
  }

  // print results
  Json::Value data;
  data["hits"] = Json::Value(Json::arrayValue);
  data["screened"] = Json::Value(bitvec_count(candidates,
                                 bitvec_num_words_for_bits(
                                     m_subStructureStorage.numFingerprints())));
  data["confirmed"] = Json::Value(static_cast<unsigned int>(result.size()));
  data["false_positives"] = Json::Value(1.0 - static_cast<double>(result.size())
      / bitvec_count(candidates,
          bitvec_num_words_for_bits(m_subStructureStorage.numFingerprints())));
  for (std::size_t i = 0; i < result.size(); ++i) {
    data["hits"][Json::ArrayIndex(i)] = Json::Value(Json::objectValue);
    Json::Value &obj = data["hits"][Json::ArrayIndex(i)];
    obj["index"] = result[i];
  }

  // deallocate fingerprint
  if (queryFingerprint)
    delete [] queryFingerprint;
  delete [] candidates;

  string content;

  if (pretty) {
    Json::StyledWriter writer;
    content = writer.write(data);
  } else {
    Json::FastWriter writer;
    content = writer.write(data);
  }

  return content;
}

void RESTfulService::init(const string &similarityIndex,
    const string &subStructureIndex, const string &moleculeFile)
{
  m_similarityStorage.load(similarityIndex);
  m_similarityIndex
    = new Helium::SimilaritySearchIndex<Helium::InMemoryRowMajorFingerprintStorage>(m_similarityStorage, 3);
  logInfo("Similarity index loaded: " + similarityIndex);
  m_subStructureStorage.load(subStructureIndex);
  logInfo("Substructure index loaded: " + subStructureIndex);
  m_moleculeFile.load(moleculeFile);
  logInfo("Molecules loaded: " + moleculeFile);
}

void RESTfulService::logInfo(const string &msg)
{
  cout << "[INFO] " << msg << endl;
}

void RESTfulService::logError(const string &msg)
{
  cerr << "[ERROR] " << msg << endl;
}

}

static int signal_number;

static void WINCDECL signal_handler(int sig_num) {
  signal_number = sig_num;
}

int main(int argc, char* argv[]) {

  // Process arguments
  int port;
  po::options_description desc("Options");
  desc.add_options()
    ("help", "produce help message")
    ("port", po::value<int>(&port)->default_value(8080), "the port to listen on")
    ("simindex", po::value<string>(), "similarity fingerprint index file")
    ("subindex", po::value<string>(), "substructure fingerprint index file")
    ("mol", po::value<string>(), "molecule file");

  po::variables_map vm;
  po::store(parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("simindex") || !vm.count("subindex")
      || !vm.count("mol")) {
    cerr << desc << endl;
    return 1;
  }

  string similarityIndex = vm["simindex"].as<string>();
  string subStructureIndex = vm["subindex"].as<string>();
  string molecules = vm["mol"].as<string>();

  struct mg_context *ctx;
  struct mg_callbacks callbacks;

  string portString = boost::lexical_cast<string>(port);

  // List of options. Last element must be NULL.
  const char *options[] = {"listening_ports", portString.c_str(),
                           NULL};

  string version = mg_version();
  service.logInfo("Mongoose version " + version);

  service.logInfo("Service listening on: " + portString);

  // Prepare callbacks structure. We have only one callback, the rest are NULL.
  memset(&callbacks, 0, sizeof(callbacks));
  callbacks.begin_request = begin_request_handler;

  // Start the web server.
  ctx = mg_start(&callbacks, NULL, options);

  if (!ctx) {
    service.logError("Error starting service on port: " + portString
        + " (is this port already is use?)");
    return 1;
  }

  // Load the indexes
  try {
    service.init(similarityIndex, subStructureIndex, molecules);
  } catch (const exception &e) {
    service.logError("Error loading index: " + string(e.what()));
    return 1;
  }

  service.logInfo("Ready to process requests");
  // Wait until user hits "enter". Server is running in separate thread.
  // Navigating to http://localhost:8080 will invoke begin_request_handler().

  signal(SIGTERM, signal_handler);
  signal(SIGINT, signal_handler);

  while (signal_number == 0)
    sleep(1);

  service.logInfo("Exiting with signal: "
      + boost::lexical_cast<string>(signal_number));

  // Stop the server.
  mg_stop(ctx);

  return 0;
}

