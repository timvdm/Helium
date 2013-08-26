# Starting the service

The service is a single executable heliumservice, it takes the following command line options:

- --port The port the service should listen on ( default is 8080 )
- --sinindex The similarity fingerprint index filepath
- --subindex The substructure fingerprint index filepath
- --mol The molecule filepath

Example usage:

    ./heliumservice --port 8089 --simindex ~/tmp/pubchem3M.simfp.hel --subindex ~/tmp/pubchem3M.subfp.hel --mol ~/tmp/pubchem3M.mols.hel

The service will take a few seconds to load the indexes at which point it will be ready to processes requests.

# Using the service

The service is accessed using HTTP GET request, for example using the command line tool curl. The URIs are of the following form:

    /helium/<similarity|substructure>?q=<query_string>&pretty=<true|false>
    
Where:

- query_string - The SMILES string being searched for.
- pretty - An optional parameter ( defaults to false ). If set to true the JSON returned will be formatted.
- limit - An optional parameter to limit the search results returned. The default limit is 100.

 
### Similarity search example:

    curl 'localhost:8089/helium/similarity?q=c1ccccc1Cl&pretty=true&limit=10'

Will return the following JSON:

    {
      "hits" : [
        {
          "index" : 11614,
          "tanimoto" : 0.7222222222222222
        },
             {
          "index" : 123476,
          "tanimoto" : 0.7222222222222222
        },

        ...
        
      ]
    }

### Substructure search example:

    curl  'localhost:8089/helium/substructure?q=c1cccc(cc(CCN)n2)c12&pretty=true'    

Will return the following JSON:

    {
      "confirmed" : 4215,
      "false_positives" : 0.9144110301135094,
      "hits" : [
      {
         "index" : 2009
      },
      {
         "index" : 3189
      },
      
      ...
      
      ],
      "screened" : 49247
    }
      
      
