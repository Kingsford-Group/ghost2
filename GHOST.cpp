#include <iostream>
#include <string>
#include <stdlib.h>
#include "config.hpp"
#include "graph.hpp"
#include "readFromGexf.hpp"
#include "computeSpectralSignatures.hpp"
#include "blastDistance.hpp"
#include "spectralToDistance.hpp"
#include "alignGraphs.hpp"

using std::string;
using std::cout;

void computeAlignment(ConfigData c)
{
  // read in graph
  Graph G = readFromGexf(c.Ggexf);
  Graph H = readFromGexf(c.Hgexf);

  // compute spectral signatures
  if(c.Gsigs == "")
  {
    computeSpectralSignatures(&G, c.hops, c.numProcessors);
    c.Gsigs = (G.getName() + ".sig.gz");
  }
  if(c.Hsigs == "")
  {
    computeSpectralSignatures(&H, c.hops, c.numProcessors);
    c.Hsigs = (H.getName() + ".sig.gz");
  }
  if(c.dumpSignatures) return; // if user only wanted sigs

  // get evalues if given
  blastMap *evals = new blastMap;
  if(c.SeqScores == "")
    evals = NULL;
  else
    *evals = getBlastMap(c.SeqScores);

  // compute distances
  vector<D_alpha> dist = 
    getDistances(c.Gsigs, c.Hsigs, (G.getName()+"_vs_"+H.getName()+".sdf"), 
                 c.alpha, evals , c.numProcessors);
  if(c.dumpDistances) return; // if user wanted just the distances...
  delete evals;

  // align graphs
  alignGraphs(G, H, dist, c.beta, c.nneighbors);
}

int main(int argc, char** argv)
{
  ConfigData c;

  // read in options
  for(int i=1; i<argc; i++)
  {
    if(string(argv[i]) == "-c") 
      if((i+1)<argc) 
        c.configure(string(argv[i+1]));
    if(string(argv[i]) == "-k")
      if((i+1)<argc)
        c.hops = atoi(argv[i+1]);
    if(string(argv[i]) == "-p")
      if((i+1)<argc)
        c.numProcessors = atoi(argv[i+1]);
  }

  if(c.Ggexf == "" || c.Hgexf == "") // required input
    { cout << "gexf files not provided\n"; return 0; }

  // and here we go!
  computeAlignment(c);
  return 0;
}

