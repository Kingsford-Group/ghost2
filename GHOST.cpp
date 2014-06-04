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
  Graph G = readFromGexf(c.Ggexf);
  cout << "done reading in G\n";
  Graph H = readFromGexf(c.Hgexf);
  cout << "done reading in H\n";
  if(c.Gsigs == "")
  {
    computeSpectralSignatures(&G, c.hops, c.numProcessors);
    c.Gsigs = (G.getName() + ".sig.gz");
  }
  cout << "done writing sigs for G\n";
  if(c.Hsigs == "")
  {
    computeSpectralSignatures(&H, c.hops, c.numProcessors);
    c.Hsigs = (H.getName() + ".sig.gz");
  }
  cout << "done writing sigs for H\n";
  blastMap *evals = new blastMap;
  if(c.SeqScores == "")
    evals = NULL;
  else
    *evals = getBlastMap(c.SeqScores);
  vector<D_alpha> dist = 
    getDistances(c.Gsigs, c.Hsigs, (G.getName()+"_vs_"+H.getName()+".sdf"), 
                 c.alpha, evals);
  delete evals;
  cout << "done writing distances\n";
  alignGraphs(G, H, dist, c.beta, c.nneighbors);
  cout << "DONE\n";
}

int main(int argc, char** argv)
{
  ConfigData c;
  for(int i=1; i<argc; i++)
  {
    if(string(argv[i]) == "-c") 
    { 
      if((i+1)<argc) 
        c.configure(string(argv[i+1])); 
    }
  }

  if(c.Ggexf == "" || c.Hgexf == "")
  {
    cout << "gexf files not provided\n";
    return 0;
  }
  computeAlignment(c);
  return 0;
}

