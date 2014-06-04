#include <iostream>
#include <string>
#include <stdlib.h>
#include "config.hpp"
#include "graph.hpp"
#include "readFromGexf.hpp"
//#include "computeSpectralSignatures.h"
#include "spectralToDistance.hpp"
#include "alignGraphs.hpp"

using std::string;
using std::cout;

void computeAlignment(ConfigData c)
{
  Graph G = readFromGexf(c.Ggexf);
  Graph H = readFromGexf(c.Hgexf);
  G.print();
  H.print();
//  computeSpectralSignatures(&G, c.hops, c.numProcessors);
//  computeSpectralSignatures(&H, c.hops, c.numProcessors);
//  spectralToDistance(G,H); // done by david
//  alignGraphs(G,H); // done by david
}

int main(int argc, char** argv)
{
  Graph G = readFromGexf("scerehc.gexf");
  Graph H = readFromGexf("scere05.gexf");
  auto dists = getDistances("scerehc.sig.gz","scere05.sig.gz", "scerehc_vs_scere05.dist", 1, NULL);
  alignGraphs(G,H, dists, 100, 3000);

  /*ConfigData c;
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
  computeAlignment(c);*/
  return 0;
}

