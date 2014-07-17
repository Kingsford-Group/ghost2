#include <iostream>
#include <string>
#include <stdlib.h>
#include "config.hpp"
#include "graph.hpp"
#include "readFromGexf.hpp"
#include "readFromNet.hpp"
#include "computeSpectralSignatures.hpp"
#include "blastDistance.hpp"
#include "spectralToDistance.hpp"
#include "alignGraphs.hpp"
#include "localImprove.hpp"
//#include "localImprove2.hpp"

using std::string;
using std::cout;

bmap readAlignment(string file)
{
  bmap result;
  ifstream fin(file);
  char c = fin.get();
  int alen = 0;
  int blen = 0;
  char a[128];
  char b[128];
  bool fillinga = true;
  while(fin.good()){
    if(c == '\n'){
      b[blen] = '\0';
      result.insert(bmap::value_type(string(a), string(b)));
      alen = blen = 0;
      fillinga = true;
    }
    else if(c == '\r'){}
    else if (c == '\t'){
      a[alen] = '\0';
      fillinga = false;
    }
    else if(fillinga) a[alen++] = c;
    else b[blen++] = c;
    c = fin.get();
  }
  return result;
}

void printMap2(vector<D_alpha> &f, string filename)
{
  ofstream fout (filename);
  auto iter = f.begin(),
  iend = f.end();
  for(; iter != iend; ++iter)
  {
    fout << iter->get_n1() << "\t" << iter->get_n2() << "\n";
  }
  fout.close();
}

void computeAlignment(ConfigData c)
{
  // read in graph
  Graph G,H;
  if(c.Ggraph.substr(c.Ggraph.size()-5) == ".gexf")
    G = readFromGexf(c.Ggraph, c.directed);
  else if(c.Ggraph.substr(c.Ggraph.size()-4) == ".net")
    G = readFromNet(c.Ggraph, c.directed);
  else {cout << "bad extension: " << c.Ggraph << "\n"; exit(0);}

  if(c.Hgraph.substr(c.Hgraph.size()-5) == ".gexf")
    H = readFromGexf(c.Hgraph, c.directed);
  else if(c.Hgraph.substr(c.Hgraph.size()-4) == ".net")
    H = readFromNet(c.Hgraph, c.directed);
  else {cout << "bad extension: " << c.Hgraph << "\n"; exit(0);}

  if(c.AlignFile != "")
  {
    bmap align = readAlignment(c.AlignFile);
    printICS(G, H, align);
    if(c.searchiter != 0){
      // get evalues if given
      blastMap *evals = new blastMap;
      if(c.SeqScores == "")
        evals = NULL;
      else
        *evals = getBlastMap(c.SeqScores);
      localImprove(G, H, evals, &align, c.searchiter, c.ratio, c.numProcessors);
      //localImprove2(G, H, &align, c.searchiter, c.numProcessors);
      printICS(G, H, align);
    }
    return;
  }

/*bmap f2 = alignQ(G, H);
  printICS(G, H, f2);
  localImprove(G, H, NULL, &f2, c.searchiter, c.ratio, c.numProcessors);
  if(c.searchiter != 0)printICS(G, H, f2);
  printMap(f2, G.getName(), H.getName());
  return;*/

  vector<D_alpha> dist;

  if(c.DistFile == ""){
  // compute spectral signatures
  if(c.Gsigs == "")
  {
    computeSpectralSignatures(&G, c.hops, c.numProcessors);
    c.Gsigs = (G.getName() + ".sig.gz");
  }
  else
  {
    unsigned start = c.Gsigs.find_last_of("/");
    if(start == string::npos) start = -1;
    string name = c.Gsigs.substr(start+1);
    if(name != (G.getName() + ".sig.gz"))
    {
      cout << "file: " << c.Gsigs << " does match the network file name.\n" <<
      "please correct before continuing.";
      exit(0);
    }
  }
  if(c.Hsigs == "")
  {
    computeSpectralSignatures(&H, c.hops, c.numProcessors);
    c.Hsigs = (H.getName() + ".sig.gz");
  }
  else 
  {
    unsigned start = c.Hsigs.find_last_of("/");
    if(start == string::npos) start = -1;
    string name = c.Hsigs.substr(start+1);
    if(name != (H.getName() + ".sig.gz"))
    {
      cout << "file: " << c.Hsigs << " does match the network file name.\n" <<
      "please correct before continuing.";
      exit(0);
    }
  }

  if(c.dumpSignatures) return; // if user only wanted sigs
  }

  // undirect for later steps
  G.direct(false); H.direct(false);

  // get evalues if given
  blastMap *evals = new blastMap;
  if(c.SeqScores == "")
    evals = NULL;
  else
    *evals = getBlastMap(c.SeqScores);

  if(c.DistFile == ""){
  // compute distances
  dist = 
    getDistances(c.Gsigs, c.Hsigs, (G.getName()+"_vs_"+H.getName()+".sdf"), 
                 c.alpha, c.beta, evals , c.numProcessors);
  // if user wanted just the distances...
  if(c.dumpDistances) { delete evals; return; }
  }else{
    dist = getDistancesFromFile(c.DistFile, c.alpha, c.beta, evals);
  }

 /* cout << "starting alignment\n";
  vector<D_alpha> ans = performMatching(dist);
  cout << "done, writing results\n";
  string file = G.getName() + "_vs_" + H.getName() + ".af";
  printMap2(ans, file);
  bmap align = readAlignment(file);
  printICS(G, H, align);*/

  // align graphs
  bmap f = alignGraphs(G, H, dist, c.nneighbors, c.seedSkip);
  localImprove(G, H, evals, &f, c.searchiter, c.ratio, c.numProcessors);
  printICS(G, H, f);
  delete evals;
  printMap(f, G.getName(), H.getName());
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

  if(c.Ggraph == "" || c.Hgraph == "") // required input
    { cout << "gexf files not provided\n"; return 0; }

  // and here we go!
  computeAlignment(c);
  return 0;
}

