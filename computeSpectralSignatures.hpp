#include <iostream>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
#include "graph.hpp"
#include "laplacian.hpp"
#include "gzWriter.hpp"

#include <boost/thread/thread.hpp>
#include "threadpool.hpp"

using boost::threadpool::pool;

void spectrum(Graph *input, int numHops, string source, 
              vector<LevelData> *output)
{
  AdjacencyMatrix m;
  for(int i=1;i<=numHops;i++)
  {
    m.extend(input, source);
    (*output).push_back(LevelData(m.getPrev(), m.getEigen(), m.getDensity()));
  }
}

void computeSpectralSignatures(Graph *input, int numHops, int numP)
{
  pool tpool(numP);
  vector<string> nodes = (*input).nodes();
  int numNodes = nodes.size();
  GzWriter w((*input).getName()+".sig.gz");
  w.writeInt(numNodes);
  w.writeInt(numHops);
  
  umap levelmap;
  vector< vector<LevelData> > data;
  data.resize(numNodes);
  for(int i=0;i<numNodes;i++)
    tpool.schedule(
      boost::bind(&spectrum, input, numHops, nodes[i], &(data[i]))
    );
  tpool.wait();

  for(int i=0;i<numNodes;i++)
    levelmap[nodes[i]] = data[i];
  w.writeData(levelmap);
}

