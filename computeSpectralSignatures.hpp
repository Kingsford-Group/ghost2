#include <iostream>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "progressBar.hpp"
#include "graph.hpp"
#include "laplacian.hpp"
#include "gzWriter.hpp"

#include <boost/thread.hpp>
#include "threadpool.hpp"

using boost::threadpool::pool;
typedef boost::posix_time::microsec_clock bclock;
typedef boost::posix_time::ptime ptime;

void spectrum(Graph *input, int numHops, string source, 
              vector<LevelData> *output, ProgressBar *p)
{
  AdjacencyMatrix m;
  for(int i=1;i<=numHops;i++)
  {
    m.extend(input, source);
    (*output).push_back(LevelData(m.getPrev(), m.getEigen(), m.getDensity()));
  }
  (*p).update();
}

void computeSpectralSignatures(Graph *input, int numHops, int numP)
{
  ptime t = bclock::local_time();
  if(numP < 1) numP = boost::thread::hardware_concurrency();
  if(numP < 1) numP = 2;
  pool tpool(numP);
  vector<string> nodes = (*input).nodes();
  int numNodes = nodes.size();
  GzWriter w((*input).getName()+".sig.gz");
  ProgressBar p(numNodes);
  cout << "creating: " << (*input).getName() << ".sig.gz\n";
  w.writeInt(numNodes);
  w.writeInt(numHops);
  
  levelMap levelmap;
  vector< vector<LevelData> > data;
  data.resize(numNodes);
  for(int i=0;i<numNodes;i++)
    tpool.schedule(
      boost::bind(&spectrum, input, numHops, nodes[i], &(data[i]), &p)
    );
  tpool.wait();
  cout << "\n";

  for(int i=0;i<numNodes;i++)
    levelmap[nodes[i]] = data[i];
  w.writeData(levelmap);
  cout << "created: " << (*input).getName() << ".sig.gz in " <<
    (bclock::local_time() - t).total_milliseconds() << " milliseconds\n";
}

