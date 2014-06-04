#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <tgmath.h>

#include <boost/unordered_map.hpp>
#include <boost/thread.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

#include "threadpool.hpp"
#include "fastonebigheader.h"

#include "gzreader.h"
#include "dalpha.h"

using std::vector;
using std::string;
using std::cout;
using std::pair;

typedef boost::posix_time::microsec_clock bclock;
typedef boost::posix_time::ptime ptime;
typedef boost::unordered_map<pair<string, string>, double> blastmap;

//Helpers for D_alpha calculation
double klDiv(double *p1, double *p2, int count)
{
  double s = 0;
  for(int i = 0; i < count; i++)
    if(p1[i] * p2[i] > 0)
      s += p1[i] * fastlog(p1[i]/p2[i]);
  return s<0?0:s;
}

double jsDist(double *p1, double *p2, int s)
{
  double avg[s];
  for(int i = 0; i < s; i++) 
    avg[i] = .5 * (p1[i] + p2[i]);
  return sqrt(.5 * klDiv(p1, avg, s) + .5 * klDiv(p2, avg, s));
}

double DTopo(LevelInfo *v1, LevelInfo *v2, int s)
{
  double dist = 0;
  for(int i = 0; i < s; i++)
    dist += jsDist(&(v1[i].signature)[0], &(v2[i].signature)[0], v1[i].signature.size());
  return dist/s;
}

//Calculates the average density for the node
/*double avgDensity(vector<LevelInfo> levels)
{
  double total=0;
  for(int i=0; i<levels.size(); i++) total += levels[i].density;
  return total / levels.size();
}*/

//Calculates the D_alphas from node n to all nodes in m2
struct distanceWorker {
  distanceWorker(spectramap::value_type n, spectramap* m2, D_alpha** r, string* o):
    it(m2->begin()), end(m2->end()), result(r), name(n.first), info(n.second), out(o){}
  void operator()()
  {
    ostringstream res;
    int i=0;
    for(; it != end; ++it){
      double topo = DTopo(&info[0], &(it->second)[0], info.size());
      res << name << "\t" << it->first << "\t" << topo << "\n";
      D_alpha* d = new D_alpha(name, it->first, topo);
      result[i++] = d;
    }
    *out = res.str();
  }
  private:
  string name; vector<LevelInfo> info; spectramap::iterator it, end; D_alpha** result; string* out;
};

void applyAlpha(double a, vector<D_alpha>& scores, blastmap* blastscores)
{
  if(blastscores==NULL) a=1;

  if(a==-1){
    vector<double> maxheap;
    for(auto it = scores.begin(); it != scores.end(); it++){
      if(it->get_da() != 0) maxheap.push_back(it->get_da());
      push_heap(maxheap.begin(), maxheap.end());
      if(maxheap.size() > (int)(.0005 * scores.size())){
        pop_heap(maxheap.begin(), maxheap.end());
        maxheap.pop_back();
      }
    }
    double structval = maxheap.front();

    maxheap.clear();
    for(auto it = scores.begin(); it != scores.end(); it++){
      double score = blastscores->at(make_pair(it->get_n1(), it->get_n2()));
      if(score != 0) maxheap.push_back(score);
      push_heap(maxheap.begin(), maxheap.end());
      if(maxheap.size() > (int)(.0005 * blastscores->size())){
        pop_heap(maxheap.begin(), maxheap.end());
        maxheap.pop_back();
      }
    }
    double seqval = maxheap.front();
    a = seqval/structval;
  }

  for(auto it = scores.begin(); it != scores.end(); it++){
    double topo = it->get_da();
    double seq = blastscores==NULL?0:blastscores->at(make_pair(it->get_n1(), it->get_n2()));
    it->update_da(a*topo + (1.0-a)*seq, seq);
  }
}

//Calculates and writes all D_topo values from the given spectral signature files
//Returns 
vector<D_alpha> getDistances(string file1, string file2, string outputname, double a, blastmap* blastscores)
{
  ptime t = bclock::local_time();
  cout << "loading sigs file: " << file1 << "\n";
  spectramap m1 = loadSigs(file1);
  cout << "loaded " << m1.size() << " vertices in " << 
    (bclock::local_time() - t).total_milliseconds() << " milliseconds\n";

  t = bclock::local_time();
  cout << "loading sigs file: " << file2 << "\n";
  spectramap m2 = loadSigs(file2);
  cout << "loaded " << m2.size() << " vertices in " << 
    (bclock::local_time() - t).total_milliseconds() << " milliseconds\n";

  ofstream out(outputname.c_str());

  t = bclock::local_time();

  vector<D_alpha**> results;
  vector<string*> outputs;
  int numthreads = boost::thread::hardware_concurrency();
  if(numthreads < 2) numthreads = 2;
  cout << "calculating D_alphas on " << numthreads << " threads\n";
  boost::threadpool::pool threads(numthreads);

  for(spectramap::iterator it1 = m1.begin(); it1 != m1.end(); ++it1){
    //results.push_back((D_alpha**)calloc(m2.size(), sizeof(D_alpha*)));
    results.push_back(new D_alpha*[m2.size()]);
    outputs.push_back(new string);
    distanceWorker worker(*it1, &m2, results.back(), outputs.back());
    threads.schedule(worker);
  }

  threads.wait();
  cout << "finished calculating " << (m1.size()*m2.size()) << " distances in " << 
    (bclock::local_time()-t).total_milliseconds() << " milliseconds\n";

  vector<D_alpha> allDistances;
  for(int i=0; i<results.size(); i++){
    out << *(outputs[i]);
    delete outputs[i];
    for(int j=0; j<m2.size(); j++){allDistances.push_back(*(results[i][j])); delete results[i][j];}
    delete[] results[i];
  }

  applyAlpha(a, allDistances, blastscores);

  /*ofstream dout((outputname + ".densities").c_str());
  for(spectramap::iterator it1 = m1.begin(); it1 != m1.end(); ++it1)
    dout << it1->first << "\t" << avgDensity(it1->second) << "\n";
  for(spectramap::iterator it2 = m2.begin(); it2 != m2.end(); ++it2) 
    dout << it2->first << "\t" << avgDensity(it2->second) << "\n";
  dout.close();*/
  out.close();
  return allDistances;
}

//int main()
//{
//  getDistances("test.sig.gz", "test.sig.gz", "distances.txt", .1, NULL);
//  return 0;
//}


