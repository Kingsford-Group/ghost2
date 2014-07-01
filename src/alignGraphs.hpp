#include <vector>
#include <algorithm>
#include <string>
#include <boost/bimap/bimap.hpp>
#include <boost/bimap/unordered_set_of.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <iostream>
#include <fstream>
#include "graph.hpp"
#include "alignmentExtender.hpp"

using std::string;
using std::ifstream;
using std::ofstream;
using std::vector;
using std::cout;

typedef boost::posix_time::microsec_clock bclock;
typedef boost::posix_time::ptime ptime;

/* helper for seedAndExtend and GreedyQAPExtend */
bool isNotUsed(string n1, string n2, bmap *f)
{
  bool flag1 = false;
  bool flag2 = false;
  try{(*f).left.at(n1);}
  catch(std::out_of_range & e) {flag1=true;}
  try{(*f).right.at(n2);}
  catch(std::out_of_range & e) {flag2=true;}
  return (flag1 && flag2);
}

distmap getDistMap(vector<D_alpha>& v)
{
  distmap m;
  for(auto it = v.begin(); it != v.end(); it++) m.insert(make_pair(make_pair(it->get_n1(), it->get_n2()), *it));
  return m;
}

/* initialize bmap f
 * while minP not empty, 
 *  pop D_alpha d
 *  if no part of d in f, GreedyQAPEXtend(G,H,d,f) */
bmap seedAndExtend(Graph& G, Graph& H, vector<D_alpha>& minP, int k)
{
  vector<D_alpha> allDists(minP);
  distmap dmap = getDistMap(allDists);
  bmap f;
  while(minP.size()>0)
  {
    D_alpha d = minP.front();
    pop_heap (minP.begin(), minP.end(),CompareD_alphaG());
    minP.pop_back();

    string n1 = d.get_n1();
    string n2 = d.get_n2();
    if(isNotUsed(n1,n2,&f))
      extendAlignment(d, G, H, &f, allDists, dmap, k);
  }
  return f;
}

/* prints all key-value pairs in f */
void printMap(bmap f, string gname, string hname)
{
  ofstream fout (gname+"_vs_"+hname+".af");
  bmap::left_const_iterator iter = f.left.begin(),
  iend = f.left.end();
  for(; iter != iend; ++iter)
  {
    fout << iter->first << "\t" << iter->second << "\n";
  }
  fout.close();
}

void toDotFile(Graph& G, Graph& H, bmap f)
{
  boost::unordered_set<string> verts;
  vector<pair<string,string>> edges1;
  vector<pair<string,string>> edges2;
  vector<pair<string,string>> edges3;
  for(auto it = f.left.begin(); it != f.left.end(); it++){
    auto adj1 = G.neighbors(it->first);
    auto adj2 = H.neighbors(it->second);
    for(auto it2 = adj1.begin(); it2 != adj1.end(); it2++){
      auto f_a = f.left.find(*it2);
      if(f_a != f.left.end()){
        if(adj2.find(f_a->second) != adj2.end()){
          if(find(edges1.begin(), edges1.end(), make_pair(f_a->first, it->first)) == edges1.end()){
            edges1.push_back(make_pair(it->first, f_a->first));
            verts.insert(it->first);
            verts.insert(f_a->first);
          }
        }else{
          if(find(edges2.begin(), edges2.end(), make_pair(f_a->first, it->first)) == edges2.end()){
            //edges2.push_back(make_pair(it->first, f_a->first));
            //verts.insert(it->first);
            //verts.insert(f_a->first);
          }
        }
      }
    }
  }
  ofstream fout ("test.dot");
  fout << "graph G {\n";
  for(auto it = verts.begin(); it != verts.end(); it++) fout << *it << "[label = \"" << *it << "\\n" << f.left.at(*it) << "\"]\n";
  for(auto it = edges1.begin(); it != edges1.end(); it++) fout << it->first << " -- " << it->second << "\n";
  for(auto it = edges2.begin(); it != edges2.end(); it++) fout << it->first << " -- " << it->second << " [color = \"red\"];\n";
  fout << "}";
  fout.close();
}

void printICS(Graph& G, Graph& H, bmap& result)
{
  int matchingEdges = 0;
  int edgesH = 0;
  for(auto it = result.left.begin(); it != result.left.end(); it++)
  {
    auto adj1 = G.neighbors(it->first);
    auto adj2 = H.neighbors(it->second);
    for(auto it2 = adj1.begin(); it2 != adj1.end(); it2++){
      auto f_a = result.left.find(*it2);
      if(f_a != result.left.end() && adj2.find(f_a->second) != adj2.end()){
        matchingEdges++;
        if(*it2 == it->first) matchingEdges++; //Double count self loops
      }
    }
    for(auto it2 = adj2.begin(); it2 != adj2.end(); it2++) 
      if(result.right.find(*it2) != result.right.end()){
        edgesH++;
        if(*it2 == it->second) edgesH++; //Double count self loops
      }
  }
  vector<string> nodesG = G.nodes();
  int edgesG = 0;
  for(auto it = nodesG.begin(); it != nodesG.end(); it++){
    auto n = G.neighbors(*it);
    edgesG += n.size();
    if(n.find(*it) != n.end()) edgesG++; //Double count self loops
  }

  double ec = ((double)matchingEdges / edgesG) * 100.0;
  double ics = ((double)matchingEdges / edgesH) * 100.0;
  cout << "Edge correctness " << matchingEdges/2 << " / " << edgesG/2 << " = " << ec << "\%\n";
  cout << "ICS = " << ics << "\%\n";
}

bmap alignGraphs(Graph& G, Graph& H, vector<D_alpha>& distances, int k)
{
  cout << "aligning graphs...\n";
  vector<D_alpha> minP;  // empty vector is a heap
  for(auto it = distances.begin(); it != distances.end(); it++)
  {
    minP.push_back(*it);
    push_heap(minP.begin(),minP.end(),CompareD_alphaG());
  }

  ptime t = bclock::local_time();
  bmap result = seedAndExtend(G, H, minP, k);
  cout << "aligned graphs in " << (bclock::local_time()-t).total_milliseconds() << " milliseconds\n";

  printICS(G, H, result);
  return result;
  //toDotFile(G, H, result);
}
