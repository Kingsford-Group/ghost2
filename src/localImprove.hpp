#include <vector>
#include <string>
#include <tgmath.h>
#include <stdlib.h>

#include "graph.hpp"
#include "spectralToDistance.hpp"

using std::string;
using std::vector;
using std::pair;

typedef boost::unordered_map<pair<string,string>, D_alpha> distmap;

struct move{
  string u, w, v, up;
  double s0, s1, s2;
  move(string u, string w, string up, string v) : u(u),w(w),up(up),v(v){} 
  void apply(bmap *f)
  {
    f->left.erase(u);
    if(up != "") f->left.erase(up);
    f->right.erase(w);
    f->right.erase(v);
    f->insert(bmap::value_type(u, v));
    if(up != "") f->insert(bmap::value_type(up, w));
  }
  void calcScores(Graph& G, Graph& H, blastmap *dists, bmap *f)
  {
    auto neighborsu = G.neighbors(u);
    auto neighborsw = H.neighbors(w);
    auto neighborsv = H.neighbors(v);
    s0 = s1 = s2 = 0;
    //Go thru all the edges that u is a part of, to see what changes
    for(auto it = neighborsu.begin(); it != neighborsu.end(); it++){
      auto mapped = f->left.find(*it);
      if(mapped != f->left.end()){
        //There is a matching edge to w that we are losing
        if(neighborsw.find(mapped->second) != neighborsw.end()) s0--;
        //There is a matching edge to v that we are adding
        if(neighborsv.find(mapped->second) != neighborsv.end()) s0++;
      }
    }
    if(up != ""){
      auto neighborsup = G.neighbors(up);
      //Go thru all the edges that up is a part of, to see what changes
      for(auto it = neighborsup.begin(); it != neighborsup.end(); it++){
        auto mapped = f->left.find(*it);
        if(mapped != f->left.end()){
          //There is a matching edge to w that we are adding
          if(neighborsw.find(mapped->second) != neighborsw.end()) s0++;
          //There is a matching edge to v that we are losing
          if(neighborsv.find(mapped->second) != neighborsv.end()) s0--;
        }
      }
    }

    if(s0 > 0 && dists){
      auto uw = dists->find(make_pair(u,w));
      auto uv = dists->find(make_pair(u,v));
      s1 = (uw != dists->end() ? uw->second : 1) - (uv != dists->end() ? uv->second : 1);

      if(up != ""){
        uw = dists->find(make_pair(up,w));
        uv = dists->find(make_pair(up,v));
        s2 = (uv != dists->end() ? uv->second : 1) - (uw != dists->end() ? uw->second : 1);
      }
    }
  }
};

int searchIter(Graph& G, Graph& H, blastmap *dists, bmap *f, double ratio)
{
  int totalDelt = 0;
  vector<string> hverts = H.nodes();
  vector<pair<string,string>> pairs;
  for(auto it = f->left.begin(); it != f->left.end(); it++) 
    pairs.push_back(make_pair(it->first, it->second));
  for(auto it = pairs.begin(); it != pairs.end(); it++){
    string u = it->first;
    string w = it->second;
    
    move *bestMove = NULL;
    bool constrained = rand() > (ratio*RAND_MAX);
    for(auto it2 = hverts.begin(); it2 != hverts.end(); it2++){
      string v = *it2;
      string up = "";
      auto upt = f->right.find(v);
      if(upt != f->right.end()) up = upt->second;
      move m(u, w, up, v);
      m.calcScores(G, H, dists, f);
      //Check if move is acceptible
      if(m.s0 > 0 && m.s1 >= 0 && (!constrained || m.s2 >= 0)){
        if(!bestMove) bestMove = new move(m);
        else if(m.s0 > bestMove->s0) *bestMove = m;
        else if(m.s0 == bestMove->s0 && m.s1 > bestMove->s1) *bestMove = m;
        else if(m.s0 == bestMove->s0 && m.s1 == bestMove->s1 && m.s2 > bestMove->s2) *bestMove = m;
      }
    }

    if(bestMove){
      bestMove->apply(f);
      totalDelt += bestMove->s0;
      delete bestMove;
    }
  }
  return totalDelt; 
}

int localImprove(Graph& G, Graph& H, blastmap *dists, bmap *f, int iters, double ratio)
{
  vector<double> ratios;
  double sum;
  for(int i=0; i < iters; i++){
    double d = exp(-i);
    sum += d;
    ratios.push_back(d);
  }
  sum = ratio/sum;
  int totalDelt = 0;
  for(int i=0; i < iters; i++){
    ptime t = bclock::local_time();
    cout << "running PISWAP iteration, " << (ratios[i] * sum * 100) << "\% unconstrained\n";
    int d = searchIter(G, H, dists, f, ratios[i] * sum);
    cout << "added " << d << " edges in " << (bclock::local_time() - t).total_milliseconds() << "ms\n";
    totalDelt += d;
    if(d == 0) break;
  }
  return totalDelt;
}
