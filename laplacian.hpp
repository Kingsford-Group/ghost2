#include <iostream>
#include <vector>
#include <math.h>
#include <boost/unordered_map.hpp>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "graph.h"

using std::vector;
using namespace Eigen;
typedef boost::unordered_map<string, int> degMap;
typedef boost::unordered_map<string, vlist> neighborhood;

class AdjacencyMatrix
{
  MatrixXd m;
  vlist nodes;
  degMap dmap;
  neighborhood neighborhoods;
  vector<string> prev;
  int nvert=0, nedge=0;
  public:
    void extend(Graph *input, string source);
    vlist getNodes() {return nodes;};
    vector<string> getPrev() {return prev;};
    vector<double> getEigen();
    double getDensity() {return (double)(nvert) / nedge;};
};

void AdjacencyMatrix::extend(Graph *input, string source)
{
  if(nodes.empty())
  {
    nodes.insert(source);
    prev.push_back(source);
    neighborhoods[source] = (*input).neighbors(source);
    dmap[source] = 0;
  }

  vector<string> next;
  for(int i=0;i<prev.size();i++)
  {
    vlist::iterator vit = neighborhoods[prev[i]].begin(), 
                    vend = neighborhoods[prev[i]].end();
    for(; vit != vend; vit++)
      if(nodes.find(*vit) == nodes.end())
      {
        nodes.insert(*vit);
        next.push_back(*vit);
        neighborhoods[*vit] = (*input).neighbors(*vit);
        dmap[*vit] = 0;
        vlist neighbors = (*input).neighbors(*vit);
        vlist::iterator it = neighbors.begin(), iend = neighbors.end();
        for(; it != iend; it++)
          if(nodes.find(*it) != nodes.end())
          {
            dmap[*vit] += 1;
            dmap[*it] += 1;
            nedge++;
          }
      }
  }
  prev = next;
  nvert = nodes.size();
  m.resize(nvert,nvert);

  if(prev.size() != 0)
  {
    vlist::iterator it1 = nodes.begin(), iend = nodes.end();
    int c1=0, c2=0;
    for(; it1 != iend; it1++)
    {
      vlist::iterator it2 = nodes.begin();
      for(; it2 != iend; it2++)
      {
        if(c1 == c2)
          m(c1,c2) = 1;
        else
        {
          if(neighborhoods[*it1].find(*it2) != neighborhoods[*it1].end())
            m(c1,c2) = -1. / sqrt(dmap[*it1] * dmap[*it2]);
          else
            m(c1,c2) = 0;
        }
        c2++;
      } 
      c1++; 
      c2=0;
    }
  }
}

vector<double> AdjacencyMatrix::getEigen()
{
  vector<double> v;
  if(nodes.size() == 0)
    return v;
  VectorXd evals = m.selfadjointView<Eigen::Upper>().eigenvalues();
  for(int i=0; i<evals.rows(); i++)
    v.push_back(evals(i));
  return v;
}
