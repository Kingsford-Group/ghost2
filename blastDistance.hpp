#include <string>
#include <boost/unordered_map.hpp>
#include <utility>
#include <iostream>
#include <fstream>

using std::string;
using std::pair;
using std::ifstream;

typedef boost::unordered_map<pair<string,string>, double> umap;

int getBlastMap(string filename)
{
  umap bevals;
  ifstream fin (filename);
  while(1)
  {
    string n1,n2;
    double d;
    fin >> n1 >> n2 >> d;
    if(fin.eof()) break;
    bevals[make_pair(n1,n2)] = d;
  }
  fin.close();
  umap::iterator iter = bevals.begin(),
  iend = bevals.end();
  return bevals;
//  for(; iter != iend; ++iter)
//  {
//    std::cout << (iter->first).first << "\t"
//              << (iter->first).second << "\t"
//              << iter->second << "\n";
//  }
}
