#include <string>
#include <boost/unordered_map.hpp>
#include <utility>
#include <iostream>
#include <fstream>

using std::string;
using std::pair;
using std::ifstream;

typedef boost::unordered_map<pair<string,string>, double> umap;

int main()
{
  umap bevals;
  ifstream fin ("AThaliana_vs_DMel.evalues");
//  ifstream fin ("test.evalues");
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
  for(; iter != iend; ++iter)
  {
//    std::cout << iter->first << "\t" << iter->second << "\n";
    std::cout << (iter->first).first << "\t"
              << (iter->first).second << "\t"
              << iter->second << "\n";
  }
  return 0;
}
