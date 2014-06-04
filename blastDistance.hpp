#include <string>
#include <boost/unordered_map.hpp>
#include <utility>
#include <iostream>
#include <fstream>

using std::string;
using std::pair;
using std::ifstream;

typedef boost::unordered_map<pair<string,string>, double> blastMap;

blastMap getBlastMap(string filename)
{
  blastMap bevals;
  ifstream fin (filename);
  std::cout << filename << " is the filename\n";
  while(1)
  {
    string n1,n2;
    double d;
    fin >> n1 >> n2 >> d;
    if(fin.eof()) break;
    bevals[make_pair(n1,n2)] = d;
  }
  fin.close();
  return bevals;
}
