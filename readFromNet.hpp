#include <string>
#include <iostream>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "graph.hpp"

using std::string;
using std::ifstream;
using std::cout;

typedef boost::posix_time::microsec_clock bclock;
typedef boost::posix_time::ptime ptime;

Graph readFromNet(string fileName)
{ 
//  ptime t = bclock::local_time();
  Graph result;
  ifstream fin(fileName);
  if(!fin.good()) 
  {
    cout << "error loading file: " << fileName << "\n";
    exit(0);
  }
  unsigned start = fileName.find_last_of("/");
  unsigned end = fileName.find(".net");
  if(start == string::npos) start = -1;
  string name = fileName.substr(start+1,end-start-1);
  result.setName(name);

  vlist V;
  while(1)
  {
    string n1,n2;
    fin >> n1 >> n2;
    if(fin.eof()) break;
    if(V.find(n1) == V.end())
      {result.addVertex(n1); V.insert(n1);}
    if(V.find(n2) == V.end())
      {result.addVertex(n2); V.insert(n2);}
    result.addEdge(n1,n2);
  }
  fin.close();
//  cout << "extracted: " << result.getName() << ".net in " <<
//    (bclock::local_time() - t).total_milliseconds() << " milliseconds\n";
  return result;
}

