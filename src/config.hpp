#include <stdio.h>
#include <stdlib.h>
#include <string>

using std::string;

struct ConfigData
{
  string Ggraph="", Hgraph="",
         Gsigs="", Hsigs="",
         SeqScores="", AlignFile="", DistFile="";
  int numProcessors=-1,
      hops=4,
      nneighbors=-1, // dummy for all
      searchiter = 10; // not sure
  double alpha=-1.0,
         beta=1.0,
         ratio=8.0; // not sure
  bool dumpSignatures=false,
       dumpDistances=false,
       directed=false;
  void use(string s);
  void configure(string filename);
  void print();
};

void ConfigData::use(string s)
{
  if(s.substr(0, 2) == "//") return;
  string pos[] = {"network1: ", "network2: ", "sigs1: ", "sigs2: ",
                  "sequencescores: ", "nneighbors: ", "searchiter: ",
                  "hops: ", "processors: ", "alpha: ", "beta: ", "ratio: ", 
                  "dumpSignatures: ", "dumpDistances: ", "directed: ",
                  "alignFile: ", "distFile: "};
  bool used=false;
  for(int i=0;i<17;i++)
    if(s.size() > pos[i].size() && s.substr(0,pos[i].size()) == pos[i])
    {
      s = s.substr(pos[i].size());
      used=true;
      switch (i)
      {
        case 0: Ggraph = s; break;
        case 1: Hgraph = s; break;
        case 2: Gsigs = s; break;
        case 3: Hsigs = s; break;
        case 4: SeqScores = s; break;
        case 5: if(s=="all") nneighbors = -1;
                else nneighbors = atoi(s.c_str()); break;
        case 6: searchiter = atoi(s.c_str()); break;
        case 7: hops = atoi(s.c_str()); break;
        case 8: numProcessors = atoi(s.c_str()); break;
        case 9: alpha = atof(s.c_str()); break;
        case 10: beta = atof(s.c_str()); break;
        case 11: ratio = atof(s.c_str()); break;
        case 12: if(s=="true") dumpSignatures=true; break;
        case 13: if(s=="true") dumpDistances=true; break;
        case 14: if(s=="true") directed=true; break;
        case 15: AlignFile = s; break;
        case 16: DistFile = s; break;
        default: break;
      }
    }
  if(!used)
    printf("WARNING: the configuration option \"%s\" may be incorrect\n\n", 
           s.c_str());
}

void ConfigData::configure(string filename)
{
  // try to open
  FILE *fin = fopen(filename.c_str(),"r");
  if(fin == NULL) 
  {
    printf("config file not found\n"); 
    exit(0);
  }
  // check for [main]
  char next[200];
  if(fscanf(fin,"%[^\n]\n",next) == 0)
  {
    printf("config file appears empty\n");
    exit(0);
  }
  if(string(next) != "[main]")
  {
    printf("no [main] in config file\n"); 
    exit(0);
  }
  // update for the rest of the lines
  while(fscanf(fin,"%[^\n]\n",next) != EOF)
    use(string(next));
  fclose(fin);
}

void ConfigData::print()
{
  printf("Ggraph: %s\nHgraph: %s\nGsigs: %s\nHsigs: %s\nSeqScores: %s\n"
          "nneighbors: %d\nsearchiter: %d\nbeta: %lf\nratio: %lf\n", 
          Ggraph.c_str(), 
          Hgraph.c_str(), 
          Gsigs.c_str(), 
          Hsigs.c_str(), 
          SeqScores.c_str(), 
          nneighbors,
          searchiter, 
          beta, 
          ratio);
}

