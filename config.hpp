#include <stdio.h>
#include <stdlib.h>
#include <string>

using std::string;

struct ConfigData
{
  string Ggexf="", Hgexf="",
         Gsigs="", Hsigs="",
         SeqScores="";
  int numProcessors=20,
      hops=4,
      nneighbors=-1, // dummy for all
      searchiter = 10; // not sure
  double alpha=-1.0,
         beta=100,
         ratio=8.0; // not sure
  bool dumpDistances=false;
  void use(string s);
  void configure(string filename);
  void print();
};

void ConfigData::use(string s)
{
  string pos[] = {"network1: ", "network2: ", "sigs1: ", "sigs2: ",
                  "sequencescores: ", "nneighbors: ", "searchiter: ",
                  "hops: ", "alpha: ", "beta: ", "ratio: "};
  for(int i=0;i<9;i++)
    if(s.size() > pos[i].size() && s.substr(0,pos[i].size()) == pos[i])
    {
      s = s.substr(pos[i].size());
      switch (i)
      {
        case 0: Ggexf = s; break;
        case 1: Hgexf = s; break;
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
        default: printf("\"%s\" may be incorrect",s.c_str()); break;
      }
    }
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
  fscanf(fin,"%[^\n]\n",next);
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
  printf("Ggexf: %s\nHgexf: %s\nGsigs: %s\nHsigs: %s\nSeqScores: %s\n"
          "nneighbors: %d\nsearchiter: %d\nbeta: %lf\nratio: %lf\n", 
          Ggexf.c_str(), 
          Hgexf.c_str(), 
          Gsigs.c_str(), 
          Hsigs.c_str(), 
          SeqScores.c_str(), 
          nneighbors,
          searchiter, 
          beta, 
          ratio);
}

