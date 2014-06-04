#include <string>
using std::string;
class D_alpha
{
  string n1,n2;
  double d;
  double s;
  public:
    D_alpha(string a,string b,double c) {n1=a; n2=b; d=c; s=0;};
    string get_n1() {return n1;};
    string get_n2() {return n2;};
    double get_da() {return d;};
    double get_ds() {return s;};
    void update_da(double n, double m) {d=n;s=m;};
};
class CompareD_alphaG
{
  public:
    bool operator()(D_alpha& x, D_alpha& y) 
   {
     return x.get_da() > y.get_da();
   }
};
class CompareD_alphaL
{
  public:
    bool operator()(D_alpha& x, D_alpha& y) 
   {
     return x.get_da() < y.get_da();
   }
};
