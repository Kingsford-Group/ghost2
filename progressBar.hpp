#include <boost/thread.hpp>
#include <iostream>

class ProgressBar
{
  int prev=0, count=0, max;
  double progress=0;
  int barWidth=50;
  boost::mutex mutex;
  public:
    ProgressBar(int i) {max=i;};
    void update();
};

void ProgressBar::update()
{
  count++;
  progress = (double)(count) / max;
  if(2*barWidth*progress > prev)
    if(mutex.try_lock())
    {
      prev++;
      std::cout << "[";
      int pos = barWidth * progress;
      for(int i=0;i<barWidth;i++)
      {
        if(i<pos) std::cout << "=";
        else if(i==pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "]" << int(progress*100.0) << "%\r";
      std::cout.flush();
      mutex.unlock();
    }
}
