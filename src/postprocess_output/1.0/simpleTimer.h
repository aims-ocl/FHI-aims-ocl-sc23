#ifndef _SIMPLE_TIMER_H_
#define _SIMPLE_TIMER_H_

#include <vtkSmartPointer.h>
#include <vtkTimerLog.h>

class simpleTimer
{
  private:
    vtkSmartPointer<vtkTimerLog> timer;
    int mpiRank;
    char myLabel[128];

  public:
    simpleTimer();
    simpleTimer(const char * label);
    simpleTimer(const char * label, int _mpiRank);
    ~simpleTimer();
    void setVerbose();
};

#endif // _SIMPLE_TIMER_H_
