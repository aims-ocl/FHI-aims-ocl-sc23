#include <string.h>
#include <iostream>
#include <simpleTimer.h>

using namespace std;

simpleTimer::simpleTimer()
{
  simpleTimer("default");
}

simpleTimer::simpleTimer(const char * label)
{
  mpiRank=-1;
  //
  strncpy( myLabel, label, 127 );
  timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();
}

simpleTimer::simpleTimer(const char * label, int _mpiRank)
{
  mpiRank=_mpiRank;
  strncpy( myLabel, label, 127 );
  timer = vtkSmartPointer<vtkTimerLog>::New();
  timer->StartTimer();
}

simpleTimer::~simpleTimer()
{
  timer->StopTimer();
  if (mpiRank==0)
    cout << "Timer: " << myLabel << ":  " << timer->GetElapsedTime() << " s" << endl;
}

void simpleTimer::setVerbose()
{
  mpiRank=0;
}
