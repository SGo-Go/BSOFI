#include <windows.h>
LARGE_INTEGER CLOCK_TIC()
{
LARGE_INTEGER  time ;
    QueryPerformanceCounter(&time);
return (time);
}
