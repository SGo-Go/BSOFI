#include <windows.h>
//#include <Winbase.h>
long long CLOCK_TIC()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long CLOCK_FREQUENCY()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
long long CLOCK_TIC_()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long CLOCK_FREQUENCY_()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
long long CLOCK_TIC__()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long CLOCK_FREQUENCY__()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
long long clock_tic()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long clock_frequency()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
long long clock_tic_()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long clock_frequency_()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
long long clock_tic__()
{
LARGE_INTEGER  time ;
long long yy;
    QueryPerformanceCounter(&time);
yy=time.QuadPart;
return (yy);
}
long long clock_frequency__()
{
LARGE_INTEGER  frequency ;
long long yy;
    QueryPerformanceFrequency(&frequency);
yy=frequency.QuadPart;
return (yy);
}
