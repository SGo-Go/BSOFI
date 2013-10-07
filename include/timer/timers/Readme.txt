integer(kind=8)    :: time_start, time_end, time
real(kind=8)      :: time_in_seconds,frequency

time_start=clock_tick( )
code to be measured
time_end=clock_tic()

time=time_end-time_start

! time are the number of cycles
! to get the time in  seconds, you have to divide time/frequency

time_in_seconds=time/frequency

