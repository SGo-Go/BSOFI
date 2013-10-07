const inline unsigned long long int clock_tic_()
{
  unsigned long long int x;
  /* __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));*/
  __asm__ volatile ("rdtsc\n" : "=A" (x));
  return x;
}
const inline unsigned long long int clock_tic__()
{
  unsigned long long int x;
  /* __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));*/
  __asm__ volatile ("rdtsc\n" : "=A" (x));
  return x;
}
const inline unsigned long long int _clock_tic()
{
  unsigned long long int x;
  /* __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));*/
  __asm__ volatile ("rdtsc\n" : "=A" (x));
  return x;
}
