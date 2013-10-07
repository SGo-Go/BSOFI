const inline unsigned long long clock_tic()
{
   unsigned int a,d;
   asm volatile("rdtsc" : "=a" (a), "=d" (d));
   return ((unsigned long long)a) | (((unsigned long long)d)<<32);
}
const inline unsigned long long clock_tic_()
{
   unsigned int a,d;
   asm volatile("rdtsc" : "=a" (a), "=d" (d));
   return ((unsigned long long)a) | (((unsigned long long)d)<<32);
}
const inline unsigned long long clock_tic__()
{
   unsigned int a,d;
   asm volatile("rdtsc" : "=a" (a), "=d" (d));
   return ((unsigned long long)a) | (((unsigned long long)d)<<32);
}
