unsigned long long walltime( )
{
#ifdef __PPC970__
  unsigned int         HightB, HightA, Low;
  unsigned long long   tt;

  do
  {
    asm volatile("mftbu %0" : "=r"(HightB) );
    asm volatile("mftb  %0" : "=r"(Low) );
    asm volatile("mftbu %0" : "=r"(HightA) );
  }
  while (HightB != HightA);

  tt = ((unsigned long long)HightA<<32) | ((unsigned long long)Low);
  return tt;
#endif
}
