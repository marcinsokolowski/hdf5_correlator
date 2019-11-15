#include <stdio.h>
#include <stdlib.h>

#include <complex>

int main()
{
   std::complex<double> z = 3;
   
   printf("mag = %.1f , abs = %.1f, arg = %.1f\n",std::norm(z),std::abs(z),std::arg(z));

}