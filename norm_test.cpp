#include <stdio.h>
#include <complex>

/* RESULTS :
   a = 2.00 / 2.00 -> a.norm = 8.00
   abs = 2.83
*/

int main()
{
   std::complex<float> a(2,2);
    
   printf("a = %.2f / %.2f -> a.norm = %.2f\n",a.real(),a.imag(),std::norm(a));
   printf("abs = %.2f\n",std::abs(a));
}