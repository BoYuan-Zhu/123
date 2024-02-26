#include <iostream>
using namespace std;
int main()
{
    int F[31];
    F[0]=0;
    F[1]=1;
    int i;
    for(i=2;i<31;i++)
      F[i]=F[i-1]+F[i-2];
    cout << "F7 = " << F[7] << ";\n";
    cout << "F15 = "<< F[15] << ";\n";
    cout << "F31 = " << F[31] << ";\n";
    return 0;
}