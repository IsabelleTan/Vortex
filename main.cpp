//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "ADIdiffusion_test.h"
#include <unistd.h>
#define GetCurrentDir getcwd


using namespace std;

void printWorkingDirectory(){
    char cCurrentPath[FILENAME_MAX];

    if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
    {
        cout << "Something went wrong when trying to obtain the current working directory" << endl;
        return;
    }

    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0';

    printf ("The current working directory is %s \n", cCurrentPath);
    return;
}

int main(){

    ADI_test_output();

    return 0;
}


/* int main(){
    value_type * a = new value_type[16*16];
    value_type * x = new value_type[16];
    value_type * b = new value_type[16]{4,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    a[0] = 5;

    for (int i = 0; i < 16; ++i) {
        if(i%4 == 0){
            a[i*16 + i] = 5;
            a[i*16 + i-1] = 0;
            a[(i-1)*16 + i] = 0;
        } else {
            a[i*16 + i] = 5;
            a[i*16 + i-1] = -2;
            a[(i-1)*16 + i] = -2;
        }
    }

    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            cout << a[i*16+j] << " ";
        }
        cout << "\n";

    }

    ThomasAlg(16, -2, 5, -2, x, b);

    for (int k = 0; k < 16; ++k) {
        cout << x[k] << ", ";
    }

    delete[] a;
    return 0;
} */