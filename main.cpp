//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "ADIdiffusion_test.h"
#include "fields.h"
#include <unistd.h>
#define GetCurrentDir getcwd

using namespace Eigen;
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

void testFunction(Ref<MatrixXd> m){
    cout << "This is a test function" << endl;
}

int main(){

    ADI_test_output();
}
