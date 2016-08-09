//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "ADIdiffusion_test.h"
#include <unistd.h>
#define GetCurrentDir getcwd


using namespace std;

/*
 * A function to print the location of the current working directory (useful to find the output files)
 */
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
    printWorkingDirectory();
    ADI_test_output();
}
