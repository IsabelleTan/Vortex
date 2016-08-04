//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "fields_test.h"
#include "datapoints.h"
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

    cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; /* not really required */

    printf ("The current working directory is %s", cCurrentPath);
    return;
}

int main(){
    int N = 100;
    value_type2 * x = new value_type2[N];

    // Make the first frame
    for (int i = 0; i < N; ++i) {
        x[i] = (i%10 - 5)*(i%10 - 5);
    }

    write_to_file("test_0", N, x);

    // Make the 2nd frame
    for (int i = 0; i < N; ++i) {
        x[i] += 5;
    }

    write_to_file("test_1", N, x);

    // Make the 3nd frame
    for (int i = 0; i < N; ++i) {
        x[i] += 5;
    }

    write_to_file("test_2", N, x);

    // Make the 4nd frame
    for (int i = 0; i < N; ++i) {
        x[i] += 5;
    }

    write_to_file("test_3", N, x);

    // Make the 5nd frame
    for (int i = 0; i < N; ++i) {
        x[i] += 5;
    }

    write_to_file("test_4", N, x);

    printWorkingDirectory();




    delete[] x;



    return 0;
}

