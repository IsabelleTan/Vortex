//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "fields_test.h"
#include "datapoints.h"


using namespace std;

int main(){
    int N;
    value_type2 * x;
    value_type2 * y;

    cout << grid_test() << endl;

    // Print the content of x and y
    for (int i = 0; i < 10; ++i) {
        cout << "x[i] = " << x[i] << endl;
    }

    return 0;
}

