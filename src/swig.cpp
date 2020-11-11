#include "titanlib.h"
#include <iostream>

using namespace titanlib;

float* titanlib::test_array(float* v, int n) {
    int count = 0;
    for(int i = 0; i < n; i++)
        count++;
    return v;
 }
void titanlib::test_not_implemented_exception() {
    throw titanlib::not_implemented_exception();
}
