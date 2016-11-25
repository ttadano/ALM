/*
 memory.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <cstdlib>
#include <new>
#include "memory.h"

// memsize calculator

unsigned long memsize_in_MB(const int size_of_one,
                            const unsigned int n1)
{
    unsigned long n = n1 * size_of_one;
    return n / 1000000;
}

unsigned long memsize_in_MB(const int size_of_one,
                            const unsigned int n1,
                            const unsigned int n2)
{
    unsigned long n = n1 * n2 * size_of_one;
    return n / 1000000;
}

unsigned long memsize_in_MB(const int size_of_one,
                            const unsigned int n1,
                            const unsigned int n2,
                            const unsigned int n3)
{
    unsigned long n = n1 * n2 * n3 * size_of_one;
    return n / 1000000;
}

unsigned long memsize_in_MB(const int size_of_one,
                            const unsigned int n1,
                            const unsigned int n2,
                            const unsigned int n3,
                            const unsigned int n4)
{
    unsigned long n = n1 * n2 * n3 * n4 * size_of_one;
    return n / 1000000;
}
