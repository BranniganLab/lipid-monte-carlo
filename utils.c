#include <stdlib.h>
#include <math.h>
#include "math.h"
#include "structures.h"
#include "utils.h"

int max(int a, int b)
{
    if (a > b)
        return a;
    else
        return b;
}

int min(int a, int b)
{
    if (a < b)
        return a;
    else
        return b;
}

int in_array(int search_for, int *array, int array_length)
{
    int i;
    for (i = 0; i < array_length; i++)
        if (array[i] == search_for)
            return YES;
    return NO;
}

double ran_theta(long *base)
{
    double x, y, theta; 
    x = ran3(base) - 0.5;
    y = ran3(base) - 0.5;
    theta = atan(y/x);
    return theta; 
}