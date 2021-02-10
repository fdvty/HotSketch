#ifndef _PARAM_H_
#define _PARAM_H_

#include <iostream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <fstream>
#include <set>
#include <map>
#include <cstdlib>

using namespace std;


#define __max(a,b)    (((a) > (b)) ? (a) : (b))
#define __min(a,b)    (((a) < (b)) ? (a) : (b))


#define CONSTANT_NUMBER 2654435761u
#define CalculateBucketPos(fp) (((fp) * CONSTANT_NUMBER) >> 15)

#define CELL_PER_BUCKET 16


struct Bucket{
public:
    uint32_t keys[CELL_PER_BUCKET];
    uint32_t freq[CELL_PER_BUCKET];
    double   timestamps[CELL_PER_BUCKET];
};


#endif // _PARAM_H_
