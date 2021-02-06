#include "param.h"
#include "trace.h"
#include "HotSketch.h"
#include <cmath>

vector<pair<uint32_t, double> > vec;
map<uint32_t, uint32_t> mp;


void test_hotsketch(int hotsize, int hot_threshold, int time_threshold){
    HotSketch hotsketch(hotsize, time_threshold);
    
    timespec time1, time2;
    uint64_t resns = 0;
    clock_gettime(CLOCK_MONOTONIC, &time1);
    int key_num = vec.size();
    for(int i = 0; i < key_num; ++i)
        hotsketch.insert(vec[i].first, vec[i].second);
    clock_gettime(CLOCK_MONOTONIC, &time2);
    resns += (long long)(time2.tv_sec - time1.tv_sec) * 1000000000LL + (time2.tv_nsec - time1.tv_nsec);
    double insertionMops = (double)1000.0 * (key_num) / resns;
    printf("# Key: %d, Insertion Mops: %lf\n", key_num, insertionMops); 

    uint32_t value = 0, all = 0, hit = 0, size = 0;
    double aae = 0, are = 0, cr = 0, pr = 0;
    int cnt = 0;
    for(auto it = mp.begin(); it != mp.end(); ++it){
        value = hotsketch.query(it->first, vec[key_num-1].second);
        if(it->second > hot_threshold){
            all++;
            if(value > hot_threshold){
                hit++;
                aae += fabs(it->second - value);
                are += fabs(it->second - value) / (double) it->second;
            }
        }
        if(value > hot_threshold)
            size++;
        cnt++;

        // if((value > hot_threshold) && (it->second <= hot_threshold)){
        //     printf("%u - %u\n", value, it->second);
        // }
    }

    aae /= hit;
    are /= hit;
    cr = hit / (double)all;
    pr = hit / (double)size;

    printf("hotsketch cells: %d\n", hotsize);
    printf("AAE:%lf, ARE:%lf, CR: %lf, PR: %lf\n", aae, are, cr, pr);
    printf("cnt:%d, size:%d, all:%d, hit: %d\n", cnt, size, all, hit);
}




int main(){
    // loadCAIDA18(vec, mp, "/Users/liuzirui/caida/130000.dat", 1, 10000000);
    loadCAIDA18(vec, mp, "/Users/liuzirui/caida/130000.dat", 1, 1000, 6000000);
    test_hotsketch(2000, 1000, 1);
    return 0;
}