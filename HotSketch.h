#ifndef _HOTSKETCH_H_
#define _HOTSKETCH_H_

#include "param.h"
#include "murmur3.h"

#include "BOBHash32.h"

class HotSketch{
public:
    Bucket* buckets;
    double time_threshold;
    uint32_t bucket_num;
    uint32_t seed;
    // BOBHash32 bhash; 

    HotSketch(uint32_t cell_num_, double time_threshold_){
        bucket_num = (cell_num_ + CELL_PER_BUCKET - 1) / CELL_PER_BUCKET;
        time_threshold = time_threshold_;
        seed = rand()%MAX_PRIME32; 
        buckets = new((std::align_val_t)64) Bucket[bucket_num];
        memset(buckets, 0, bucket_num*sizeof(Bucket));

        // bhash.initialize(seed);
    }

    ~HotSketch(){
        delete [] buckets;
    }

    uint32_t query(uint32_t key, double nowtime){
        int pos = CalculatePos(key); 
        for(int i = 0; i < CELL_PER_BUCKET; ++i){
            if(buckets[pos].keys[i] == key){
                if(nowtime - buckets[pos].timestamps[i] >= time_threshold)
                    return 0;
                else
                    return buckets[pos].freq[i];
            }
        }
        return 0;
    }

    bool insert(uint32_t key, double nowtime){
        int pos = CalculatePos(key); 
        int empty_i = -1;
        int min_freq = __INT32_MAX__; 
        int min_i = -1;
        for(int i = 0; i < CELL_PER_BUCKET; ++i){
            if(buckets[pos].keys[i] == key){
                if(nowtime - buckets[pos].timestamps[i] >= time_threshold){
                    // buckets[pos].freq[i] = 1;
                    buckets[pos].freq[i] = __max(1, buckets[pos].negvotes / CELL_PER_BUCKET);
                    buckets[pos].timestamps[i] = nowtime;
                    buckets[pos].negvotes = 0;
                }
                else{
                    buckets[pos].freq[i]++;
                    buckets[pos].timestamps[i] = nowtime;
                }
                return true;
            }
            else if(buckets[pos].keys[i] == 0 || (nowtime - buckets[pos].timestamps[i] >= time_threshold)){
                if(empty_i == -1)
                    empty_i = i;
            }
            else if(buckets[pos].freq[i] < min_freq){
                min_freq = buckets[pos].freq[i];
                min_i = i;
            }
        }
        if(empty_i != -1){
            buckets[pos].keys[empty_i] = key;
            // buckets[pos].freq[empty_i] = 1;
            buckets[pos].freq[empty_i] = __max(1, buckets[pos].negvotes / CELL_PER_BUCKET);
            buckets[pos].timestamps[empty_i] = nowtime;
            buckets[pos].negvotes = 0;
            return true;
        }
        if(min_freq < buckets[pos].negvotes){
            buckets[pos].keys[min_i] = key;
            // buckets[pos].freq[min_i] = 1;
            buckets[pos].freq[min_i] = __max(1, buckets[pos].negvotes / CELL_PER_BUCKET);
            buckets[pos].timestamps[min_i] = nowtime;
            buckets[pos].negvotes = 0;
            return true;
        }
        buckets[pos].negvotes++;
        return false;
    }

private:
    inline int CalculatePos(uint32_t key){
        return CalculateBucketPos(key) % bucket_num; 
        // return MurmurHash3_x86_32((const char*)&key, sizeof(uint32_t), seed) % bucket_num; 
        // return bhash.run((const char*)&key, sizeof(uint32_t)) % bucket_num;
    }
};



#endif // _HOTSKETCH_H_