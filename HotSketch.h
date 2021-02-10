#ifndef _HOTSKETCH_H_
#define _HOTSKETCH_H_


#define USING_SIMD

#include "param.h"
#include <immintrin.h>

// #include "murmur3.h"
// #include "BOBHash32.h"

class HotSketch{
public:
    Bucket* buckets;
    double time_threshold;
    uint32_t bucket_num;
    // uint32_t seed;
    // BOBHash32 bhash; 

    HotSketch(uint32_t cell_num_, double time_threshold_){
        bucket_num = (cell_num_ + CELL_PER_BUCKET - 1) / CELL_PER_BUCKET;
        time_threshold = time_threshold_;
        buckets = new(std::align_val_t{64}) Bucket[bucket_num];
        memset(buckets, 0, bucket_num*sizeof(Bucket));
        // seed = rand()%MAX_PRIME32; 
        // bhash.initialize(seed);
    }

    ~HotSketch(){
        delete [] buckets;
    }

    uint32_t query(uint32_t key, double nowtime){
#ifdef USING_SIMD
        int pos = CalculatePos(key); 
        
        const __m512i item = _mm512_set1_epi32(key);                                                                                                                                                                                             
        __m512i *keys_p = (__m512i*)(buckets[pos].keys);
        int matched = 0;
        matched = _mm512_cmpeq_epi32_mask(item, *keys_p);
        if(matched != 0){ 
            int matched_index = __tzcnt_u32((uint32_t)matched);
            if(matched_index < CELL_PER_BUCKET-1 && nowtime - buckets[pos].timestamps[matched_index] < time_threshold)
                return buckets[pos].freq[matched_index];
        }
        return 0;

#else
        int pos = CalculatePos(key); 
        
        for(int i = 0; i < CELL_PER_BUCKET-1; ++i){
            if(buckets[pos].keys[i] == key){
                if(nowtime - buckets[pos].timestamps[i] >= time_threshold)
                    return 0;
                else
                    return buckets[pos].freq[i];
            }
        }
        return 0;
#endif
    }

    bool insert(uint32_t key, double nowtime){
#ifdef USING_SIMD
        int pos = CalculatePos(key);

        const __m512i item = _mm512_set1_epi32(key);                                                                                                                                                                                             
        __m512i *keys_p = (__m512i*)(buckets[pos].keys);
        int matched = 0;
        matched = _mm512_cmpeq_epi32_mask(item, *keys_p);


        if(matched != 0){ 
            int matched_index = __tzcnt_u32((uint32_t)matched);
            if(nowtime - buckets[pos].timestamps[matched_index] >= time_threshold){
                buckets[pos].freq[matched_index] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
                buckets[pos].timestamps[matched_index] = nowtime;
                buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
            }
            else{
                buckets[pos].freq[matched_index]++;
                buckets[pos].timestamps[matched_index] = nowtime;
            }
            return true;
        }

        const __m512i zero_item = _mm512_set1_epi32(0);
        matched = _mm512_cmpeq_epi32_mask(zero_item, *keys_p);
        if(matched != 0){
            int matched_index = __tzcnt_u32((uint32_t)matched);
            if(matched_index < CELL_PER_BUCKET-1){
                buckets[pos].keys[matched_index] = key;
                buckets[pos].freq[matched_index] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
                buckets[pos].timestamps[matched_index] = nowtime;
                buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
                return true;
            }
        }


        int timeout_index = -1;
        for(int i = 0; i < CELL_PER_BUCKET-1; ++i){
            if(nowtime - buckets[pos].timestamps[i] >= time_threshold){
                timeout_index = i;
                break;
            }
        }
        if(timeout_index != -1){
            buckets[pos].keys[timeout_index] = key;
            buckets[pos].freq[timeout_index] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
            buckets[pos].timestamps[timeout_index] = nowtime;
            buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
            return true;
        }

        const uint32_t mask_base = 0x7fffffff;
        __m512 mask = (__m512)_mm512_set_epi32(mask_base, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

        const __m512 *freqs = (__m512*)(buckets[pos].freq);
        __m512 results = _mm512_or_ps(*freqs, mask);

        __m256i low_256 = _mm256_castps_si256(_mm512_extractf32x8_ps(results, 0));
        __m256i high_256 = _mm256_castps_si256(_mm512_extractf32x8_ps(results, 1));
        __m256i x_256 = _mm256_min_epi32(low_256, high_256);

        __m128i low_128 = _mm_castps_si128(_mm256_extractf128_ps((__m256)x_256, 0));
        __m128i high_128 = _mm_castps_si128(_mm256_extractf128_ps((__m256)x_256, 1));
        __m128i x_128 = _mm_min_epi32(low_128, high_128);

        __m128i min1 = _mm_shuffle_epi32(x_128, _MM_SHUFFLE(0,0,3,2));
	    __m128i min2 = _mm_min_epi32(x_128, min1);
	    __m128i min3 = _mm_shuffle_epi32(min2, _MM_SHUFFLE(0,0,0,1));
	    __m128i min4 = _mm_min_epi32(min2, min3);
	    int min_freq = _mm_cvtsi128_si32(min4);
        
        const __m512i freq_item = _mm512_set1_epi32(min_freq);
        int freq_matched = _mm512_cmpeq_epi32_mask(freq_item, (__m512i)results);
        int min_index = __tzcnt_u32((uint32_t)freq_matched);

        if(min_freq < buckets[pos].freq[CELL_PER_BUCKET-1]){
            buckets[pos].keys[min_index] = key;
            buckets[pos].freq[min_index] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
            buckets[pos].timestamps[min_index] = nowtime;
            buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
            return true;
        }
        buckets[pos].freq[CELL_PER_BUCKET-1]++;
        return false;

#else
        int pos = CalculatePos(key); 
        int empty_i = -1;
        int min_freq = __INT32_MAX__; 
        int min_i = -1;
        for(int i = 0; i < CELL_PER_BUCKET-1; ++i){
            if(buckets[pos].keys[i] == key){
                if(nowtime - buckets[pos].timestamps[i] >= time_threshold){
                    // buckets[pos].freq[i] = 1;
                    buckets[pos].freq[i] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
                    buckets[pos].timestamps[i] = nowtime;
                    buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
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
            buckets[pos].freq[empty_i] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
            buckets[pos].timestamps[empty_i] = nowtime;
            buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
            return true;
        }
        if(min_freq < buckets[pos].freq[CELL_PER_BUCKET-1]){
            buckets[pos].keys[min_i] = key;
            // buckets[pos].freq[min_i] = 1;
            buckets[pos].freq[min_i] = __max(1, buckets[pos].freq[CELL_PER_BUCKET-1] / CELL_PER_BUCKET);
            buckets[pos].timestamps[min_i] = nowtime;
            buckets[pos].freq[CELL_PER_BUCKET-1] = 0;
            return true;
        }
        buckets[pos].freq[CELL_PER_BUCKET-1]++;
        return false;
#endif
    }

private:
    inline int CalculatePos(uint32_t key){
        return CalculateBucketPos(key) % bucket_num; 
        // return MurmurHash3_x86_32((const char*)&key, sizeof(uint32_t), seed) % bucket_num; 
        // return bhash.run((const char*)&key, sizeof(uint32_t)) % bucket_num;
    }
};



#endif // _HOTSKETCH_H_
