#ifndef _TRACE_H_
#define _TRACE_H_

#include "param.h"


void loadCAIDA18(vector<pair<uint32_t, double> > &vec, map<uint32_t, uint32_t> &mp, const char* filename, 
        double time_threshold, int hot_threshold, int read_num = 0){
    set<uint32_t> idset;
    map<uint32_t, double> tmap; 
	map<uint32_t, uint32_t> idmap;
    printf("Open %s \n", filename);
    FILE* pf = fopen(filename, "rb");
    if(!pf){
        printf("%s not found!\n", filename);
        exit(-1);
    }
    vec.clear();
    char trace[30];
    uint32_t ret = 0;
    int cnt = 0;

    uint32_t tkey; 
    double ttime; 
    while(fread(trace, 1, 21, pf)){
        tkey = *(uint32_t*) (trace+5); 
        ttime = *(double*) (trace+13); 
        vec.push_back(pair<uint32_t, double>(tkey, ttime)); 
        idset.insert(tkey);
        map<uint32_t, double>::iterator iter; 
        if((iter = tmap.find(tkey)) == tmap.end()){
            tmap[tkey] = ttime;
            mp[tkey] = 1;
        }
        else{
            if(ttime - iter->second >= time_threshold){
                mp[tkey] = 1;
                iter->second = ttime;
            }
            else{
                mp[tkey]++;
                iter->second = ttime;
            }
        }
        if(idmap.find(tkey) == idmap.end())
            idmap[tkey] = 1;
        else 
            idmap[tkey]++;
        if(++cnt == read_num)
            break;
        if(cnt % 1000000 == 1){
            printf("%u : %lf\n", tkey, ttime);
        }
    }
    fclose(pf);

    int flowlet_num = 0;
    int hotlet_num = 0;
    for(auto it = mp.begin(); it != mp.end(); ++it){
        if(ttime - tmap[it->first] >= time_threshold)
            it->second = 0;
        if(it->second > 0)
            flowlet_num++;
        if(it->second > hot_threshold)
            hotlet_num++;
    }


    printf("load %d items, %lu distinct items. \n", cnt, idset.size());
    printf("flowlet number: %d, hotlet number: %d\n", flowlet_num, hotlet_num);
    printf("mp.size: %lu\n", mp.size());
}


#endif // _TRACE_H_