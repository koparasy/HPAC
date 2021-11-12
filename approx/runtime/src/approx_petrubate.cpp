//===--- approx_petrubate.cpp - runtime implementation of TAF ----------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
/// \file
/// This files is the driver of the approximate petrubation.
///
//===----------------------------------------------------------------------==//
#include <cfloat>
#include <fstream> 
#include <iostream>
#include <iomanip>
#include <memory>
#include <cstring>
#include <unordered_map>
#include <cmath>
#include <omp.h>
#include <random>

#include "approx_data_util.h"
#include "approx_internal.h"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

using namespace std;

char *petrubate_file = nullptr;

#define IGNORE 0
#define RECORD 1
#define PETRUBATE 2

enum PetrubationMethod{
 ABS = 1,
 REL,
 DIST,
};

struct approx_petru_var{
    long numElements;
    const char *type;
    char *varName;
    double error;
    uint8_t dir;
    uint8_t error_type;
    int position;
    approx_petru_var() : numElements(0), 
            type(nullptr), varName(nullptr), 
            error(0.0), dir(0), position(-1) {}
};

class PetrubateMap{
    private:
        // Key is the region pointer
        // Vector contains region names and indexes to approx_petru_var 
        unordered_map<uintptr_t, std::vector<std::pair<char*, int>>> appSpace;
        approx_petru_var **region;
        int *regionNumVars;
        char **regionNames;
        std::string descr_file;
        int execution_type;
        long region_mem_size;
        long last_index_size;
    public:
        PetrubateMap(const char *fName, const char *type) {
            descr_file.assign(fName,std::strlen(fName));
            region_mem_size = 1000; 
            last_index_size = 0;
            region = new approx_petru_var*[region_mem_size]();
            regionNumVars = new int[region_mem_size]();
            regionNames = new char*[region_mem_size]();
            if (std::strcmp(type,"PETRUBATE") == 0){
                execution_type = PETRUBATE;
            }
            else if ( std::strcmp(type, "RECORD") == 0){
                execution_type = RECORD;
            }
            if (execution_type == PETRUBATE){
                std::ifstream in_file(descr_file);
                json JDescr;
                in_file >> JDescr;
                for (json::iterator it = JDescr.begin(); it != JDescr.end(); ++it) {
                    char *region_name = new char[it.key().length()+1];
                    std::strcpy(region_name, it.key().c_str());
                    regionNames[last_index_size] = region_name;
                    regionNumVars[last_index_size] = it.value().size();
                    region[last_index_size] = new approx_petru_var[it.value().size()];
                    for (json::iterator it2 = it.value().begin(); 
                            it2 != it.value().end(); ++it2) {
                        int numElem = it2.value()[0];
                        double error = it2.value()[1];
                        string typeName = it2.value()[2].get<std::string>(); 
                        uint8_t error_type = it2.value()[3];
                        int position = it2.value()[4];
                        region[last_index_size][position].numElements = numElem; 
                        region[last_index_size][position].error = error; 
                        region[last_index_size][position].position = position; 
                        region[last_index_size][position].error_type = error_type; 
                        char *var_name = new char[it2.key().length()+1];
                        std::strcpy(var_name, it2.key().c_str());
                        region[last_index_size][position].varName = var_name;
                        char* cstr = new char [typeName.length()+1];
                        std::strcpy (cstr, typeName.c_str());
                        region[last_index_size][position].type= cstr;
                    }
                    last_index_size +=1;
                }
            }
        }
        ~PetrubateMap(){
            if (execution_type == RECORD){
                json JDescr;
                for (auto &it: appSpace){
                    for ( auto &vecIt : it.second){
                        std::string regionName(vecIt.first);
                        int index = vecIt.second;
                        int numVars = regionNumVars[index];
                        for (int i = 0; i < numVars; i++){
                            std::string tmpName(region[index][i].varName);
                            JDescr[regionName][tmpName] = {region[index][i].numElements, 
                                    region[index][i].error, region[index][i].type, region[index][i].error_type, i};
                            region[index][i].position = i;
                        }
                    }
                }
                std::ofstream record_out(descr_file);
                record_out << std::setw(4) << JDescr << std::endl;
            }
            for (int i = 0; i < last_index_size ; i++){
                delete [] regionNames[i];
                for (int k = 0; k < regionNumVars[i]; k++){
                    delete [] region[i][k].varName; 
                }
                delete [] region[i];
                region[i] = nullptr;
            }
            delete [] region;
            delete [] regionNames;
            region = nullptr;
            regionNames = nullptr;
            for ( auto it : appSpace ){
                for ( auto v : it.second ){
                    v.first = nullptr;
                }
            }
        }
    void petrubate(uintptr_t address, approx_var_info_t *arguments, int num_arguments, const char *label);
    void record(uintptr_t address, approx_var_info_t *arguments, int num_arguments, const char *label);
    int getExecutionType() { return execution_type; }
};

void 
PetrubateMap::record(uintptr_t address, approx_var_info_t *arguments, 
       int num_arguments, const char *label){
    auto &record = appSpace[address];
    int rid = -1;
    for (int index = 0; index < record.size(); index++){
        if (strcmp(record[index].first, label) == 0){
            rid = index;
        }
    }
    // We are in the record and we found our region.
    if (rid != -1)
        return;

    approx_petru_var *tmp = new approx_petru_var[num_arguments];
    for (int p = 0; p < num_arguments; p++){
        tmp[p].numElements = arguments[p].num_elem;
        tmp[p].error = 0;
        tmp[p].error_type = 1;
        tmp[p].dir = arguments[p].dir;
        tmp[p].type = getTypeName((ApproxType) arguments[p].data_type);
        int name_size = strlen(arguments[p].name) + 1;
        tmp[p].varName = new char[name_size];
        std::strncpy(tmp[p].varName, arguments[p].name, name_size);
    }

    if (last_index_size >= region_mem_size){
        approx_petru_var **tmp_region = new approx_petru_var*[2*region_mem_size];
        int *tmp_region_mem_size = new int[2*region_mem_size]();
        std::memcpy(tmp_region, region, sizeof(approx_petru_var*)*region_mem_size);
        std::memcpy(tmp_region_mem_size, regionNumVars, sizeof(int)*region_mem_size);
        region_mem_size *= 2;
        delete [] region;
        delete [] regionNumVars;
        region = tmp_region;
        regionNumVars = tmp_region_mem_size;
    }
    int len = strlen(label) +1;
    char* region_name = new char[len];
    std::strncpy(region_name, label, len);
    region[last_index_size] = tmp;
    regionNumVars[last_index_size] = num_arguments;
    int entry_index = last_index_size;
    last_index_size+=1;
    record.push_back(std::make_pair(region_name,entry_index));
}

void 
PetrubateMap::petrubate(uintptr_t accurate, approx_var_info_t *arguments, 
        int num_arguments, const char *label){
        int found = 0;
    int rId = -1;
    thread_local int prev_id = 0;
    for (int i = 0; i < last_index_size; i++){
        if (strcmp(label, regionNames[prev_id]) == 0){
            rId = i;
            break;
        }
        prev_id += 1;
        prev_id = prev_id % last_index_size;
    }
    prev_id = (prev_id +1) % last_index_size;
    if ( rId == -1 ){
        std::cout<< "Error region is not recorded however I am in petrubate-mode\n";
        exit(-1);
    }

    for (int p = 0; p < num_arguments; p++){
        approx_petru_var& var = region[rId][p];
        if ( var.error_type == PetrubationMethod::REL){
            petrubate_var(arguments[p].ptr, (ApproxType) arguments[p].data_type, 
                arguments[p].num_elem, var.error);
        }
        else if ( var.error_type == PetrubationMethod::ABS ){
            abs_error_var(arguments[p].ptr, (ApproxType) arguments[p].data_type, 
                arguments[p].num_elem, var.error);
        }
        else if ( var.error_type == PetrubationMethod::DIST){
            dist_error_var(arguments[p].ptr, (ApproxType) arguments[p].data_type, 
                arguments[p].num_elem, var.error);

        }
    }
}

PetrubateMap** allPetru = nullptr;
int maxPetruThreads = 36;
const char *pType;
const char *pName;

void deinitPetrubate(){
    if (allPetru){
        for (int i = 0; i < maxPetruThreads; i++){
            if ( allPetru[i] )
                delete allPetru[i];
            allPetru[i] = nullptr;
        }
    }
    delete[] allPetru;
    allPetru = nullptr;
}

void register_petrubate(const char *lName, const char *ltype){
    pName = lName;
    pType = ltype;
    if (!allPetru){
        allPetru = new PetrubateMap*[maxPetruThreads];
        for (int i = 0; i < maxPetruThreads; i++)
            allPetru[i] = nullptr;
    }
}

void petrubate(void (*accurate)(void *), approx_var_info_t *arguments, int num_arguments, const char *label){
    uintptr_t address = (uintptr_t) accurate;
    static thread_local int threadId = -1;
    static thread_local PetrubateMap *curr = nullptr;

    if (threadId == -1){
        if (omp_in_parallel())
            threadId = omp_get_thread_num();
        else
            threadId = 0;
    }

    if ( !curr ){
        curr = allPetru[threadId];
        if ( !curr ){
            curr = new PetrubateMap(pName, pType);
            allPetru[threadId] = curr;
        }
    }

    if ( curr->getExecutionType() == RECORD ){
        curr->record(address, arguments, num_arguments, label);
    }
    else if ( curr->getExecutionType() == PETRUBATE ){
        curr->petrubate(address, arguments, num_arguments, label);
    }
}
