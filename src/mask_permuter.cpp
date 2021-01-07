//
// Written by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <cmath>

#include "mask_permuter.h"




mask_permuter::mask_permuter() = default;

mask_permuter::mask_permuter(int nspan, int nweight) {
    reset_sw(nspan, nweight);
    //init span and weight from arguments
    init_mask(); //initialize mask as defined by span and weight
}

mask_permuter::~mask_permuter() = default;

/*initializes mask as smallest lexicographical
 * sequence.*/
void mask_permuter::init_mask(){
    int infloc = 0;
    base_mask = {};
    for (int i = 0; i < span/2; i++){
        if(infloc >= (span/2) - (weight/2)) {
            base_mask.push_back(1);
            infloc++;
        }
        else {
            base_mask.push_back(0);
            infloc++;
        }
    }
}

void mask_permuter::reset_sw(int nspan, int nweight) {
    span = nspan;
    weight = nweight;
    init_mask();
}

/*Public function that returns the coco mask
 * of the next permutation*/
bool mask_permuter::get_next(std::vector<int> &ovec){
    bool check = update_permpos();
    ovec = curr_permpos;
    return check;
}

bool mask_permuter::get_next(unsigned char* msk, std::vector<int> &ovec){
    bool check = update_permpos();
    for(int i=0; i<curr_permpos.size(); i++){
        msk[i] = curr_permpos[i];
    }
    ovec = curr_permpos;
    return check;
}

unsigned int mask_permuter::get_permNum(){
    int hspan = span/2;
    int hweight = weight/2;
    return tgamma(hspan+1)/(tgamma(hweight+1)*tgamma((hspan-hweight)+1));
}

bool mask_permuter::update_permpos(){
    bool check;
    if (!init) {
        check = permuter(); //update curr_permpos with next permutation
    }
    else {
        mask_mkr(base_mask);
        init = false;
        check = true;
    }
    return check;
}

void mask_permuter::mask_mkr(std::vector<int> mask){
    int mod; //modifies generation of symmetric part
    curr_permpos = {};
    int i = 0;
    for (const int &t: mask) {
        if (t == 1) {
            curr_permpos.push_back(i);
            i++;
        }
        else {
            i++;
        }
    }
    if (span%2 != 0) {
        curr_permpos.push_back((span/2));
        mod = 0;
    }
    else{
        mod = 1;
    }
    for (int k=(weight/2)-1; k>=0; k--){
        curr_permpos.push_back((span/2)+((span/2)-mod)-curr_permpos[k]);
    }
}


/* Generates next permutation of base_mask
 * until lexicographically largest sequence
 * is met.*/
bool mask_permuter::permuter(){
    if(std::next_permutation(base_mask.begin(), base_mask.end())) {
        mask_mkr(base_mask);
        return true;
    }
    else{
        return false;
    }
}

/*debug function to print out a vector*/
void mask_permuter::show_vec(std::vector<int> vec) {
    for (const int tt: vec) {
        std::cout << tt << " ";
    }
    std::cout << std::endl;
}
