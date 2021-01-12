//
// Written by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

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
//TODO: Make start and end of mask always informative
void mask_permuter::init_mask(){
    //informative location counter (1)
    int infloc = 0;
    base_mask = {};
    //generate mask
    for(int i = 0; i < span/2; i++){
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

void mask_permuter::set_rand(unsigned int start, unsigned int stop, unsigned int maskNum){
    unsigned int maxPerm = get_permNum();
    //check if stop is set and if not set it to maximum possible
    if(stop == 0){
        stop = maxPerm;
    }
    //number of unique masks to draw can't be larger than specified range
    if (maskNum>(stop-start+1)){
        throw "Error";
    }
    //random number is generated from hardware
    std::random_device rd;
    //generator is seeded
    std::mt19937 gen(rd());
    //range of randomly generated numbers set
    std::uniform_int_distribution<> distr(start,stop);
    //vector for sampled mask IDs
    std::vector<unsigned int> masks;
    unsigned int rnum;
    //draw masks and make sure that each one is unique in vector "masks".
    for (int i=0; i<maskNum; ++i){
        while(true) {
            rnum = distr(gen);
            if(!(std::find(masks.begin(), masks.end(), rnum) != masks.end())){
                masks.push_back(rnum);
                break;
            }
        }
    }
    //sort mask IDs since they are generated in increasing order
    std::sort(masks.begin(), masks.end());
    rand_masks = masks;
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
    unsigned long long int perms = tgamma(hspan+1)/(tgamma(hweight+1)*tgamma((hspan-hweight)+1));
    return perms;
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

/*converts mask (gapped kmer) composed of 1s and 0s into
 * positions indicating the informative locations (1)*/
void mask_permuter::mask_mkr(std::vector<int> mask){
    //modifies generation of symmetric part for even and uneven masks
    int mod;
    curr_permpos = {};
    //init counter for informative positions
    int i = 0;
    /*iterate through permuted left half of mask
     * and add positions of informative locations
     * to curr_permpos*/
    for (const int &t: mask) {
        if (t == 1) {
            curr_permpos.push_back(i);
            i++;
        }
        else {
            i++;
        }
    }
    /*add middle informative position if span is uneven
     * and set mod to calculate symmetric positions accordingly */
    if (span%2 != 0) {
        curr_permpos.push_back((span/2));
        mod = 0;
    }
    else{
        mod = 1;
    }
    //add symmetric right side of mask (gapped kmer)
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
