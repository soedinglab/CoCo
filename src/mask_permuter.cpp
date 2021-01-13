//
// Created by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#include <bits/stdc++.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>

#include "mask_permuter.h"



//public functions
mask_permuter::mask_permuter() = default;

mask_permuter::mask_permuter(int nspan, int nweight) {
    //init span and weight from arguments
    reset_sw(nspan, nweight);
    max_perm = calc_maxPerm();
    //initialize gapped kmer as defined by span and weight
    init_mask();
}

mask_permuter::~mask_permuter() = default;

void mask_permuter::set_rand(int start, int stop, int maskNum){
    unsigned int maxPerm =get_maxPerm();
    //check if permutation is still first
    if(perm_count != 0){
        throw "Error: 'set_rand' can not be used if the mask was already permuted!\n";
    }
    //check if stop is set and if not set it to maximum possible
    if(stop == 0){
        stop = maxPerm;
    }
    //number of unique gapped kmers to draw can't be larger than specified range
    if (maskNum>(stop-start+1)){
        throw "Error: 'mask_permuter.cpp' number of unique masks to draw can't be larger than specified range\n";
    }
    //random number is generated from hardware
    std::random_device rd;
    //generator is seeded
    std::mt19937 gen(rd());
    //range of randomly generated numbers set
    std::uniform_int_distribution<> distr(start,stop);
    //vector for sampled gapped kmer IDs
    std::vector<int> masks;
    int rnum;
    //draw gapped kmers and make sure that each one is unique in vector "masks".
    for (int i=0; i<maskNum; ++i){
        while(true) {
            rnum = distr(gen);
            if(!(std::find(masks.begin(), masks.end(), rnum) != masks.end())){
                masks.push_back(rnum);
                break;
            }
        }
    }
    //sort gapped kmer IDs since they are generated in increasing order
    std::sort(masks.begin(), masks.end());
    set_random = true;
    /*set intitial_mask = false so that next_permutation is not called
     * for the first gapped kmer (0) */
    if (masks[0] != 0) {
        initial_mask = false;
    }
    //debug output
    std::cerr << "Randomly chosen gapped kmers: ";
    for (int &ii: masks){
        std::cerr << ii << " ";
    }
    std::cerr << "\n";
    rand_masks = masks;

}

/*Public function that returns the coco mask
 * of the next permutation either only as a vector
 * or as well an array (overloaded)*/
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

//return number of possible permutations
int mask_permuter::get_maxPerm(){
    return max_perm;
}

//returns current permutation number
int mask_permuter::get_permCount(){
    return perm_count;
}


//private functions
void mask_permuter::reset_sw(int nspan, int nweight) {
    span = nspan;
    weight = nweight;
    init_mask();
}

/*initializes gapped kmer as smallest lexicographical
 * sequence.*/
//TODO: Make start and end of mask always informative
void mask_permuter::init_mask(){
    //informative location counter
    int uinfloc = 0;
    //always first position informative
    base_mask = {1};
    //generate gapped kmer (minus one position)
    for(int i = 0; i < (span/2)-1; i++){
        if(uinfloc >= (span/2) - (weight/2)) {
            base_mask.push_back(1);
            uinfloc++;
        }
        else {
            base_mask.push_back(0);
            uinfloc++;
        }
    }
}

/*converts mask (gapped kmer) composed of 1s and 0s into
 * positions indicating the informative locations (1)*/
void mask_permuter::mask_mkr(std::vector<int> mask){
    //modifies generation of symmetric part for even and uneven gapped kmers
    int mod;
    curr_permpos = {};
    //first position is always informative
    //curr_permpos.push_back(0);

    //init counter for informative positions
    int i = 0;
    /*iterate through permuted left half of gapped kmer
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
bool mask_permuter::permuter() {
    //skip next_permutation() call for the first gapped kmer (pos 0)
    if(initial_mask){
        mask_mkr(base_mask);
        initial_mask = false;
        return true;
    }
    //get next permutation of base mask (gapped kmer) and apply mask_mkr
    else{
        if(std::next_permutation(base_mask.begin()+1, base_mask.end())){
            mask_mkr(base_mask);
            perm_count += 1;
            return true;
        }
        else{
            return false;
        }
    }
}

//skip a permutation (necessary to draw random gapped kmers)
bool mask_permuter::skip_perm(){
    initial_mask = false;
    if(std::next_permutation(base_mask.begin(), base_mask.end())) {
        perm_count += 1;
        return true;
    }
    else{
        return false;
    }
}

//gets called by get_next() to generate next
bool mask_permuter::update_permpos(){
    bool check;
    unsigned int rmask;
    //if random is set all gapped kmers that are not drawn by set_rand are skipped
    if(set_random){
        if (rand_pos == rand_masks.size()){
            return false;
        }
        rmask = rand_masks[rand_pos];
        rand_pos += 1;
        while(true){
            if (rmask-1 == perm_count || (rmask == 0 && rmask == perm_count)){
                break;
            }
            skip_perm();
        }
    }
    check=permuter(); //update curr_permpos with next permutation
    return check;
}

//calculate maximum number of permutations
int mask_permuter::calc_maxPerm(){
    int hspan = (span/2)-1;
    int hweight = (weight/2)-1;
    int perms = tgamma(hspan+1)/(tgamma(hweight+1)*tgamma((hspan-hweight)+1));
    return perms;
}

/*debug function to print out a vector*/
void mask_permuter::show_vec(std::vector<int> vec) {
    for (const int tt: vec) {
        std::cerr << tt << " ";
    }
    std::cerr << std::endl;
}
