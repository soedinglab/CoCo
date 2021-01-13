//
// Created by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#ifndef MASK_PERMUTER_MASK_PERMUTER_H
#define MASK_PERMUTER_MASK_PERMUTER_H

#include <vector>


class mask_permuter {

    std::vector<int> base_mask;
    std::vector<int> curr_permpos;
    std::vector<int> rand_masks;
    int rand_pos = 0;
    int span;
    int weight;
    int perm_count = 0;
    int max_perm = 0;
    bool initial_mask = true;
    bool set_random = false;

public:
    mask_permuter();

    mask_permuter(int span, int weight);

    ~mask_permuter();

    void set_rand(int start, int stop, int maskNum);

    bool get_next(std::vector<int> &ovec);

    bool get_next(unsigned char* arr, std::vector<int> &ovec);

    int get_maxPerm();

    int get_permCount();

private:
    void reset_sw(int span, int weight);

    void init_mask();

    void mask_mkr(std::vector<int> mask);

    bool permuter();

    bool skip_perm();

    bool update_permpos();

    int calc_maxPerm();

    //debug function
    static void show_vec(std::vector<int> vec);


};


#endif //MASK_PERMUTER_MASK_PERMUTER_H
