//
// Written by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#ifndef MASK_PERMUTER_MASK_PERMUTER_H
#define MASK_PERMUTER_MASK_PERMUTER_H

#include <vector>


class mask_permuter {

    std::vector<int> base_mask;
    std::vector<int> curr_permpos;
    std::vector<unsigned int> rand_masks;
    int span;
    int weight;
    bool init = true;
    bool randomize = false;

public:
    mask_permuter();

    mask_permuter(int span, int weight);

    ~mask_permuter();

    bool get_next(std::vector<int> &ovec);

    bool get_next(unsigned char* arr, std::vector<int> &ovec);

    unsigned int get_permNum();

    void reset_sw(int span, int weight);

    void set_rand(unsigned int start, unsigned int stop, unsigned int maskNum);

private:
    void init_mask();

    void mask_mkr(std::vector<int> mask);

    bool permuter();

    bool update_permpos();

    //debug function
    static void show_vec(std::vector<int> vec);


};


#endif //MASK_PERMUTER_MASK_PERMUTER_H
