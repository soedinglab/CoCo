//
// Written by Anton Robert Georg Farr on 01.12.20. <anton.farr@gmail.com>
//

#ifndef MASK_PERMUTER_MASK_PERMUTER_H
#define MASK_PERMUTER_MASK_PERMUTER_H

#include <vector>


class mask_permuter {

    std::vector<int> base_mask;
    std::vector<int> curr_permpos;
    int span;
    int weight;
    bool init = true;

public:
    mask_permuter();

    mask_permuter(int span, int weight);

    ~mask_permuter();

    bool get_next(std::vector<int> &ovec);

    bool get_next(unsigned char* arr, std::vector<int> &ovec);

    long int get_permNum();

    void reset_sw(int span, int weight);

private:
    void init_mask();

    void mask_mkr(std::vector<int> mask);

    bool permuter();

    bool update_permpos();

    long long int fact(long long int num);

    //debug function
    static void show_vec(std::vector<int> vec);


};


#endif //MASK_PERMUTER_MASK_PERMUTER_H
