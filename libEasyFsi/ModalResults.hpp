#pragma once
#include <vector>
#include "DynamicArray.hpp"
#include "TinyVector.hpp"

namespace EasyLib {

    struct ModalResults
    {
        int                                ngrid{ 0 };
        int                                nmode{ 0 };
        std::vector<int>                   ids;
        std::vector<TinyVector<double, 3>> coords;
        std::vector<double>                freq; // frequency in units of rad/s
        std::vector<double>                mass; // 
        DynamicArray<double, 3>            phi;  // [nmode,ngrid,6]

        void clear();
        void load(const char* file);
        void save(const char* file)const;
    };
}