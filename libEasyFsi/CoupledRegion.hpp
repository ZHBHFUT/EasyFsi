#pragma once
/*-------------------------------------------------------------------------
This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from
the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

2. Altered source versions must be plainly marked as such, and must not
   be misrepresented as being the original software.

3. This notice may not be removed or altered from any source distribution.
-------------------------------------------------------------------------*/

//!-------------------------------------------------------------
//! @file       CoupledRegion.hpp
//!             The definition of CoupledRegion class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include "ModelInterface.hpp"

namespace EasyLib {

    //! @brief Coupled region class
    class CoupledRegion : public ModelInterface
    {
    public:
        virtual ~CoupledRegion() = default;

        void load(const char* /*file*/)override {/* TBD */ }
        void save(const char* /*file*/)const override {/* TBD */ }

        int_l nnode()const final { /* TBD */ return 0; }
        int_l nelem()const final { /* TBD */ return 0; }

    private:
        // TBD

    };
}
