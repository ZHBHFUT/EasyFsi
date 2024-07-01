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
//! @file       Memory.cpp
//!             The implement of memory management functions.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#include "Memory.hpp"
#include <memory>

namespace EasyLib {

    //Memory::Memory(const Memory& other)
    //    :base_(other.size_in_bytes_ > 0 ? ::operator new(other.size_in_bytes_) : nullptr),
    //    size_in_bytes_(other.size_in_bytes_)
    //{
    //    if (base_)std::memcpy(base_, other.base_, size_in_bytes_);
    //}
    //Memory::Memory(Memory&& other)noexcept
    //    :base_(other.base_),
    //    size_in_bytes_(other.size_in_bytes_)
    //{
    //    other.base_ = nullptr;
    //    other.size_in_bytes_ = 0;
    //}
    //Memory::Memory(size_t size)
    //    :base_(size>0 ? ::operator new(size) : nullptr),
    //    size_in_bytes_(size)
    //{
    //    if (base_)std::memset(base_, 0, size_in_bytes_);
    //}
    //Memory::~Memory()
    //{
    //    if (base_)::operator delete(base_);
    //}
    //Memory& Memory::operator=(const Memory& other)
    //{
    //    if (&other == this)return *this;
    //
    //    if (size_in_bytes_ != other.size_in_bytes_) {
    //        if (base_)::operator delete(base_);
    //        base_ = other.size_in_bytes_ > 0 ? ::operator new(other.size_in_bytes_) : nullptr;
    //        size_in_bytes_ = other.size_in_bytes_;
    //        if(base_)std::memcpy(base_, other.base_, size_in_bytes_);
    //    }
    //    return *this;
    //}
    //Memory& Memory::operator=(Memory&& other)noexcept
    //{
    //    if (&other == this)return *this;
    //    if (base_)::operator delete(base_);
    //    base_ = other.base_;
    //    size_in_bytes_ = other.size_in_bytes_;
    //    other.base_ = nullptr;
    //    other.size_in_bytes_ = 0;
    //    return *this;
    //}
    //
    //void Memory::reallocate(size_t new_size)
    //{
    //    if (size_in_bytes_ == new_size)return;
    //    void* p = new_size > 0 ? ::operator new(new_size) : nullptr;
    //    if (base_ && p)std::memcpy(p, base_, size_in_bytes_);
    //    if (base_)::operator delete(base_);
    //    base_ = p;
    //    size_in_bytes_ = new_size;
    //}
    //void Memory::deallocate()
    //{
    //    if (base_)::operator delete(base_);
    //    base_ = nullptr;
    //    size_in_bytes_ = 0;
    //}
}

static size_t counter          = 0;
static func_allocate   alloc   = ::operator new;
static func_deallocate dealloc = ::operator delete;

extern "C" void set_allocator(func_allocate falloc, func_deallocate fdealloc)
{
    alloc   = falloc ? falloc : (func_allocate)::operator new;
    dealloc = fdealloc ? fdealloc : (func_deallocate)::operator delete;
}

extern "C" void* allocate(size_t size_in_bytes)
{
    return
        size_in_bytes > 0
        ? (alloc ? alloc(size_in_bytes) : ::operator new(size_in_bytes))
        : nullptr;
}

extern "C" void  deallocate(void* pointer)
{
    if (pointer) {
        if (dealloc)dealloc(pointer);
        else ::operator delete(pointer);
    }
}
