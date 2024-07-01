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
//! @file       ModelInterface.cpp
//!             The implement of ModelInterface class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------
#include "ModelInterface.hpp"

namespace EasyLib {

    ModelInterface::ModelInterface(ModelInterface&& mi)noexcept
        :fields_(std::move(mi.fields_))
    {}

    ModelInterface& ModelInterface::operator = (ModelInterface&& mi)noexcept
    {
        fields_ = std::move(mi.fields_);
        mi.fields_.clear();
        return *this;
    }

    void ModelInterface::set_name(std::string sname)
    {
        name_ = sname;
        name_.erase(0, name_.find_first_not_of("\r\t\n ")); // trim left
        name_.erase(name_.find_last_not_of("\r\t\n ") + 1); // trim right
    }

    void ModelInterface::set_user_id(int id)
    {
        user_id_ = id;
    }

    void ModelInterface::remove_all_field()
    {
        fields_.clear();
    }

    void ModelInterface::register_field(const FieldInfo& fd)
    {
        for (auto& f : fields_)if (f.info == &fd)return;
#if __cplusplus >= 201703L
        auto& f = fields_.emplace_back(Field{});
#else
        fields_.push_back(Field{});
        auto& f = fields_.back();
#endif
        f.info = &fd;
        f.data.resize(fd.location == NodeCentered ? nnode() : nelem(), fd.ncomp, 0);
    }

    Field& ModelInterface::field(const char* field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field \"%s\" is not found!", field_name);
        return *(Field*)nullptr;
    }
    const Field& ModelInterface::field(const char* field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field \"%s\" is not found!", field_name);
        return *(const Field*)nullptr;
    }
    Field& ModelInterface::field(const std::string& field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field \"%s\" is not found!", field_name.c_str());
        return *(Field*)nullptr;
    }
    const Field& ModelInterface::field(const std::string& field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)return f;
        error("field \"%s\" is not found!", field_name.c_str());
        return *(Field*)nullptr;
    }
    std::pair<bool, Field*> ModelInterface::find_field(const char* field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, Field*>(true, &f);
        return std::pair<bool, Field*>{false, nullptr};
    }
    std::pair<bool, const Field*> ModelInterface::find_field(const char* field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::make_pair(true, & f);
        return std::pair<bool, const Field*>{false, nullptr};
    }
    std::pair<bool, Field*> ModelInterface::find_field(const std::string& field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, Field*>(true, & f);
        return std::pair<bool, Field*>{false, nullptr};
    }
    std::pair<bool, const Field*> ModelInterface::find_field(const std::string& field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, const Field*>(true, & f);
        return std::pair<bool, const Field*>{false, nullptr};
    }

}
