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
        auto& f = fields_.emplace_back(Field{});
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
                return std::pair<bool, Field*>{true, & f};
        return std::pair<bool, Field*>{false, nullptr};
    }
    std::pair<bool, const Field*> ModelInterface::find_field(const char* field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, const Field*>{true, & f};
        return std::pair<bool, const Field*>{false, nullptr};
    }
    std::pair<bool, Field*> ModelInterface::find_field(const std::string& field_name)
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, Field*>{true, & f};
        return std::pair<bool, Field*>{false, nullptr};
    }
    std::pair<bool, const Field*> ModelInterface::find_field(const std::string& field_name)const
    {
        for (auto& f : fields_)
            if (f.info->name == field_name)
                return std::pair<bool, const Field*>{true, & f};
        return std::pair<bool, const Field*>{false, nullptr};
    }

}
