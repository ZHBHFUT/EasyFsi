#pragma once
#include "Index.hpp"
#include "Field.hpp"

namespace EasyLib {

    class Communicator;

    class ModelInterface
    {
    public:
        ModelInterface() = default;
        ModelInterface(const ModelInterface&) = default;
        ModelInterface& operator = (const ModelInterface&) = default;

        ModelInterface(ModelInterface&& mi)noexcept;
        ModelInterface& operator = (ModelInterface&& mi)noexcept;

        virtual ~ModelInterface() = default;

        virtual void load(const char* file) = 0;
        virtual void save(const char* file)const = 0;

        virtual int_l nnode()const = 0;
        virtual int_l nelem()const = 0;

        //! @brief Get name of the ModelInterface.
        inline auto& name()const noexcept { return name_; }

        //! @brief Get user ID of this ModelInterface.
        inline int user_id()const noexcept { return user_id_; }

        void set_name(std::string sname);

        void set_user_id(int id);

        void register_field(const FieldInfo& fd);

        void remove_all_field();

        inline auto& fields()const noexcept { return fields_; }

        Field&       field(const char* field_name);
        const Field& field(const char* field_name)const;
        Field&       field(const std::string& field_name);
        const Field& field(const std::string& field_name)const;

        std::pair<bool, Field*>       find_field(const char* field_name);
        std::pair<bool, const Field*> find_field(const char* field_name)const;
        std::pair<bool, Field*>       find_field(const std::string& field_name);
        std::pair<bool, const Field*> find_field(const std::string& field_name)const;

        friend class Communicator;
    private:
        std::string name_;         //! name of this model interface.
        int         user_id_{ 0 }; //! user ID of this model interface.
        Fields      fields_;       //! fields list.
    };
}
