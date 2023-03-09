#pragma once

#include "Communicator.hpp"

namespace EasyLib {

    class Application;

    class MPXContext
    {
    public:

        void add_application(Application& app);

        void init();

    private:
        Communicator* comm_{ nullptr };
    };
}
