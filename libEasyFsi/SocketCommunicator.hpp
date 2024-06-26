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
//! @file       SocketCommunicator.hpp
//!             The definition SocketCommunicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#ifdef __linux__
#include <sys/types.h>
#include <sys/socket.h>
typedef int SOCKET;
#elif _WIN32
#define NOMINMAX
#define _WINSOCK_DEPRECATED_NO_WARNINGS
#include <winsock.h>
#else
#error "Platform not supported"
#endif

#include <vector>

#include "Communicator.hpp"

namespace EasyLib {

    class SocketCommunicator : public Communicator
    {
    public:
        using Communicator::send;
        using Communicator::recv;

        //! @brief Default constructor.
        SocketCommunicator();

        SocketCommunicator(const SocketCommunicator&) = delete;
        SocketCommunicator& operator = (const SocketCommunicator&) = delete;

        SocketCommunicator(SocketCommunicator&& comm)noexcept;
        SocketCommunicator& operator = (SocketCommunicator&&)noexcept;

        //! @brief Destructor.
        virtual ~SocketCommunicator();

        void set_constant(const char* /*name*/, int /*value*/)final {}
        void set_constant(const char* /*name*/, void* /*value*/)final {}
        void set_function(const char* /*name*/, void* /*func*/)final {}

        //! @brief Initialize the socket communicator.
        //! @param as_root     Whether to run as the root node
        //! @param np          Number of applications that participate in communication.
        //! @param master_ip   The IP address of the root node, default using the local IP address：127.0.0.1
        //! @param port        Port number of the root node，default is 50001.
        //! @param timeout_sec Timeout value in seconds for connecting to server, default is 60s.
        void init(bool as_root, int np, const char* master_ip = "127.0.0.1", unsigned short port = 50001, int timeout_sec = 60);

        void init(int argc, const char** argv);

        //! @brief Get rank of current application, zero-based value.
        //! @return 
        int rank()const noexcept final { return rank_; }

        //! @brief Get the number of applications that participate in communication.
        int size()const noexcept final { return size_; }

        void send(const void* data, int count, DataType type, int dest_rank, int tag)final;
        void recv(void* data, int count, DataType type, int src_rank, int tag)final;
        void async_send(const void* data, int count, DataType type, int dest_rank, int tag)final;
        void wait()final;

        void disconnect() final;

    private:
        bool send_(int n_bytes, const void* data, int dest_rank, int tag)const;
        bool recv_(int n_bytes, void* data, int src_rank, int tag);

        bool init_server0_();
        bool init_client0_(const char* ip_addr, unsigned short server_port, int timeout_sec);
        bool init_serverX_();
        bool init_clientX_(int timeout_sec);

    private:
        struct Connection
        {
            //std::string    remote_app_name;
            std::string    remote_host_name;
            std::string    remote_host_ip;
            unsigned short remote_host_port{ 0 };
            unsigned short remote_is_big_endian{ 0 };
            long           remote_pid{ 0 };
            SOCKET         socket{ INVALID_SOCKET };
        };
        using cvec = std::vector<Connection>;

        //std::string    app_name_;                      // 当前应用名称
        std::string    host_name_;                     // 当前应用所在主机名称
        std::string    host_ip_;                       // 当前应用所在主机IP地址
        unsigned short host_port_{ 50001 };            // 当前应用服务端端口号
        bool           is_connected_{ false };         // 是否已建立链接
        SOCKET         server_sock_{ INVALID_SOCKET }; // 作为服务端的套接字，对最后一个参与应用无效
        cvec           connections_;
        int            size_{ 1 };
        int            rank_{ 0 };
        int            timeout_{ 60 };
    };
}
