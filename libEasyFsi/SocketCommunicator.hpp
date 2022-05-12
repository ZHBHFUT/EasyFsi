#pragma once
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

        SocketCommunicator();

        virtual ~SocketCommunicator();

        //! @brief 初始化套接字通信子
        //! @param app_name    应用程序名称
        //! @param as_master   是否作为主节点运行
        //! @param np          参与通信的应用个数
        //! @param master_ip   主节点IP地址，默认为本机IP：127.0.0.1
        //! @param start_port  主节点端口号，默认为 50001
        //! @param timeout_sec 超时秒数，如果客户端在此规定时间内无法链接到参与的应用则退出，默认为60秒
        void init(const char* app_name, bool as_master, int np, const char* master_ip = "127.0.0.1", unsigned short start_port = 50001, int timeout_sec = 60);

        int rank()const noexcept final { return rank_; }

        int size()const noexcept final { return size_; }

        bool send(const int_l* data, int count, int dest_rank, int tag)final;
        bool send(const int_g* data, int count, int dest_rank, int tag)final;
        bool send(const double* data, int count, int dest_rank, int tag)final;
        bool send(const float* data, int count, int dest_rank, int tag)final;
        bool send(const char* data, int count, int dest_rank, int tag)final;

        bool recv(int_l* data, int count, int src_rank, int tag)final;
        bool recv(int_g* data, int count, int src_rank, int tag)final;
        bool recv(double* data, int count, int src_rank, int tag)final;
        bool recv(float* data, int count, int src_rank, int tag)final;
        bool recv(char* data, int count, int src_rank, int tag)final;

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
            std::string    remote_app_name;
            std::string    remote_host_name;
            std::string    remote_host_ip;
            unsigned short remote_host_port{ 0 };
            unsigned short remote_is_big_endian{ 0 };
            long           remote_pid{ 0 };
            SOCKET         socket{ INVALID_SOCKET };
        };
        using cvec = std::vector<Connection>;

        std::string    app_name_; // 当前应用名称
        std::string    host_name_;// 当前应用所在主机名称
        std::string    host_ip_;  // 当前应用所在主机IP地址
        unsigned short host_port_{ 50001 }; // 当前应用服务端端口号
        bool           is_connected_{ false };// 是否已建立链接
        SOCKET         server_sock_{ INVALID_SOCKET }; // 作为服务端的套接字，对最后一个参与应用无效
        cvec           connections_;
        int            size_{ 1 };
        int            rank_{ 0 };
        int            timeout_{ 60 };
    };
}
