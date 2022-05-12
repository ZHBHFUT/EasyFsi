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

        //! @brief ��ʼ���׽���ͨ����
        //! @param app_name    Ӧ�ó�������
        //! @param as_master   �Ƿ���Ϊ���ڵ�����
        //! @param np          ����ͨ�ŵ�Ӧ�ø���
        //! @param master_ip   ���ڵ�IP��ַ��Ĭ��Ϊ����IP��127.0.0.1
        //! @param start_port  ���ڵ�˿ںţ�Ĭ��Ϊ 50001
        //! @param timeout_sec ��ʱ����������ͻ����ڴ˹涨ʱ�����޷����ӵ������Ӧ�����˳���Ĭ��Ϊ60��
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

        std::string    app_name_; // ��ǰӦ������
        std::string    host_name_;// ��ǰӦ��������������
        std::string    host_ip_;  // ��ǰӦ����������IP��ַ
        unsigned short host_port_{ 50001 }; // ��ǰӦ�÷���˶˿ں�
        bool           is_connected_{ false };// �Ƿ��ѽ�������
        SOCKET         server_sock_{ INVALID_SOCKET }; // ��Ϊ����˵��׽��֣������һ������Ӧ����Ч
        cvec           connections_;
        int            size_{ 1 };
        int            rank_{ 0 };
        int            timeout_{ 60 };
    };
}
