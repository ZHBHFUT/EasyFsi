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
//! @file       SocketCommunicator.cpp
//!             The implement of SocketCommunicator class.
//! @author     ZHANG Bing, zhangbing@hfut.edu.cn
//! @version    1.0
//! @copyright  2023, all rights reserved.
//! @data       2023-06-08
//!-------------------------------------------------------------

#ifdef __linux__
#include <unistd.h> // gethostname, getpid
#elif _WIN32
#include <Windows.h> // GetCurrentProcessId 
#endif
#include <sstream>
#include <thread>
#include <mutex>
#include <memory>  // unique_ptr
#include <bit>     // C++20, endian

#include "Assert.hpp"
#include "Logger.hpp"
#include "ByteSwap.hpp"
#include "SocketCommunicator.hpp"

#ifdef _WIN32
#pragma comment(lib, "Ws2_32.lib")
#endif

namespace EasyLib {

    static std::string GetHostName()
    {
        char buf[256] = { '\0' };
        if (gethostname(buf, sizeof(buf)) != 0) {
            error("failed get host name!");
        }
        return std::string(buf);
    }
    static std::string GetHostIP(const std::string& host_name)
    {
        auto host_entry = gethostbyname(host_name.c_str());
        if (!host_entry) {
            error("failed get host by name!");
            return std::string();
        }
        else {
            return std::string(inet_ntoa(*((struct in_addr*)host_entry->h_addr_list[0])));
        }
    }
    static long GetThisProcessId()
    {
#ifdef __linux__
        return static_cast<long>(getpid());
#elif _WIN32
        return static_cast<long>(GetCurrentProcessId());
#else
#error "Platform not supported"
#endif
    }

    static int GetLastSocketError()
    {
#ifdef __linux__
        return errno;
#elif _WIN32
        return WSAGetLastError();
#else
#error "Platform not supported"
#endif
    }
    static bool InitSocket()
    {
        static bool init_ = false;
        if (init_)return true;

#ifdef __linux__

#elif _WIN32
        //initialize socket resource
        WSADATA wsaData;

        //0--sucessfully, !=0 -- failed
        if (WSAStartup(MAKEWORD(1, 1), &wsaData) != 0) {
            error("InitSocket(): WSAStartup() failed. ERROR CODE = %d.", GetLastSocketError());
            //return false;
        }

        if (LOBYTE(wsaData.wVersion) != 1 ||
            HIBYTE(wsaData.wVersion) != 1) {
            WSACleanup();
            error("InitSocket(): invalid WSADATA");
            //return false;
        }
#else
#error "Platform not supported"
#endif
        init_ = true;
        return true;
    }
    static bool IsErrorReturn(int ret)
    {
#ifdef __linux__
        return ret < 0;
#elif _WIN32
        return ret == SOCKET_ERROR;
#else
#error "Platform not supported"
#endif
    }
    //static bool IsDissconnect()
    //{
    //    int err = GetLastSocketError();
    //#ifdef __linux__
    //    return (err == ECONNRESET) || (err == ENOTCONN);
    //#elif _WIN32
    //    return (err == WSAECONNRESET) || (err == WSAENOTCONN);
    //#else
    //#error "Platform not supported"
    //#endif
    //}
    static bool IsValidSocket(SOCKET sock)noexcept
    {
#ifdef __linux__
        return sock > 0;
#elif _WIN32
        return !(sock == INVALID_SOCKET);
#else
#error "Platform not supported"
#endif
    }
    static void CloseSocket(SOCKET sock)noexcept
    {
        if (IsValidSocket(sock)) {
#ifdef __linux__
            close(sock);
#elif _WIN32
            closesocket(sock);
#else
#error "Platform not supported"
#endif
        }
    }
    static void MakeServerAddress(const char* name, unsigned short port, sockaddr_in& addr)
    {
        memset(&addr, 0, sizeof(addr));

#ifdef __linux__
        hostent* server = gethostbyname(name);
        memcpy(
            (char*)server->h_addr,
            (char*)&addr.sin_addr.s_addr,
            server->h_length
        );
#elif _WIN32
        addr.sin_addr.s_addr = inet_addr(name);//htonl(INADDR_ANY);
#else
#error "Platform not supported"
#endif
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);
    }
    static void MakeLocalAddress(unsigned short port, sockaddr_in& addr)
    {
        memset(&addr, 0, sizeof(addr));

#ifdef __linux__
        addr.sin_addr.s_addr = INADDR_ANY;
#elif _WIN32
        addr.sin_addr.s_addr = htonl(INADDR_ANY); //inet_addr("localhost");//htonl(INADDR_ANY);
#else
#error "Platform not supported"
#endif
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port);
    }
    static int ReceiveData(SOCKET sock, void* data, int n_bytes)
    {
        int nrecv = 0;
        for (; nrecv < n_bytes;) {
            int n = ::recv(sock, (char*)data + nrecv, n_bytes - nrecv, 0);
            if (n > 0)nrecv += n;
            else return nrecv;
        }
        return nrecv;
    }

    static constexpr const unsigned short is_big_endian = std::endian::native == std::endian::big ? 1 : 0;

    struct DataHead
    {
        int tag{ 0 };
        int len{ 0 };
    };

    //---------------------------------------------
    // implement of SocketCommunicator
    //---------------------------------------------

    static constexpr const int sock_family   = AF_INET;
    static constexpr const int sock_type     = SOCK_STREAM;
    static constexpr const int sock_protocol = IPPROTO_TCP;

    SocketCommunicator::SocketCommunicator()
        :Communicator()
    {
        InitSocket();
    }
    SocketCommunicator::SocketCommunicator(SocketCommunicator&& comm)noexcept
        ://app_name_(std::move(comm.app_name_)),
        host_name_(std::move(comm.host_name_)),
        host_ip_(std::move(comm.host_ip_)),
        host_port_(comm.host_port_),
        is_connected_(comm.is_connected_),
        server_sock_(std::move(comm.server_sock_)),
        connections_(std::move(comm.connections_)),
        size_(comm.size_),
        rank_(comm.rank_),
        timeout_(comm.timeout_)
    {
        comm.server_sock_ = INVALID_SOCKET;
        comm.connections_.clear();
    }
    SocketCommunicator& SocketCommunicator::operator = (SocketCommunicator&& comm)noexcept
    {
        if (&comm != this) {
            for (auto& c : connections_)if (IsValidSocket(c.socket))CloseSocket(c.socket);
            if (IsValidSocket(server_sock_))CloseSocket(server_sock_);
            connections_.clear();
            server_sock_ = INVALID_SOCKET;

            //app_name_     = std::move(comm.app_name_);
            host_name_    = std::move(comm.host_name_);
            host_ip_      = std::move(comm.host_ip_);
            host_port_    = comm.host_port_;
            is_connected_ = comm.is_connected_;
            server_sock_  = std::move(comm.server_sock_);
            connections_  = std::move(comm.connections_);
            size_         = comm.size_;
            rank_         = comm.rank_;
            timeout_      = comm.timeout_;

            comm.server_sock_ = INVALID_SOCKET;
            comm.connections_.clear();
        }
        return *this;
    }
    SocketCommunicator::~SocketCommunicator()
    {
        for (auto& c : connections_)if (IsValidSocket(c.socket))CloseSocket(c.socket);
        if (IsValidSocket(server_sock_))CloseSocket(server_sock_);
    }

    void SocketCommunicator::init(int argc, const char** argv)
    {
        bool as_master = false;
        //const char* app_name = nullptr;
        const char* master_ip = nullptr;
        unsigned short port = 50001;
        int np = 0;
        int timeout_sec = 60;
        for (int i = 0; i < argc; ++i) {
            if (strcmp(argv[i], "-master") == 0)as_master = true;
            else if (strcmp(argv[i], "-np") == 0) {
                ++i;
                if (i >= argc) {
                    error("participator number is missing");
                    return;
                }
                np = atoi(argv[i]);
            }
            else if (strcmp(argv[i], "-ip") == 0) {
                ++i;
                if (i >= argc) {
                    error("master IP address is missing");
                    return;
                }
                master_ip = argv[i];
            }
            else if (strcmp(argv[i], "-port") == 0) {
                ++i;
                if (i >= argc) {
                    error("master IP address is missing");
                    return;
                }
                port = static_cast<unsigned short>(atoi(argv[i]));
            }
        }

        this->init(as_master, np, master_ip, port, timeout_sec);
    }

    void SocketCommunicator::init(bool as_master, int np, const char* master_ip/* = "127.0.0.1"*/, unsigned short port/* = 50001*/, int timeout_sec/* = 60*/)
    {
        debug("  initialize socket communicator\n");

        if (as_master && np < 1)error("invalid participator number: %d", np);

        for (auto& c : connections_)if (IsValidSocket(c.socket))CloseSocket(c.socket);
        connections_.clear();

        //app_name_  = app_name;
        host_name_ = GetHostName();
        host_ip_   = GetHostIP(host_name_);

        rank_ = as_master ? 0 : -1;
        size_ = as_master ? np : 0;
        timeout_ = timeout_sec;

        // master: gather information

        // init root server
        std::unique_ptr<std::thread> t;
        if (as_master) {
            host_ip_ = master_ip && *master_ip != '\0' ? master_ip : "127.0.0.1";
            host_port_ = port;
            size_ = np;
            t = std::unique_ptr<std::thread>(new std::thread([this]() {
                this->init_server0_();
                }));
        }

        // init client for root(master)
        init_client0_(master_ip, port, timeout_sec);

        // waiting for server initializing.
        if (t.get() && t->joinable())t->join();

        // modify server port
        for (int i = 0; i < size_; ++i) {
            auto& p = connections_[i];
            for (int j = i + 1; j < size_; ++j) {
                auto& q = connections_[j];
                if (p.remote_host_name == q.remote_host_name) {
                    q.remote_host_ip = p.remote_host_ip;
                    q.remote_host_port = p.remote_host_port + static_cast<unsigned short>(j - i);
                }
            }
        }
        host_name_ = connections_[rank_].remote_host_name;
        host_ip_ = connections_[rank_].remote_host_ip;
        host_port_ = connections_[rank_].remote_host_port;

        // init server 
        if (rank_ != 0) {
            t = std::unique_ptr<std::thread>(new std::thread([this]() {
                this->init_serverX_();
                }));
        }

        // connecting
        init_clientX_(timeout_sec);

        // waiting for server initializing.
        if (t.get() && t->joinable())t->join();

        if (server_sock_ != INVALID_SOCKET) {
            info("\nThis server = %s:%d\n", host_ip_.c_str(), (int)host_port_);
        }

        info(
            "\n%4s %20s %12s %5s %6s %5s\n",
            "RANK", "HOST", "IP", "PORT", "PID", "ENDIAN"
        );
        for (int i = 0; i < size_; ++i) {
            auto& c = connections_.at(i);
            info(
                "%4d %20s %12s %5d %6ld %5d %lld\n",
                i,
                //c.remote_app_name.c_str(),
                c.remote_host_name.c_str(),
                c.remote_host_ip.c_str(),
                (int)c.remote_host_port,
                c.remote_pid,
                (int)c.remote_is_big_endian,
                c.socket
            );
        }
    }

    void SocketCommunicator::send(const void* data, int count, DataType type, int dest_rank, int tag)
    {
        send_(nbyte_of_type(type) * count, data, dest_rank, tag);
    }
    void SocketCommunicator::recv(void* data, int count, DataType type, int src_rank, int tag)
    {
        recv_(nbyte_of_type(type) * count, data, src_rank, tag);

        // bit swap
        if (connections_.at(src_rank).remote_is_big_endian != is_big_endian) {
            if      (nbyte_of_type(type) == 2) {
                auto p = reinterpret_cast<uint16_t*>(data);
                for (int i = 0; i < count; ++i)p[i] = byteswap(p[i]);
            }
            else if (nbyte_of_type(type) == 4) {
                auto p = reinterpret_cast<uint32_t*>(data);
                for (int i = 0; i < count; ++i)p[i] = byteswap(p[i]);
            }
            else if (nbyte_of_type(type) == 8) {
                auto p = reinterpret_cast<uint64_t*>(data);
                for (int i = 0; i < count; ++i)p[i] = byteswap(p[i]);
            }
        }
    }
    
    void SocketCommunicator::disconnect()
    {
        if (!is_connected_)return;

        for (auto& c : connections_)if (IsValidSocket(c.socket))CloseSocket(c.socket);
        if (IsValidSocket(server_sock_))CloseSocket(server_sock_);

        connections_.clear();
        server_sock_ = INVALID_SOCKET;

        is_connected_ = false;
        //info("disconnected!\n");
    }

    bool SocketCommunicator::init_server0_()
    {
        if (rank_ != 0)return true;

        // close previous sockets
        for (auto& c : connections_)if (IsValidSocket(c.socket))CloseSocket(c.socket);
        connections_.clear();
        if (IsValidSocket(server_sock_))CloseSocket(server_sock_);
        server_sock_ = INVALID_SOCKET;
        
        //info("\ninitialize master: %s:%d\n", host_ip_.c_str(), (int)host_port_);

        // create socket
        server_sock_ = socket(sock_family, sock_type, sock_protocol);
        if (!IsValidSocket(server_sock_)) {
            error("Communicator_t::init(): socket() failed with error = %d", GetLastSocketError());
            return false;
        }

        // make address for server
        sockaddr_in addrSrv;
        MakeLocalAddress(host_port_, addrSrv);

        // Bind the socket.
        int iResult = bind(server_sock_, (sockaddr*)&addrSrv, sizeof(addrSrv));
        if (IsErrorReturn(iResult)) {
            int ierr = GetLastSocketError();
            CloseSocket(server_sock_);
            error("SocketServer::run(), bind() failed. ERROR CODE = %d", ierr);
            return false;
        }

        // listen
        iResult = listen(server_sock_, SOMAXCONN);
        if (IsErrorReturn(iResult)) {
            int ierr = GetLastSocketError();
            CloseSocket(server_sock_);
            error("SocketServer::run(), listen() failed. ERROR CODE = %d", ierr);
            return false;
        }

        // allocate connections
        connections_.resize(size_);
        auto& c0 = connections_.front();
        //c0.remote_app_name      = app_name_;
        c0.remote_host_name     = host_name_;
        c0.remote_host_ip       = host_ip_;
        c0.remote_host_port     = host_port_;
        c0.remote_is_big_endian = is_big_endian;
        c0.remote_pid           = GetThisProcessId();
        
        // Accept the connection.
        //debug("\n[%d] waiting for connecting ...", rank_);

        // 
        std::ostringstream os_participators;
        os_participators
            //<< c0.remote_app_name << '\n'  // app name
            << c0.remote_host_name << '\n' // host name
            << c0.remote_host_ip << '\n'   // ip
            << c0.remote_host_port << '\n' // port
            << c0.remote_is_big_endian << '\n'
            << c0.remote_pid << '\n' // pid
            ;

        //info("Participators:\n    [0]: %s\n", c0.remote_app_name.c_str());

        for (int client_rank = 1; client_rank < size_; ++client_rank) {
            auto& c = connections_[client_rank];
            c.socket = accept(server_sock_, NULL, NULL);
            if (!IsValidSocket(c.socket)) {
                info("failed\n");
                int ierr = GetLastSocketError();
                CloseSocket(server_sock_);
                error("accept() failed. ERROR CODE = %d", ierr);
                return false;
            }

            // send and receive endian
            char big_endian = static_cast<char>(is_big_endian);
            ::send(c.socket, &big_endian, sizeof(big_endian), 0);
            big_endian = '\0';
            ::recv(c.socket, &big_endian, sizeof(big_endian), 0);

            // send size and rank from master,
            int buf[2] = { size_, client_rank };
            ::send(c.socket, (const char*)buf, sizeof(buf), 0);

            // send server IP
            int len = static_cast<int>(host_ip_.length());
            ::send(c.socket, (const char*)&len, sizeof(len), 0);
            ::send(c.socket, host_ip_.c_str(), len, 0);

            // send server name
            len = static_cast<int>(host_name_.length());
            ::send(c.socket, (const char*)&len, sizeof(len), 0);
            ::send(c.socket, host_name_.c_str(), len, 0);

            // receive participator info:
            //   app_name
            //   host_name
            //   ip
            //   port
            //   is_big_endian
            //   pid

            std::string str;
            len = 0;
            if (sizeof(len) != ReceiveData(c.socket, (char*)&len, sizeof(len))) {
                int ierr = GetLastSocketError();
                CloseSocket(c.socket);
                error("failed receive data. ERROR CODE = %d", ierr);
                return false;
            }
            if (big_endian != is_big_endian)len = byteswap(len); //? byte swap
            str.resize(len, '\0');
            if (len != ReceiveData(c.socket, (char*)str.data(), len)) {
                int ierr = GetLastSocketError();
                CloseSocket(c.socket);
                error("failed receive data. ERROR CODE = %d", ierr);
                return false;
            }

            os_participators << str;

            std::istringstream iss(str);
            //std::getline(iss, c.remote_app_name);
            std::getline(iss, c.remote_host_name);
            std::getline(iss, c.remote_host_ip);
            iss >> c.remote_host_port >> c.remote_is_big_endian >> c.remote_pid;

            //info("    [%d]: %s\n", client_rank, c.remote_app_name.c_str());
        }

        // send participator info from master
        //  app_name0
        //  host_name0
        //  host_ip0
        //  host_port0
        //  is_big_endian0
        //  pid0
        //  app_name1
        //  host_name1
        //  host_ip1
        //  host_port1
        //  is_big_endian1
        //  pid1
        //  ...
        //
        auto str = os_participators.str();
        int len = static_cast<int>(str.length());
        for (int client_rank = 1; client_rank < size_; ++client_rank) {
            auto& c = connections_.at(client_rank);
            ::send(c.socket, (const char*)&len, sizeof(len), 0);
            ::send(c.socket, str.c_str(), len, 0);
        }

        is_connected_ = true;
        //info("OK\n");

        return true;
    }
    bool SocketCommunicator::init_client0_(const char* ip_addr, unsigned short server_port, int timeout_sec)
    {
        if (rank_ == 0)return true;

        // create socket
        SOCKET sock = socket(sock_family, sock_type, sock_protocol);
        if (!IsValidSocket(sock)) {
            error("Communicator_t::init(): socket() failed with error = %d", GetLastSocketError());
            return false;
        }

        if (!ip_addr || *ip_addr == '\0')ip_addr = "127.0.0.1";

        // make server address
        sockaddr_in addrSrv;
        MakeServerAddress(ip_addr, server_port, addrSrv);

        info("\nconnecting %s:%d ", ip_addr, (int)server_port);

        // Connect to server.
        if (timeout_sec < 1)timeout_sec = 1;
        for (int i = 1; i <= timeout_sec; ++i) {
            int iResult = connect(sock, (sockaddr*)&addrSrv, sizeof(addrSrv));
            if (IsErrorReturn(iResult)) {
                if (i == timeout_sec) {
                    CloseSocket(sock);
                    error("Communicator::init(), connect() failed. ERROR CODE = %d", GetLastSocketError());
                    return false;
                }
                Sleep(1000);
                info("try connecting again %d/%d\n", i, 60);
            }
            else {
                break;
            }
        }

        // receive and send endian
        char big_endian = '\0';
        ::recv(sock, &big_endian, sizeof(big_endian), 0);
        char big = static_cast<char>(is_big_endian);
        ::send(sock, &big, sizeof(big), 0);

        // receive size and rank from root
        //  size, rank
        int buf[2] = { 0 };
        if (sizeof(buf) != ReceiveData(sock, buf, sizeof(buf))) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive size info. ERROR CODE = %d", ierr);
            return false;
        }
        if (big_endian != big) {
            buf[0] = byteswap(buf[0]);//? byte swap
            buf[0] = byteswap(buf[0]);
        }
        size_ = buf[0];
        rank_ = buf[1];

        // allocate participator data
        connections_.resize(size_);
        auto& c0 = connections_.at(0);
        
        // receive server ip
        int len = 0;
        if (sizeof(len) != ReceiveData(sock, (char*)&len, sizeof(len))) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive server IP length. ERROR CODE = %d", ierr);
            return false;
        }
        if (big_endian != big)len = byteswap(len); //? byte swap
        c0.remote_host_ip.resize(len);
        if (len != ReceiveData(sock, c0.remote_host_ip.data(), len)) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive server IP. ERROR CODE = %d", ierr);
            return false;
        }

        // receive server name
        if (sizeof(len) != ReceiveData(sock, (char*)&len, sizeof(len))) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive server name length. ERROR CODE = %d", ierr);
            return false;
        }
        if (big_endian != big)len = byteswap(len); //? byte swap
        c0.remote_host_name.resize(len);
        if (len != ReceiveData(sock, c0.remote_host_name.data(), len)) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive server name. ERROR CODE = %d", ierr);
            return false;
        }

        c0.remote_host_port = server_port;
        c0.socket = sock;

        // send participator info:
        //   app_name
        //   host name
        //   ip
        //   port
        //   is_big_endian
        //   pid

        std::ostringstream os_participators;
        os_participators
            //<< app_name_ << '\n'
            << host_name_ << '\n'
            << host_ip_ << '\n'
            << host_port_ << '\n'
            << is_big_endian << '\n'
            << GetThisProcessId() << '\n'
            ;
        auto s = os_participators.str();
        len = static_cast<int>(s.length());
        if (big_endian != is_big_endian)len = byteswap(len);
        ::send(sock, (const char*)&len, sizeof(len), 0);
        ::send(sock, s.c_str(), len, 0);

        // receive participators info from root
        if (sizeof(len) != ReceiveData(sock, &len, sizeof(len))) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive participator info length. ERROR CODE = %d", ierr);
            return false;
        }
        if (big_endian != big)len = byteswap(len); //? byte swap
        std::string str(len, '\0');
        if (len != ReceiveData(sock, str.data(), len)) {
            int ierr = GetLastSocketError();
            CloseSocket(sock);
            error("failed receive participator info. ERROR CODE = %d", ierr);
            return false;
        }

        // parse info:
        //  app_name
        //  host name
        //  ip
        //  port
        //  is_big_endian
        //  pid
        std::istringstream is_participators(str);
        for (auto& c : connections_) {
            //std::getline(is_participators, c.remote_app_name);
            std::getline(is_participators, c.remote_host_name);
            std::getline(is_participators, c.remote_host_ip);
            is_participators
                >> c.remote_host_port
                >> c.remote_is_big_endian
                >> c.remote_pid;
            std::getline(is_participators, s); // next line
        }

        //info("Participators:\n");
        //for (int i = 0; i < size_; ++i)
        //    info("    [%d]: %s\n", i, connections_[i].remote_app_name.c_str());
        //info("\n");

        is_connected_ = true;
        info("!!!OK!!!\n");

        return true;
    }
    bool SocketCommunicator::init_serverX_()
    {
        if (rank_ == 0)return true;

        //info("\ninitialize server: %s:%d\n", host_ip_.c_str(), (int)host_port_);

        // close previous sockets
        for (int i = rank_ + 1; i < size_ && i < connections_.size(); ++i) {
            auto& c = connections_.at(i);
            if (IsValidSocket(c.socket)) {
                CloseSocket(c.socket);
                c.socket = INVALID_SOCKET;
            };
        }
        if (IsValidSocket(server_sock_))CloseSocket(server_sock_);
        server_sock_ = INVALID_SOCKET;
        
        // create socket
        server_sock_ = socket(sock_family, sock_type, sock_protocol);
        if (!IsValidSocket(server_sock_)) {
            error("init_serverX_(): socket() failed with error = %d", GetLastSocketError());
            return false;
        }

        // make address for server
        sockaddr_in addrSrv;
        MakeLocalAddress(host_port_, addrSrv);

        // Bind the socket.
        int iResult = bind(server_sock_, (sockaddr*)&addrSrv, sizeof(addrSrv));
        if (IsErrorReturn(iResult)) {
            int ierr = GetLastSocketError();
            CloseSocket(server_sock_);
            error("SocketServer::run(), bind() failed. ERROR CODE = %d", ierr);
            return false;
        }

        // listen
        iResult = listen(server_sock_, SOMAXCONN);
        if (IsErrorReturn(iResult)) {
            int ierr = GetLastSocketError();
            CloseSocket(server_sock_);
            error("SocketServer::run(), listen() failed. ERROR CODE = %d", ierr);
            return false;
        }

        // Accept the connection.
        //debug("\n[%d] waiting for connecting ...", rank_);
        std::ostringstream oss;
        for (int client_rank = rank_ + 1; client_rank < size_; ++client_rank) {
            auto& c = connections_[client_rank];
            c.socket = accept(server_sock_, NULL, NULL);
            if (!IsValidSocket(c.socket)) {
                //info("failed\n");
                int ierr = GetLastSocketError();
                CloseSocket(server_sock_);
                error("accept() failed. ERROR CODE = %d", ierr);
                return false;
            }

            //info("\n%s connected", connections_[client_rank].remote_app_name.c_str());
        }

        is_connected_ = true;
        //info("\n!!!OK!!!\n");

        return true;
    }
    bool SocketCommunicator::init_clientX_(int timeout_sec)
    {
        if (rank_ == 0)return true;

        for (int server = 1; server < rank_; ++server) {
            auto& c = connections_.at(server);

            info("\n  connecting [%d] ...", server);

            // create socket
            c.socket = socket(sock_family, sock_type, sock_protocol);
            if (!IsValidSocket(c.socket)) {
                error("init_clientX_(): socket() failed with error = %d", GetLastSocketError());
                return false;
            }

            // make server address
            sockaddr_in addrSrv;
            MakeServerAddress(c.remote_host_ip.c_str(), c.remote_host_port, addrSrv);

            // Connect to server.
            if (timeout_sec < 1)timeout_sec = 1;
            for (int i = 1; i <= timeout_sec; ++i) {
                int iResult = connect(c.socket, (sockaddr*)&addrSrv, sizeof(addrSrv));
                if (IsErrorReturn(iResult)) {
                    if (i == timeout_sec) {
                        CloseSocket(c.socket);
                        error("init_clientX_(), connect() failed. ERROR CODE = %d", GetLastSocketError());
                        return false;
                    }
                    Sleep(1000);
                    info("\ntry connecting again %d/%d", i, timeout_sec);
                }
                else {
                    break;
                }
            }

            is_connected_ = true;
            info("!!!OK!!!\n");
        }

        return true;
    }

    bool SocketCommunicator::send_(int n_bytes, const void* data, int dest_rank, int tag)const
    {
        // check
        if (!is_connected_) {
            error("SocketCommunicator::send() failed. Invalid connection!");
            return false;
        }
        if (dest_rank < 0 || dest_rank >= size_) {
            error("SocketCommunicator::recv_() failed. Invalid source rank!");
            return false;
        }
        if (n_bytes <= 0)return true;
        if (dest_rank == rank_) {
            error("SocketCommunicator::send() failed. Unable send data to self!");
            return false;
        }

        //debug("[%d] send %d bytes to [%d]\n", rank_, n_bytes, dest_rank);

        auto& c = connections_.at(dest_rank);
        if (!IsValidSocket(c.socket)) {
            error("SocketCommunicator::send() failed. Invalid socket!");
            return false;
        }

        // send header: src, des, tag, len
        DataHead head{ tag, n_bytes };
        if (sizeof(head) != ::send(c.socket, (const char*)&head, sizeof(head), 0)) {
            warn("failed sending header!\n");
            ASSERT(false);
            return false;
        }

        // send data
        if (n_bytes != ::send(c.socket, (const char*)data, n_bytes, 0)) {
            warn("failed sending data!\n");
            ASSERT(false);
            return false;
        }

        return true;
    }
    bool SocketCommunicator::recv_(int n_bytes, void* data, int src_rank, int tag)
    {
        // check
        if (!is_connected_) {
            error("SocketCommunicator::recv() failed. Invalid connection!");
            return false;
        }
        if (src_rank < 0 || src_rank >= size_) {
            error("SocketCommunicator::recv() failed. Source rank out of range!");
            return false;
        }
        if (src_rank == rank_) {
            error("SocketCommunicator::recv() failed. Unable receive data from self!");
            return false;
        }
        auto& c = connections_.at(src_rank);
        if (!IsValidSocket(c.socket)) {
            error("SocketCommunicator::recv() failed. Invalid socket!");
            return false;
        }
        if (n_bytes <= 0)return true;

        //debug("[%d] recv %d bytes from [%d]\n", rank_, n_bytes, src_rank);

        // receive head
        DataHead head;
        if (sizeof(head) != ReceiveData(c.socket, &head, sizeof(head))) {
            warn("failed receiving header!\n");
            ASSERT(false);
            return false;
        }
        if (c.remote_is_big_endian != is_big_endian) {
            head.tag = byteswap(head.tag);
            head.len = byteswap(head.len);
        }
        //debug("  Head = %d %d\n", head.tag, head.len);
        if (head.tag != tag) {
            warn("Tag mismatch! Please check send/recv order!\n");
            ASSERT(false);
            return false;
        }
        if (head.len != n_bytes) {
            warn("Length not agree!\n");
            ASSERT(false);
            return false;
        }

        // receive data
        if (n_bytes != ReceiveData(c.socket, data, n_bytes)) {
            warn("failed receiving data!\n");
            ASSERT(false);
            return false;
        }

        return true;
    }
}
