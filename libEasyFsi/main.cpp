#include <iostream>

#include "SocketCommunicator.hpp"

int main(int argc, char** argv)
{
    //std::cout << "Hello World!" << std::endl;

    int np = 1;
    const char* name = "unnamed";
    bool master = false;

    if (argc <= 1)return -1;

    for (int i = 1; i < argc; ++i) {
        if      (strcmp(argv[i], "-np") == 0) {
            ++i;
            np = atoi(argv[i]);
        }
        else if (strcmp(argv[i], "-name") == 0) {
            ++i;
            name = argv[i];
        }
        else if (strcmp(argv[i], "-master") == 0) {
            master = true;
        }
    }
    
    EasyLib::SocketCommunicator comm;
    comm.init(name, master, np);

    std::cout << "np = " << comm.size() << ", rank = " << comm.rank() << '\n';
    std::cout.flush();

    std::string msg_send = "Hello";
    std::string msg_recv;

    for (int i = 0; i < comm.size(); ++i) {
        if (i == comm.rank())continue;
        comm.send(msg_send, i, 0);
        std::cout << "send \""<< msg_send << "\" to [" << i << "]\n";
    }

    for (int i = 0; i < comm.size(); ++i) {
        if (i == comm.rank())continue;
        comm.recv(msg_recv, i, 0);
        std::cout << "recv from [" << i << "] " << msg_recv << '\n';
        std::cout.flush();
    }

    comm.disconnect();

    return 0;
}
