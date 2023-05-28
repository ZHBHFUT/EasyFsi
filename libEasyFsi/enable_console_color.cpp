#if defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64
#include <Windows.h>
#endif

void enable_console_color()
{
#if defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64
    static bool init = false;
    if (init)return;

    HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD dwMode = 0;
    GetConsoleMode(hOut, &dwMode);
    dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
    SetConsoleMode(hOut, dwMode);

    // References:
    //SetConsoleMode() and ENABLE_VIRTUAL_TERMINAL_PROCESSING?
    //https://stackoverflow.com/questions/38772468/setconsolemode-and-enable-virtual-terminal-processing

    // Windows console with ANSI colors handling
    // https://superuser.com/questions/413073/windows-console-with-ansi-colors-handling
    init = true;
#else
    return;
#endif
}
