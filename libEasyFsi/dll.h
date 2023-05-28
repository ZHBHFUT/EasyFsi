#pragma once

#if defined(USE_DLL)
#if defined(LIBRARY_EXPORT) // inside DLL
#   define DLLAPI __declspec(dllexport)
#else // outside DLL
#   define DLLAPI __declspec(dllimport)
#endif  // LIBRARY_EXPORT
#else
#   define DLLAPI
#endif

