#pragma once

#if defined(__clang__) || defined(__GNUC__)
#define _force_inline_ __attribute__((always_inline)) inline
#define _force_inline_ [[gnu::always_inline]] inline
#elif defined(_MSC_VER)
#ifdef _DEBUG
#define _force_inline_ inline
#else
#pragma warning(error: 4714)
#define _force_inline_ __forceinline
#endif
#else
#error Unsupported compiler
#endif
