#pragma once

#ifdef __cplusplus
extern "C" {
#endif
    typedef void* (*func_allocate)(size_t nbytes);
    typedef void  (*func_deallocate)(void* pointer);

    void set_allocator(func_allocate falloc, func_deallocate fdealloc);

    void* allocate(size_t size_in_bytes);

    void  deallocate(void* pointer);

#ifdef __cplusplus
}
#endif
