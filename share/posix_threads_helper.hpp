#ifndef POSIX_THREADS_HELPER_HPP
#define POSIX_THREADS_HELPER_HPP

#include "defines.hpp"

#ifdef _WIN32
#error POSIX threading only
#else // POSIX
#include <pthread.h>
#include <sched.h>
#endif // _WIN32

typedef void*(*PThreadFn)(void*);

extern void getProcessAffinity(cpu_set_t &cpuset) throw();
extern void setProcessAffinity(const cpu_set_t &cpuset) throw();

extern void getThreadAffinity(cpu_set_t &cpuset) throw();
extern void setThreadAffinity(const cpu_set_t &cpuset) throw();

extern void createThread(pthread_t &tid, const PThreadFn tfn, void *const arg, const int tix, const unsigned Long aff) throw();

#ifndef MAIN_THREAD
#define MAIN_THREAD -1
#else // MAIN_THREAD
#error MAIN_THREAD not definable externally
#endif // !MAIN_THREAD

#ifndef INVALID_THREAD
#define INVALID_THREAD INT_MIN
#else // INVALID_THREAD
#error INVALID_THREAD not definable externally
#endif // !INVALID_THREAD

extern int getThreadPix() throw();
extern int getThreadTix() throw();

extern bool barrierWait(pthread_barrier_t &bar) throw();

#endif // !POSIX_THREADS_HELPER_HPP
