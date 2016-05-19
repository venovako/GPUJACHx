// posix_threads_helper.cu: POSIX threads helper.

#include "posix_threads_helper.hpp"

#include "my_utils.hpp"

static TLS
struct
{
  int pix;
  int tix;
} threadVars = { INVALID_THREAD, MAIN_THREAD };

void getProcessAffinity(cpu_set_t &cpuset) throw()
{
  CPU_ZERO(&cpuset);
  SYSI_CALL(sched_getaffinity(static_cast<pid_t>(0), sizeof(cpuset), &cpuset));
}

void setProcessAffinity(const cpu_set_t &cpuset) throw()
{
  SYSI_CALL(sched_setaffinity(static_cast<pid_t>(0), sizeof(cpuset), &cpuset));  
}

void getThreadAffinity(cpu_set_t &cpuset) throw()
{
  CPU_ZERO(&cpuset);
  SYSI_CALL(pthread_getaffinity_np(pthread_self(), sizeof(cpuset), &cpuset));
}

void setThreadAffinity(const cpu_set_t &cpuset) throw()
{
  SYSI_CALL(pthread_setaffinity_np(pthread_self(), sizeof(cpuset), &cpuset));
}

struct ThreadAttrs
{
  PThreadFn tfn;
  void *arg;
  int pix;
  int tix;
  unsigned Long aff;
};

static void* threadStarter(ThreadAttrs *tha) throw()
{
  PThreadFn tfn = tha->tfn;
  void *arg = tha->arg;

  cpu_set_t
    &setdst = *(new cpu_set_t),
    &setsrc = *(new cpu_set_t);

  threadVars.pix = tha->pix;
  threadVars.tix = tha->tix;

  getThreadAffinity(setsrc);
  if (tha->aff) {
    if (*static_cast<unsigned Long*>(memcpy(&setdst, &setsrc, sizeof(cpu_set_t))) &= tha->aff)
      setThreadAffinity(setdst);
    else {
      (void)snprintf(err_msg, err_msg_size, "Affinity mask empty for thread %d", threadVars.tix);
      DIE(err_msg);
    }
  }
  else {
    const int cnt = CPU_COUNT(&setsrc);
    int tix = threadVars.tix;
    if (cnt > tix) {
      int cpu = 0;
      do {
        if (CPU_ISSET(cpu, &setsrc) && !tix--)
          break;
      } while (++cpu);
      CPU_ZERO(&setdst);
      CPU_SET(cpu, &setdst);
      setThreadAffinity(setdst);
    }
    else {
      (void)snprintf(err_msg, err_msg_size, "Thread %d oversubscribing", threadVars.tix);
      DIE(err_msg);
    }
  }

  delete &setsrc;
  delete &setdst;
  delete tha;

  return tfn(arg);
}

void createThread
(
 pthread_t &tid,
 const PThreadFn tfn,
 void *const arg,
 const int tix,
 const unsigned Long aff
) throw()
{
  if (!tfn)
    DIE("NULL tfn");
  if (0 > tix)
    DIE("tix < 0");

  ThreadAttrs &tha = *(new ThreadAttrs);
  tha.tfn = tfn;
  tha.arg = arg;
  tha.pix = getThreadTix();
  tha.tix = tix;
  tha.aff = aff;

  sigset_t sig_all, sig_old;
  SYSI_CALL(sigfillset(&sig_all));
  SYSI_CALL(sigemptyset(&sig_old));
  SYSI_CALL(pthread_sigmask(SIG_BLOCK, &sig_all, &sig_old));
  SYSI_CALL(pthread_create(&tid, static_cast<const pthread_attr_t*>(NULL), reinterpret_cast<PThreadFn>(threadStarter), &tha));
  SYSI_CALL(pthread_sigmask(SIG_SETMASK, &sig_old, &sig_all));
}

int getThreadPix() throw()
{
  return threadVars.pix;
}

int getThreadTix() throw()
{
  return threadVars.tix;
}

bool barrierWait(pthread_barrier_t &bar) throw()
{
  const int err = pthread_barrier_wait(&bar);
  const bool ret = (PTHREAD_BARRIER_SERIAL_THREAD == err);
  SYSI_CALL(ret ? 0 : err);
  return ret;
}
