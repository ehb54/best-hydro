// BEST H-matrix solver: a watch on the node's available memory.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
//
// While a MemoryGuard exists, a thread of its own reads MemAvailable (/proc/meminfo) every 50 ms
// and raises a flag once it falls below a floor.  The long loops of the build, the H-LU and GMRES
// check the flag (memory_low(), one atomic load) and stop early; the caller then frees what it
// built and reports the error.  So a job that needs more memory than the estimate said stops
// cleanly before the node runs out, instead of making it swap or calling the OOM killer, which on
// a shared node hits everyone.
#pragma once

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

namespace bh {

class MemoryGuard {
 public:
   // floor in bytes; floor <= 0, or no /proc/meminfo: no thread, and the flag stays down.  The first
   // reading is made here.
   explicit MemoryGuard( double floor );
   ~MemoryGuard();
   MemoryGuard( const MemoryGuard& ) = delete;
   MemoryGuard& operator=( const MemoryGuard& ) = delete;

 private:
   void watch();
   double floor_;
   bool quit_ = false;
   std::mutex m_;
   std::condition_variable cv_;
   std::thread th_;
};

namespace detail {
extern std::atomic<bool> g_memory_low;
}

// True once a guard has seen MemAvailable below its floor (until the next guard starts).
inline bool memory_low() { return detail::g_memory_low.load( std::memory_order_relaxed ); }
// The MemAvailable reading, in bytes, that raised the flag.
double memory_low_seen();

}  // namespace bh
