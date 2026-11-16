// BEST H-matrix solver: the watch on available memory (memory.h).
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "memory.h"

#include <chrono>

#include "mat.h"

namespace bh {

namespace detail {
std::atomic<bool> g_memory_low( false );
}

namespace {
std::atomic<double> g_low_seen( 0 );
}

double memory_low_seen() { return g_low_seen.load(); }

MemoryGuard::MemoryGuard( double floor ) : floor_( floor ) {
   detail::g_memory_low = false;
   g_low_seen = 0;
   if ( floor_ <= 0 ) {
      return;
   }
   // The first reading here, so that a node already below the floor stops the solver at once.
   const double avail = double( mem_available() );
   if ( avail > 0 ) {
      if ( avail < floor_ ) {
         g_low_seen = avail;
         detail::g_memory_low = true;
      }
      th_ = std::thread( &MemoryGuard::watch, this );
   }
}

MemoryGuard::~MemoryGuard() {
   if ( th_.joinable() ) {
      {
         std::lock_guard<std::mutex> lk( m_ );
         quit_ = true;
      }
      cv_.notify_all();
      th_.join();
   }
}

void MemoryGuard::watch() {
   std::unique_lock<std::mutex> lk( m_ );
   while ( !quit_ ) {
      const double avail = double( mem_available() );
      if ( avail > 0 && avail < floor_ && !detail::g_memory_low ) {
         g_low_seen = avail;
         detail::g_memory_low = true;
      }
      cv_.wait_for( lk, std::chrono::milliseconds( 50 ), [&] { return quit_; } );
   }
}

}  // namespace bh
