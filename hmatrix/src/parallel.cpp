// BEST H-matrix solver: the thread pool behind parallel.h.
// Copyright (c) 2026 E. H. Brookes.  MIT License; see LICENSE.
#include "parallel.h"

#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace bh {

namespace {

thread_local bool t_in_parallel = false;

int hardware_threads() {
   unsigned n = std::thread::hardware_concurrency();
   return n > 0 ? int( n ) : 1;
}

class Pool {
 public:
   Pool() { start( hardware_threads() ); }
   ~Pool() { stop(); }

   int size() const { return n_; }

   void resize( int n ) {
      std::lock_guard<std::mutex> run( run_m_ );
      stop();
      start( n );
   }

   void run( const std::function<void( int )>& job ) {
      if ( n_ <= 1 ) {
         call( job, 0 );
      } else {
         std::lock_guard<std::mutex> run( run_m_ );
         {
            std::lock_guard<std::mutex> lk( m_ );
            job_ = &job;
            pending_ = n_ - 1;
            error_ = nullptr;
            ++generation_;
         }
         cv_job_.notify_all();
         call( job, 0 );
         std::unique_lock<std::mutex> lk( m_ );
         cv_done_.wait( lk, [&] { return pending_ == 0; } );
         job_ = nullptr;
      }
      if ( error_ ) {
         std::exception_ptr e = error_;
         error_ = nullptr;
         std::rethrow_exception( e );
      }
   }

 private:
   void call( const std::function<void( int )>& job, int thread ) {
      t_in_parallel = true;
      try {
         job( thread );
      } catch ( ... ) {
         std::lock_guard<std::mutex> lk( m_ );
         if ( !error_ ) {
            error_ = std::current_exception();
         }
      }
      t_in_parallel = false;
   }

   void worker( int thread ) {
      long seen = 0;
      for ( ;; ) {
         const std::function<void( int )>* job;
         {
            std::unique_lock<std::mutex> lk( m_ );
            cv_job_.wait( lk, [&] { return quit_ || generation_ != seen; } );
            if ( quit_ ) {
               return;
            }
            seen = generation_;
            job = job_;
         }
         call( *job, thread );
         {
            std::lock_guard<std::mutex> lk( m_ );
            if ( --pending_ == 0 ) {
               cv_done_.notify_one();
            }
         }
      }
   }

   void start( int n ) {
      n_ = std::max( 1, n );
      quit_ = false;
      generation_ = 0;
      for ( int t = 1; t < n_; ++t ) {
         workers_.emplace_back( &Pool::worker, this, t );
      }
   }

   void stop() {
      {
         std::lock_guard<std::mutex> lk( m_ );
         quit_ = true;
      }
      cv_job_.notify_all();
      for ( std::thread& w : workers_ ) {
         w.join();
      }
      workers_.clear();
   }

   int n_ = 1;
   std::vector<std::thread> workers_;
   std::mutex m_, run_m_;
   std::condition_variable cv_job_, cv_done_;
   const std::function<void( int )>* job_ = nullptr;
   long generation_ = 0;
   int pending_ = 0;
   bool quit_ = false;
   std::exception_ptr error_;
};

Pool& pool() {
   static Pool p;
   return p;
}

// The work-stealing pool behind fork_join.  Slot 0 belongs to the thread outside the pool that
// calls fork_join (the one running BEST), slots 1 .. n-1 to the pool's threads.  Each slot holds
// the tasks its thread has offered: the owner adds and takes back at the end, other threads steal
// from the front, so the oldest (largest) tasks are stolen first.
struct Task {
   const std::function<void()>* fn = nullptr;
   std::atomic<int> state{ 0 };  // 0 offered, 1 running, 2 done
   std::exception_ptr error;
};

thread_local int t_slot = -1;  // this thread's slot, -1 outside the pool

class Stealer {
 public:
   Stealer() { start( hardware_threads() ); }
   ~Stealer() { stop(); }

   void resize( int n ) {
      stop();
      start( n );
   }

   void fork_join( const std::function<void()>& a, const std::function<void()>& b ) {
      int self = t_slot;
      if ( self < 0 ) {
         if ( n_ <= 1 || in_parallel() ) {
            a();
            b();
            return;
         }
         self = 0;
      }
      Task tb;
      tb.fn = &b;
      Slot& s = *slots_[ self ];
      {
         std::lock_guard<std::mutex> lk( s.m );
         s.q.push_back( &tb );
         s.count.store( int( s.q.size() ) );
      }
      ++offered_;
      if ( sleepers_.load() > 0 ) {
         { std::lock_guard<std::mutex> lk( sleep_m_ ); }
         sleep_cv_.notify_one();
      }
      std::exception_ptr ea;
      try {
         a();
      } catch ( ... ) {
         ea = std::current_exception();
      }
      bool mine = false;
      {
         std::lock_guard<std::mutex> lk( s.m );
         if ( !s.q.empty() && s.q.back() == &tb ) {
            s.q.pop_back();
            s.count.store( int( s.q.size() ) );
            mine = true;
         }
      }
      if ( mine ) {
         --offered_;
         run( &tb );
      } else {
         // Taken by another thread: run offered tasks until it is done.
         unsigned seed = unsigned( self ) * 2654435761u + 1u;
         while ( tb.state.load( std::memory_order_acquire ) != 2 ) {
            if ( !steal_and_run( self, seed ) ) {
               std::this_thread::yield();
            }
         }
      }
      if ( ea ) {
         std::rethrow_exception( ea );
      }
      if ( tb.error ) {
         std::rethrow_exception( tb.error );
      }
   }

 private:
   struct Slot {
      std::mutex m;
      std::deque<Task*> q;
      std::atomic<int> count{ 0 };
   };

   static void run( Task* t ) {
      t->state.store( 1 );
      try {
         ( *t->fn )();
      } catch ( ... ) {
         t->error = std::current_exception();
      }
      t->state.store( 2, std::memory_order_release );
   }

   // Steals the oldest task of another slot (or, last, this one's) and runs it.
   bool steal_and_run( int self, unsigned& seed ) {
      seed = seed * 1664525u + 1013904223u;
      const int start = int( ( seed >> 8 ) % unsigned( n_ ) );
      for ( int k = 0; k <= n_; ++k ) {
         const int v = k < n_ ? ( start + k ) % n_ : self;
         if ( k < n_ && v == self ) {
            continue;
         }
         Slot& s = *slots_[ v ];
         if ( s.count.load() == 0 ) {
            continue;
         }
         Task* t = nullptr;
         {
            std::lock_guard<std::mutex> lk( s.m );
            if ( !s.q.empty() ) {
               t = s.q.front();
               s.q.pop_front();
               s.count.store( int( s.q.size() ) );
            }
         }
         if ( t ) {
            --offered_;
            run( t );
            return true;
         }
      }
      return false;
   }

   void worker( int slot ) {
      t_slot = slot;
      unsigned seed = unsigned( slot ) * 2654435761u + 1u;
      int idle = 0;
      while ( !quit_.load() ) {
         if ( steal_and_run( slot, seed ) ) {
            idle = 0;
            continue;
         }
         if ( ++idle < 100 ) {
            std::this_thread::yield();
            continue;
         }
         std::unique_lock<std::mutex> lk( sleep_m_ );
         ++sleepers_;
         sleep_cv_.wait( lk, [&] { return quit_.load() || offered_.load() > 0; } );
         --sleepers_;
         idle = 0;
      }
   }

   void start( int n ) {
      n_ = std::max( 1, n );
      quit_ = false;
      slots_.clear();
      for ( int i = 0; i < n_; ++i ) {
         slots_.push_back( std::make_unique<Slot>() );
      }
      for ( int i = 1; i < n_; ++i ) {
         threads_.emplace_back( &Stealer::worker, this, i );
      }
   }

   void stop() {
      {
         std::lock_guard<std::mutex> lk( sleep_m_ );
         quit_ = true;
      }
      sleep_cv_.notify_all();
      for ( std::thread& t : threads_ ) {
         t.join();
      }
      threads_.clear();
   }

   int n_ = 1;
   std::vector<std::unique_ptr<Slot>> slots_;
   std::vector<std::thread> threads_;
   std::atomic<bool> quit_{ false };
   std::atomic<long> offered_{ 0 };
   std::atomic<int> sleepers_{ 0 };
   std::mutex sleep_m_;
   std::condition_variable sleep_cv_;
};

Stealer& stealer() {
   static Stealer s;
   return s;
}

}  // namespace

int threads() { return pool().size(); }

void set_threads( int n ) {
   pool().resize( n < 1 ? hardware_threads() : n );
   stealer().resize( threads() );
}

bool in_parallel() { return t_in_parallel; }

void run_on_threads( const std::function<void( int )>& job ) { pool().run( job ); }

void fork_join( const std::function<void()>& a, const std::function<void()>& b ) { stealer().fork_join( a, b ); }

}  // namespace bh
