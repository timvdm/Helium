/*
 * Copyright (c) 2013, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HELIUM_CONCURRENT_H
#define HELIUM_CONCURRENT_H

#include <Helium/contract.h>

#include <vector>
#include <thread>

namespace Helium {

  /**
   * @brief Class to run taks in parallel.
   *
   * The Concurrent class can be used to execute tasks that can be put in a
   * std::vector. The results will also be placed in a std::vector of the same
   * size.
   *
   * @tparam Callable A functor or function to do the acutal work.
   * @tparam Task The type of an individual task.
   * @tparam Result The type of the result of an individual task.
   */
  template<typename Callable, typename Task, typename Result>
  class Concurrent
  {
    private:
      /**
       * @brief Wrapper to call callable.
       *
       * This wrapper calls the callable and sets the done flag once finished.
       *
       * @param callable The functor of function to do the actual work.
       * @param task An individual task item.
       * @param result Output for the result of the task item.
       * @param done The done flag, entry from m_done.
       */
      static void runTask(Callable callable, const Task &task, Result &result, char &done)
      {
        callable(task, result);
        done = 1;
      }

    public:
      /**
       * @pre @p tasks and @p results should have the same size (i.e. results
       *      should be resized before calling this function).
       *
       * @param callable The functor of function to do the actual work.
       * @param tasks The tasks to be run in parallel.
       * @param results Output variable to put the results into.
       */
      void run(Callable callable, const std::vector<Task> &tasks, std::vector<Result> &results)
      {
        PRE(tasks.size() == results.size());

        // determine number of threads
        unsigned numThreads = std::thread::hardware_concurrency();
        if (!numThreads)
          numThreads = 2;
        if (numThreads > tasks.size())
          numThreads = tasks.size();

        // keep track of the number of processed tasks
        std::size_t next = numThreads;
        std::size_t processed = 0;
        m_done.resize(numThreads);

        // create initial threads
        for (int i = 0 ; i < numThreads; ++i)
          m_threads.push_back(std::thread(runTask, callable, std::ref(tasks[i]), std::ref(results[i]), std::ref(m_done[i])));

        // process all
        while (true) {

          // check for threads that have finished
          for (int i = 0 ; i < numThreads; ++i) {
            if (m_done[i]) {
              processed++;
              m_done[i] = 0;
              m_threads[i].join();

              // launch a new thread if there are more tasks remaining
              if (next < tasks.size()) {
                m_threads[i] = std::thread(runTask, callable, std::ref(tasks[next]), std::ref(results[next]), std::ref(m_done[i]));
                ++next;
              }
            }
          }

          // stop if all tasks are processed
          if (processed == tasks.size())
            break;

          // yield this thread to let the worker threads do their work
          std::this_thread::yield();
        }
      }

    private:
      std::vector<std::thread> m_threads; //!< The threads.
      std::vector<char> m_done; //!< This flag is set when a thread is done processing it's task.
  };

}

#endif
