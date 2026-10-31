#ifndef PGPHASE_TASK_POOL_HPP
#define PGPHASE_TASK_POOL_HPP

/**
 * @file task_pool.hpp
 * @brief Work the workers of one batch of chunks share.
 *
 * A chunk's independent pieces of work (its noisy regions' MSAs, collect_phase_noisy.cpp) go to the pool, and the
 * chunk's worker runs them itself too; a worker with no chunks left runs any chunk's until the batch is done
 * (run_graph_chunk_batch). Only the batch's workers ever run, so --threads holds. Without it the slowest chunk's
 * regions ran one after another while the other workers sat idle: chr20's 31 Mb chunk, ten MSAs, 104 s of a 131 s run.
 */

#include <condition_variable>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <mutex>

namespace pgphase_collect {

class TaskPool {
public:
    /// Runs fn(0) .. fn(n - 1) and returns once all have: this thread runs them, and idle workers (serve) take them
    /// too. This thread runs no other caller's, so a chunk resumes as soon as its own are done. Rethrows the first
    /// exception one threw.
    void run_all(std::size_t n, const std::function<void(std::size_t)>& fn) {
        if (n == 0) return;
        Group group;
        group.left = n;
        std::unique_lock<std::mutex> lock(mutex_);
        for (std::size_t i = 0; i < n; ++i) queue_.push_back(Task{&group, &fn, i});
        ready_.notify_all();
        while (group.left > 0) {
            if (!run_one(lock, &group)) ready_.wait(lock);
        }
        if (group.error) std::rethrow_exception(group.error);
    }

    /// Runs one queued task, if there is one: a worker between chunks takes the queue's work first, a chunk already
    /// under way waiting on it.
    bool help() {
        std::unique_lock<std::mutex> lock(mutex_);
        return run_one(lock, nullptr);
    }

    /// A worker with no chunks left: runs queued tasks until close().
    void serve() {
        std::unique_lock<std::mutex> lock(mutex_);
        while (true) {
            if (run_one(lock, nullptr)) continue;
            if (closed_) return;
            ready_.wait(lock);
        }
    }

    /// No task will come any more (every worker is out of chunks): serve() returns once the queue is empty.
    void close() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            closed_ = true;
        }
        ready_.notify_all();
    }

private:
    struct Group {
        std::size_t left = 0;
        std::exception_ptr error;
    };
    struct Task {
        Group* group;
        const std::function<void(std::size_t)>* fn;
        std::size_t i;
    };

    // Runs the first queued task of `only` (any task, when null) with the lock released; false when there is none.
    bool run_one(std::unique_lock<std::mutex>& lock, const Group* only) {
        auto it = queue_.begin();
        while (it != queue_.end() && only != nullptr && it->group != only) ++it;
        if (it == queue_.end()) return false;
        const Task task = *it;
        queue_.erase(it);
        lock.unlock();
        std::exception_ptr error;
        try {
            (*task.fn)(task.i);
        } catch (...) {
            error = std::current_exception();
        }
        lock.lock();
        if (error && !task.group->error) task.group->error = error;
        if (--task.group->left == 0) ready_.notify_all();
        return true;
    }

    std::mutex mutex_;
    std::condition_variable ready_;
    std::deque<Task> queue_;
    bool closed_ = false;
};

}  // namespace pgphase_collect

#endif
