#include "task_pool.hpp"

#include <atomic>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace pgphase_collect;

static bool check(bool condition, const std::string& message) {
    if (!condition) std::cerr << "FAIL: " << message << "\n";
    return condition;
}

// The caller alone (no worker serving): every task runs, on the caller.
static bool check_caller_alone() {
    TaskPool pool;
    std::vector<int> out(100, 0);
    const std::thread::id caller = std::this_thread::get_id();
    std::atomic<bool> elsewhere{false};
    pool.run_all(out.size(), [&](size_t i) {
        out[i] = static_cast<int>(i) * 2;
        if (std::this_thread::get_id() != caller) elsewhere = true;
    });
    bool ok = true;
    for (size_t i = 0; i < out.size(); ++i) ok = check(out[i] == static_cast<int>(i) * 2, "task " + std::to_string(i)) && ok;
    ok = check(!elsewhere, "a task ran off the caller with no worker serving") && ok;
    pool.run_all(0, [](size_t) { throw std::runtime_error("no task to run"); });
    pool.close();
    pool.serve();  // returns at once: closed, nothing queued
    return ok;
}

// Two callers and three serving workers: each caller gets every one of its own tasks run exactly once; serve()
// returns after close().
static bool check_callers_and_workers() {
    TaskPool pool;
    std::vector<std::thread> workers;
    for (int w = 0; w < 3; ++w) workers.emplace_back([&pool] { pool.serve(); });
    std::vector<std::atomic<int>> a(500), b(300);
    std::thread other([&] { pool.run_all(b.size(), [&](size_t i) { ++b[i]; }); });
    pool.run_all(a.size(), [&](size_t i) { ++a[i]; });
    other.join();
    pool.close();
    for (std::thread& w : workers) w.join();
    bool ok = true;
    for (size_t i = 0; i < a.size(); ++i) ok = check(a[i] == 1, "caller a, task " + std::to_string(i)) && ok;
    for (size_t i = 0; i < b.size(); ++i) ok = check(b[i] == 1, "caller b, task " + std::to_string(i)) && ok;
    return ok;
}

// A task that throws: the others still run, and run_all rethrows once all have.
static bool check_exception() {
    TaskPool pool;
    std::thread worker([&pool] { pool.serve(); });
    std::atomic<int> ran{0};
    bool thrown = false;
    try {
        pool.run_all(50, [&](size_t i) {
            ++ran;
            if (i == 7) throw std::runtime_error("task 7");
        });
    } catch (const std::runtime_error& e) {
        thrown = std::string(e.what()) == "task 7";
    }
    pool.close();
    worker.join();
    return check(thrown, "run_all rethrows the task's exception") && check(ran == 50, "every task ran");
}

// help() runs one queued task and reports whether there was one.
static bool check_help() {
    TaskPool pool;
    bool ok = check(!pool.help(), "help() with nothing queued");
    std::atomic<bool> started{false}, release{false};
    std::atomic<int> ran{0};
    // The caller runs task 0, which waits until the main thread has helped with task 1.
    std::thread caller([&] {
        pool.run_all(2, [&](size_t i) {
            if (i == 0) {
                started = true;
                while (!release) std::this_thread::yield();
            }
            ++ran;
        });
    });
    while (!started) std::this_thread::yield();
    while (!pool.help()) std::this_thread::yield();
    release = true;
    caller.join();
    ok = check(ran == 2, "both tasks ran") && ok;
    return ok;
}

int main() {
    bool ok = check_caller_alone();
    ok = check_callers_and_workers() && ok;
    ok = check_exception() && ok;
    ok = check_help() && ok;
    if (ok) std::cout << "ALL PASS\n";
    return ok ? 0 : 1;
}
