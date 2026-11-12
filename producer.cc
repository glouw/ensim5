#include "ensim.hh"

#include <queue>
#include <array>
#include <atomic>
#include <thread>
#include <mutex>

namespace ensim {

static constexpr size_t g_producer_block_size = 200;
static constexpr size_t g_producer_queue_size = 4096;
static constexpr size_t g_producer_history_size = 4096;
static constexpr int g_producer_poll_us = 200;

template<typename T, size_t N>
class ring
{
    std::array<T, N> self = {};
    size_t head = 0;
    size_t tail = 0;
    size_t elems = 0;

public:
    void clear()
    {
        head = tail = elems = 0;
    }

    size_t size() const
    {
        return elems;
    }

    size_t capacity() const
    {
        return N;
    }

    void push_back(const T& value)
    {
        self[tail++] = value;
        tail %= N;
        if(elems == N)
        {
            head++;
            head %= N;
        }
        else
        {
            elems++;
        }
    }

    bool empty() const
    {
        return elems == 0;
    }

    void pop_front()
    {
        if(not empty())
        {
            head++;
            head %= N;
            elems--;
        }
    }

    size_t index(const size_t i) const
    {
        return (head + i) % N;
    }

    T& operator[](const size_t i)
    {
        return self[index(i)];
    }

    const T& operator[](const size_t i) const
    {
        return self[index(i)];
    }
};

class producer
{
    ensim::engine* engine = nullptr;
    ring<float, g_producer_queue_size> queue = {};
    ring<size_t, g_producer_history_size> history = {};
    std::vector<size_t> snapshot = {};
    std::mutex engine_mutex = {};
    std::mutex queue_mutex = {};
    std::mutex history_mutex = {};
    std::atomic<bool> done = false;
    std::thread thread = std::thread(&producer::run, this);
    std::vector<float> ready = {};

public:
    producer()
    {
        snapshot.reserve(g_producer_history_size);
    }

    std::vector<float>& consume(const size_t samples)
    {
        std::lock_guard lock1(queue_mutex);
        ready.clear();
        for(size_t i = 0; i < samples; i++)
        {
            ready.push_back(queue[0]);
            queue.pop_front();
        }
        std::lock_guard lock2(history_mutex);
        history.push_back(queue.size());
        return ready;
    }

    bool must_produce()
    {
        std::lock_guard lock(queue_mutex);
        return queue.size() + g_producer_block_size < queue.capacity();
    }

    bool produce()
    {
        if(must_produce())
        {
            std::lock_guard lock1(engine_mutex);
            if(engine)
            {
                engine->run(g_producer_block_size);
                const std::vector<float>& slice = engine->get_audio_signal();
                std::lock_guard lock2(queue_mutex);
                for(float value: slice)
                {
                    queue.push_back(value);
                }
                std::lock_guard lock3(history_mutex);
                history.push_back(queue.size());
                return true;
            }
        }
        return false;
    }

    void run()
    {
        while(not done)
        {
            if(not produce())
            {
                const auto delay = std::chrono::microseconds(g_producer_poll_us);
                std::this_thread::sleep_for(delay);
            }
        }
    }

    void set(ensim::engine* other)
    {
        std::lock_guard lock(engine_mutex);
        engine = other;
    }

    void kill()
    {
        done = true;
        thread.join();
    }

    const std::vector<size_t>& get_history()
    {
        std::lock_guard lock1(history_mutex);
        snapshot.clear();
        for(size_t i = 0; i < history.size(); i++)
        {
            snapshot.push_back(history[i]);
        }
        return snapshot;
    }
};

producer g_producer = {};

const std::vector<float>& consume(const size_t frames)
{
    return g_producer.consume(frames);
}

const std::vector<size_t>& get_buffer_size_history()
{
    return g_producer.get_history();
}

void attach(ensim::engine* engine)
{
    g_producer.set(engine);
}

void kill()
{
    g_producer.kill();
}

}
