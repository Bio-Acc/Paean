#pragma once

#include <iostream>
#include <string>
#include <chrono>
#include <utility>

using namespace std::chrono;

class Clock {
public:
    Clock(std::string event)
        : _event(std::move(event)) {}

    Clock() = delete;

    inline void start() {
        _ts = high_resolution_clock::now();
    }
    inline void stop() {
        _te = high_resolution_clock::now();
    }

    void log() {
        double sec = duration_cast<milliseconds>(_te - _ts).count();
        std::cout << _event << " costs " << (sec / 1000) << "s" << std::endl;
    }

private:
    std::string _event;
    high_resolution_clock::time_point _ts;
    high_resolution_clock::time_point _te;
};

#define PROFILE_START(event)  \
    Clock clk##event(#event); \
    clk##event.start();

#define PROFILE_STOP(event) \
    clk##event.stop();      \
    clk##event.log();
