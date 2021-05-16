//
// Created by FeuFeve on 16/05/2021.
//

#ifndef RAYTRACER_CHRONO_H
#define RAYTRACER_CHRONO_H

#include <chrono>
#include <iomanip>
#include <iostream>
#include <utility>

using namespace std;
using namespace chrono;


class Chrono {

private:
    string name;

    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;
    double duration = 0;

    bool paused = false;
    high_resolution_clock::time_point pauseStart;
    double totalPauseDuration = 0;

    static double calculateDurationBetween(high_resolution_clock::time_point start, high_resolution_clock::time_point end) {
        double duration = duration_cast<nanoseconds>(end - start).count();
        duration *= 1e-9;
        return duration;
    }

public:
    explicit Chrono(const string& name) {
        this->name = name;
        start = high_resolution_clock::now();
    };

    bool pause() {
        if (paused)
            return false;

        paused = true;
        pauseStart = high_resolution_clock::now();
        return true;
    }

    bool unpause() {
        if (!paused)
            return false;

        paused = false;
        auto pauseEnd = high_resolution_clock::now();
        totalPauseDuration += calculateDurationBetween(pauseStart, pauseEnd);
        return true;
    }

    double getCurrentDuration() {
        double pauseDuration = 0;
        if (paused) {
            auto pauseEnd = high_resolution_clock::now();
            pauseDuration = calculateDurationBetween(pauseStart, pauseEnd);
        }

        end = high_resolution_clock::now();
        duration = calculateDurationBetween(start, end) - (totalPauseDuration + pauseDuration);
        return duration;
    }

    double getLastDuration() const {
        return duration - totalPauseDuration;
    }

    void displayCurrentDuration() {
        printDuration(name + ": ", getCurrentDuration());
    }

    void displayLastDuration() {
        printDuration(name + ": ", duration);
    }

    double getETR(double currentPercentage) {
        getCurrentDuration();
        return (100/currentPercentage) * duration - duration;
    }

    static void printDuration(const string& message, double duration) {
        cout << "# " << message << fixed << duration << setprecision(9) << "s" << endl;
    }
};


#endif //RAYTRACER_CHRONO_H
