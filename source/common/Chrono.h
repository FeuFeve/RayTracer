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
    double duration;

public:
    explicit Chrono(const string& name) {
        this->name = name;
        start = high_resolution_clock::now();
//        cout << "# " << name << " started." << endl;
    };

    ~Chrono() = default;

    Chrono* stop() {
        end = high_resolution_clock::now();
        duration = duration_cast<nanoseconds>(end - start).count();
        duration *= 1e-9;
        return this;
    }

    void display() {
        printDuration(name + " ended: ", duration);
    }

    double getDuration() const {
        return duration;
    }

    static void printDuration(const string& message, double duration) {
        cout << message << fixed << duration << setprecision(9) << "s" << endl;
    }
};


#endif //RAYTRACER_CHRONO_H
