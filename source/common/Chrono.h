//
// Created by FeuFeve on 16/05/2021.
//

#ifndef RAYTRACER_CHRONO_H
#define RAYTRACER_CHRONO_H

#include <chrono>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace chrono;


class Chrono {

private:
    high_resolution_clock::time_point start;
    high_resolution_clock::time_point end;

public:
    Chrono() {
        start = high_resolution_clock::now();
    };

    ~Chrono() = default;

    void stop(const string& message) {
        end = high_resolution_clock::now();
        double time_taken = duration_cast<nanoseconds>(end - start).count();
        time_taken *= 1e-9;

        cout << message << ": " << fixed << time_taken << setprecision(9) << "s" << endl;
    }
};


#endif //RAYTRACER_CHRONO_H
