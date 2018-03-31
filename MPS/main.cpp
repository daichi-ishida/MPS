#include <iostream>
#include <sstream>
#include <iomanip>
#include "Scene.h"
#include "Simulator.h"
#include "constants.h"


int main() {
    double time = 0;
    // ファイル出力インターバル用カウンタ
    int interval_counter = 0;

    Scene scene(time);

    // シーンから粒子配列作成
    scene.makeParticleSet();
    
    Simulator simulator(scene.getParticleSetPtr());

    std::cout << "\n*** START PARTICLE-SIMULATION ***\n";

    scene.writeData();

    while (1)
    {
        ++interval_counter;
        time += DELTA_TIME;
        simulator.update();

        if (interval_counter % OUTPUT_INTERVAL == 0)
        {
            std::ostringstream sout;
            sout << std::setfill('0') << std::setw(4) << std::right << interval_counter;
            std::cout << "TimeStepNumber: " << sout.str();
            sout.str("");
            sout.clear();
            sout << std::fixed << std::setprecision(3) << time;
            std::cout << "   Time: " << sout.str() << "(s)   NumberOfParticles: " << simulator.getNumberOfParticles() << std::endl;

            scene.writeData();
        }
        if (time >= FINISH_TIME) { break; }
    }

    std::cout << "*** END ***\n\n";
    return 0;
}