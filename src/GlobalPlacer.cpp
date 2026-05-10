#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>
#include <cstdlib>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

GlobalPlacer::GlobalPlacer(Placement &placement)
    : _placement(placement) {
}

void GlobalPlacer::place() {

    ////////////////////////////////////////////////////////////////////
    // Global placement algorithm
    ////////////////////////////////////////////////////////////////////

    // TODO: Implement your global placement algorithm here.
    const size_t num_modules = _placement.numModules();  // You may modify this line.
    std::vector<Point2<double>> positions(num_modules);  // Optimization variables (positions of modules). You may modify this line.
    for (unsigned i = 0; i < num_modules; i++) {
        if (_placement.module(i).isFixed()) continue;
        int left = _placement.boundryLeft();
        int width = _placement.boundryRight() - left;
        int bottom = _placement.boundryBottom();
        int height = _placement.boundryTop() - bottom;
        positions[i].x = ((double)rand() / RAND_MAX - 0.5) * width * 0.1 + left + 0.5 * width;
        positions[i].y = ((double)rand() / RAND_MAX - 0.5) * height * 0.1 + bottom + 0.5 * height;
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    
    double gamma = 1;
    // Wirelength wirelength(_placement, gamma);
    Density density(_placement);

    const double kAlpha = 0.001;                         // Constant step size
    SimpleConjugateGradient optimizer(density, positions, kAlpha);

    optimizer.Initialize();
    
    auto minmax_x = minmax_element(positions.begin(), positions.end(), [](Point2<double> a, Point2<double> b) {
        return a.x < b.x;
    });
    auto minmax_y = minmax_element(positions.begin(), positions.end(), [](Point2<double> a, Point2<double> b) {
        return a.y < b.y;
    });
    printf("initial: f = %9.4f, min = (%f, %f), max = (%f, %f)\n", density(positions), minmax_x.first->x, minmax_y.first->y, minmax_x.second->x, minmax_y.second->y);

    for (size_t i = 0; i < 1000; ++i) {
        optimizer.Step();
        if (i % 1 == 0) {
            auto minmax_x = minmax_element(positions.begin(), positions.end(), [](Point2<double> a, Point2<double> b) {
                return a.x < b.x;
            });
            auto minmax_y = minmax_element(positions.begin(), positions.end(), [](Point2<double> a, Point2<double> b) {
                return a.y < b.y;
            });
            printf("iter = %3lu, f = %9.4f, min = (%f, %f), max = (%f, %f)\n", i, density.value(), minmax_x.first->x, minmax_y.first->y, minmax_x.second->x, minmax_y.second->y);
        }
        for (size_t i = 0; i < num_modules; i++) {
            if (_placement.module(i).isFixed()) continue;
            _placement.module(i).setPosition(positions[i].x, positions[i].y);
        }
    }


    ////////////////////////////////////////////////////////////////////
    // Write the placement result into the database. (You may modify this part.)
    ////////////////////////////////////////////////////////////////////
    size_t fixed_cnt = 0;
    for (size_t i = 0; i < num_modules; i++) {
        // If the module is fixed, its position should not be changed.
        // In this programing assignment, a fixed module may be a terminal or a pre-placed module.
        if (_placement.module(i).isFixed()) {
            fixed_cnt++;
            continue;
        }
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    printf("INFO: %lu / %lu modules are fixed.\n", fixed_cnt, num_modules);
}

void GlobalPlacer::plotPlacementResult(const string outfilename, bool isPrompt) {
    ofstream outfile(outfilename.c_str(), ios::out);
    outfile << " " << endl;
    outfile << "set title \"wirelength = " << _placement.computeHpwl() << "\"" << endl;
    outfile << "set size ratio 1" << endl;
    outfile << "set nokey" << endl
            << endl;
    outfile << "plot[:][:] '-' w l lt 3 lw 2, '-' w l lt 1" << endl
            << endl;
    outfile << "# bounding box" << endl;
    plotBoxPLT(outfile, _placement.boundryLeft(), _placement.boundryBottom(), _placement.boundryRight(), _placement.boundryTop());
    outfile << "EOF" << endl;
    outfile << "# modules" << endl
            << "0.00, 0.00" << endl
            << endl;
    for (size_t i = 0; i < _placement.numModules(); ++i) {
        Module &module = _placement.module(i);
        plotBoxPLT(outfile, module.x(), module.y(), module.x() + module.width(), module.y() + module.height());
    }
    outfile << "EOF" << endl;
    outfile << "pause -1 'Press any key to close.'" << endl;
    outfile.close();

    if (isPrompt) {
        char cmd[200];
        sprintf(cmd, "gnuplot %s", outfilename.c_str());
        if (!system(cmd)) {
            cout << "Fail to execute: \"" << cmd << "\"." << endl;
        }
    }
}

void GlobalPlacer::plotBoxPLT(ofstream &stream, double x1, double y1, double x2, double y2) {
    stream << x1 << ", " << y1 << endl
           << x2 << ", " << y1 << endl
           << x2 << ", " << y2 << endl
           << x1 << ", " << y2 << endl
           << x1 << ", " << y1 << endl
           << endl;
}
