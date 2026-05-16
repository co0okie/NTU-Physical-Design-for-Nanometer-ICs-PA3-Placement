#include "GlobalPlacer.h"

#include <cstdio>
#include <vector>
#include <set>
#include <cstdlib>
#include <numeric>

#include "ObjectiveFunction.h"
#include "Optimizer.h"
#include "Point.h"

template <typename T>
void printDistribution(const std::vector<Point2<T>>& points);

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
        int width = _placement.rectangleChip().width();
        int height = _placement.rectangleChip().height();
        
        positions[i].x = 0.05 * ((double)rand() / RAND_MAX - 0.5) * width + _placement.rectangleChip().centerX();
        positions[i].y = 0.05 * ((double)rand() / RAND_MAX - 0.5) * height + _placement.rectangleChip().centerY();
        _placement.module(i).setPosition(positions[i].x, positions[i].y);
    }
    
    ObjectiveFunction cost_function(_placement);

    SimpleConjugateGradient optimizer(cost_function, positions);

    optimizer.Initialize();

    cost_function(positions);
    cost_function.Backward(positions);

    printf("initial: wirelength %8.0f, density %8.0f, lambda %.6f, overflow %.3f\n", cost_function.wirelengthCost(), cost_function.densityCost(), cost_function.lambda, cost_function.overflowRatio());
    // plotPlacementResult( "init.plt", 1 );

    for (size_t i = 0; i < 10000; ++i) {
        optimizer.Step();
        for (size_t i = 0; i < num_modules; i++) {
            if (!_placement.module(i).isFixed()) 
                _placement.module(i).setPosition(positions[i].x, positions[i].y);
            else
                positions[i] = {_placement.module(i).x(), _placement.module(i).y()};
        }
        cost_function.lambda *= 0.999;
        if (i % 25 == 0) {
            double overflow_ratio = cost_function.overflowRatio();
            printf("iter %4lu, wirelength %8.0f, density %8.0f, lambda %.6f, overflow %.3f\n", i, cost_function.wirelengthCost(), cost_function.densityCost(), cost_function.lambda, overflow_ratio);
            // plotPlacementResult( "init.plt", 1 );
            // printDistribution(positions);
            // printDistribution(cost_function.grad());
            if (overflow_ratio < 0.1) break;
        }
        if (i % 500 == 0) {
            // plotPlacementResult( "init.plt", 1 );
        }
    }
    
    // plotPlacementResult( "init.plt", 1 );

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

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>

template <typename T>
void printDistribution(const std::vector<Point2<T>>& points) {
    if (points.empty()) {
        std::cout << "Vector is empty.\n";
        return;
    }

    // 1. 找出實際的 X 與 Y 範圍 (Min & Max)
    T min_x = points[0].x, max_x = points[0].x;
    T min_y = points[0].y, max_y = points[0].y;

    for (const auto& p : points) {
        if (p.x < min_x) min_x = p.x;
        if (p.x > max_x) max_x = p.x;
        if (p.y < min_y) min_y = p.y;
        if (p.y > max_y) max_y = p.y;
    }

    constexpr int n_div = 11;

    // 2. 切成 7 等分，計算每個區塊的寬度
    double width_x = static_cast<double>(max_x - min_x) / n_div;
    double width_y = static_cast<double>(max_y - min_y) / n_div;

    // 防止最大最小值一樣導致除以零
    if (width_x == 0) width_x = 1.0; 
    if (width_y == 0) width_y = 1.0;

    // 3. 計算每個點落在哪個格子內 (7x7)
    int counts[n_div][n_div] = {0};
    for (const auto& p : points) {
        int idx_x = static_cast<int>(std::floor((p.x - min_x) / width_x));
        int idx_y = static_cast<int>(std::floor((p.y - min_y) / width_y));
        
        // 確保剛好等於最大值時不會越界 (限制在 0~6)
        counts[std::clamp(idx_y, 0, n_div - 1)][std::clamp(idx_x, 0, n_div - 1)]++;
    }

    // 4. 印出分佈圖
    int total_points = points.size();
    const int col_width = 6; // 控制排版的欄寬

    // Y軸從上到下印出 (index 6 遞減到 0)
    for (int y = n_div; y --> 0;) {
        // 計算該區間的中間平均值
        double center_y = min_y + (y + 0.5) * width_y;
        
        // 印出 Y 軸標籤 (左側)
        printf("%6.0f", std::round(center_y));
        
        // 印出 7 個格子的千分比 (直接四捨五入)
        for (int x = 0; x < n_div; ++x) {
            printf("%6.0f", std::round(counts[y][x] * 1000.0 / total_points));
        }
        std::cout << "\n";
    }

    // 最後一行：印出 X 軸標籤 (底部)
    std::cout << std::setw(col_width) << " "; // 左下角對齊用的空格
    for (int x = 0; x < n_div; ++x) {
        double center_x = min_x + (x + 0.5) * width_x;
        printf("%6.0f", std::round(center_x));
    }
    std::cout << "\n";
}