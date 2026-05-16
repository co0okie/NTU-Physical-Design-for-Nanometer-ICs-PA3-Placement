#define _GLIBCXX_USE_CXX11_ABI 0  // Align the ABI version to avoid compatibility issues with `Placment.h`
#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H

#include <vector>

#include "Placement.h"
#include "Point.h"

/**
 * @brief Base class for objective functions
 */
class BaseFunction {
   public:
    /////////////////////////////////
    // Conssutructors
    /////////////////////////////////

    BaseFunction(const size_t &input_size) : grad_(input_size) {}

    /////////////////////////////////
    // Accessors
    /////////////////////////////////

    const std::vector<Point2<double>> &grad() const { return grad_; }
    const double &value() const { return value_; }

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    // Forward pass, compute the value of the function
    virtual const double &operator()(const std::vector<Point2<double>> &input = {}) = 0;

    // Backward pass, compute the gradient of the function
    virtual const std::vector<Point2<double>> &Backward(const std::vector<Point2<double>> &input = {}) = 0;

   protected:
    /////////////////////////////////
    // Data members
    /////////////////////////////////

    std::vector<Point2<double>> grad_;  // Gradient of the function
    double value_;                      // Value of the function
};

/**
 * @brief Wirelength function
 */
class Wirelength : public BaseFunction {
    // TODO: Implement the wirelength function, add necessary data members for caching
   public:
    Wirelength(Placement &placement, double gamma);

    /////////////////////////////////
    // Methods
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input = {}) override;
    const std::vector<Point2<double>> &Backward(const std::vector<Point2<double>> &input = {}) override;

    double gamma; // the smaller gamma is, the more accurate to HPWL

   private:
    Placement &placement_;
    vector<Point2<double>> max_of_net;
    vector<Point2<double>> min_of_net;
    vector<Point2<double>> max_wa_of_net;
    vector<Point2<double>> min_wa_of_net;
    vector<Point2<double>> denominator_max_of_net;
    vector<Point2<double>> denominator_min_of_net;
};

class DensitySigmoid : public BaseFunction {
   public:
    DensitySigmoid(Placement &placement);

    const double &operator()(const std::vector<Point2<double>> &input = {}) override;
    const std::vector<Point2<double>> &Backward(const std::vector<Point2<double>> &input = {}) override;

    Point2<double> binDim() const { return dim_bin; }

    double overflowRatio() const;

   private:
    Placement &placement_;
    vector<vector<double>> bins; // total density of a bin
    Point2<double> dim_bin; // x: width, y: height

    // get [col, row] from point p
    Point2<int> binIndex(Point2<double> p) {
        Point2<double> chip(placement_.boundryLeft(), placement_.boundryBottom());
        return Floor((p - chip) / dim_bin);
    }
    Point2<double> binCenter(unsigned col, unsigned row) {
        Point2<double> chip(placement_.boundryLeft(), placement_.boundryBottom());
        return chip + (0.5 + Point2<double>(col, row)) * dim_bin;
    }
    pair<Point2<int>, Point2<int>> getBeginEndBinIndex(Module& module);
};

/**
 * @brief Objective function for global placement
 */
class ObjectiveFunction : public BaseFunction {
    // TODO: Implement the objective function for global placement, add necessary data
    // members for caching
    //
    // Hint: The objetive function of global placement is as follows:
    //       f(t) = wirelength(t) + lambda * density(t),
    // where t is the positions of the modules, and lambda is the penalty weight.
    // You may need an interface to update the penalty weight (lambda) dynamically.
   public:
    ObjectiveFunction(Placement &placement);
    /////////////////////////////////
    // Methods
    /////////////////////////////////

    const double &operator()(const std::vector<Point2<double>> &input = {}) override;
    const std::vector<Point2<double>> &Backward(const std::vector<Point2<double>> &input = {}) override;

    const double wirelengthCost() const { return wirelength_.value(); }
    const double densityCost() const { return density_.value(); }
    Point2<double> binDim() const { return density_.binDim(); }
    double overflowRatio() const { return density_.overflowRatio(); }

    double lambda;

   private:
    Placement &placement_;
    Wirelength wirelength_;
    DensitySigmoid density_;
};

#endif  // OBJECTIVEFUNCTION_H
