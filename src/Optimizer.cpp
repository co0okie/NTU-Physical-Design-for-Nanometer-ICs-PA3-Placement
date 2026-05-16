#include "Optimizer.h"

#include <cmath>

SimpleConjugateGradient::SimpleConjugateGradient(BaseFunction &obj,
                                                 std::vector<Point2<double>> &var)
: BaseOptimizer(obj, var),
  grad_prev_(var.size()),
  dir_prev_(var.size()),
  step_(0) {
    fill(dir_prev_.begin(), dir_prev_.end(), Point2<double>(0, 0));
}

void SimpleConjugateGradient::Initialize() {
    // Before the optimization starts, we need to initialize the optimizer.
    step_ = 0;
}

/**
 * @details Update the solution once using the conjugate gradient method.
 */
void SimpleConjugateGradient::Step() {
    const size_t &kNumModule = var_.size();

    // Compute the gradient direction
    obj_(var_);       // Forward, compute the function value and cache from the input
    obj_.Backward(var_);  // Backward, compute the gradient according to the cache

    for (size_t i = 0; i < kNumModule; ++i) {
        Point2<double> dir = -obj_.grad().at(i) + 0.9 * dir_prev_.at(i);
        var_[i] = var_[i] + dir;
        dir_prev_[i] = dir;
    }

    step_++;
}
