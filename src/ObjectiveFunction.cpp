#include "ObjectiveFunction.h"

#include "cstdio"
#include "cfloat"
#include "numeric"

Wirelength::Wirelength(Placement& placement, double gamma) 
: BaseFunction(placement.numModules()), gamma(gamma), placement_(placement), 
  max_of_net(placement.numNets()), min_of_net(placement.numNets()),
  max_wa_of_net(placement.numNets()), min_wa_of_net(placement.numNets()),
  denominator_max_of_net(placement.numNets()), denominator_min_of_net(placement.numNets()) {
}

const double& Wirelength::operator()(const std::vector<Point2<double>>& input) {
    value_ = 0;
    for (unsigned i = 0; i < placement_.numNets(); i++) {
        Net& net = placement_.net(i);
        max_of_net[i] = -DBL_MAX;
        min_of_net[i] = DBL_MAX;
        Point2<double>& max_ = max_of_net[i];
        Point2<double>& min_ = min_of_net[i];
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            max_ = Max(max_, {pin.x(), pin.y()});
            min_ = Min(min_, {pin.x(), pin.y()});
        }
        Point2<double> gamma = (max_of_net[i] - min_of_net[i]) / 5;
        Point2<double> cache;
        Point2<double> num_max(0);
        Point2<double> num_min(0);
        Point2<double>& den_max = denominator_max_of_net[i] = 0;
        Point2<double>& den_min = denominator_min_of_net[i] = 0;
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            Point2<double> p = {pin.x(), pin.y()};
            num_max += p * (cache = Exp((p - max_) / gamma));
            den_max += cache;
            num_min += p * (cache = Exp((-p - -min_) / gamma));
            den_min += cache;
        }
        max_wa_of_net[i] = num_max / den_max;
        min_wa_of_net[i] = num_min / den_min;
        value_ += max_wa_of_net[i].x - min_wa_of_net[i].x +
                  max_wa_of_net[i].y - min_wa_of_net[i].y;
    }
    return value_;
}

const std::vector<Point2<double>>& Wirelength::Backward(
const std::vector<Point2<double>>& input) {
    fill(grad_.begin(), grad_.end(), Point2<double>(0, 0));
    for (unsigned i = 0; i < placement_.numNets(); i++) {
        Net& net = placement_.net(i);
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            Point2<double> p = {pin.x(), pin.y()};
            Point2<double> gamma = (max_of_net[i] - min_of_net[i]) / 5;
            grad_[pin.moduleId()] +=
                Exp((p - max_of_net[i]) / gamma) / denominator_max_of_net[i] *
                (1 + (p - max_wa_of_net[i]) / gamma)
                -
                Exp((-p - -min_of_net[i]) / gamma) / denominator_min_of_net[i] *
                (1 - (p - min_wa_of_net[i]) / gamma);
        }
    }
    return grad_;
}

DensitySigmoid::DensitySigmoid(Placement& placement)
: BaseFunction(placement.numModules()), placement_(placement){
    // creating the bin grid
    vector<unsigned> modules(placement.numModules());
    Point2<double> module_dim_min = {DBL_MAX, DBL_MAX};
    Point2<double> module_dim_avg = {0, 0};
    for (unsigned i = 0; i < placement.numModules(); i++) {
        Module& module = placement.module(i);
        if (module.isFixed()) continue;
        if (module.width() != 0) module_dim_min.x = min(module_dim_min.x, module.width());
        if (module.height() != 0) module_dim_min.y = min(module_dim_min.y, module.height());
        module_dim_avg += {module.width(), module.height()};
    }
    module_dim_avg /= placement.numModules();
    printf("avg module width: %f, height: %f\n", module_dim_avg.x, module_dim_avg.y);
    // printf("width_min: %f, height_min: %f\n", module_dim_min.x, module_dim_min.y);

    double chip_width = placement.boundryRight() - placement.boundryLeft();
    double chip_height = placement.boundryTop() - placement.boundryBottom();
    unsigned column_count = max(1, (int)round(chip_width / (2 * module_dim_avg.x)));
    unsigned row_count = max(1, (int)round(chip_height / (2 * module_dim_avg.y)));
    bins.assign(column_count, vector<double>(row_count, 0));
    dim_bin.x = chip_width / column_count;
    dim_bin.y = chip_height / row_count;

    printf("total bins: %d * %d = %d\n", column_count, row_count, column_count * row_count);
}

pair<Point2<int>, Point2<int>> DensitySigmoid::getBeginEndBinIndex(Module& module) {
    Point2<int> n_bins(bins.size(), bins[0].size());
    auto [dim_min, dim_max] = MinMax({module.width(), module.height()}, dim_bin);
    auto d_max = 0.5 * dim_max + 1.5 * dim_min;
    Point2<double> m(module.centerX(), module.centerY());
    Point2<int> begin_index = Clamp(binIndex(m - d_max), {0, 0}, n_bins - 1);
    Point2<int> end_index = Clamp(binIndex(m + d_max), {0, 0}, n_bins - 1);
    return make_pair(begin_index, end_index);
}

double sigmoid(double x) {
    if (x < -100.0) return 0.0;
    if (x > 100.0) return 1.0;
    return 1. / (1. + exp(-4 * x));
}

/// @brief 1D density model using sigmoid function
/// @param x difference between position of module and bin center
/// @param wm width/height of module
/// @param wb width/height of bin
/// @return 
double sigmoidDensity(double x, double wm, double wb) {
    auto [w_min, w_max] = minmax(wm, wb);
    double left = sigmoid((x + 0.5 * w_max) / w_min);
    double right = sigmoid((x - 0.5 * w_max) / w_min);
    return w_min / wb * (left - right);
}

const double& DensitySigmoid::operator()(const std::vector<Point2<double>>& input) {
    value_ = 0;
    for (auto& col : bins) fill(col.begin(), col.end(), 0);
    int n_cols = bins.size();
    int n_rows = n_cols ? bins[0].size() : 0;
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        Module& module = placement_.module(i);

        Point2<double> m(module.centerX(), module.centerY());
        Point2<double> dim_m(module.width(), module.height());

        auto [begin_index, end_index] = getBeginEndBinIndex(module);

        for (int col = begin_index.x; col <= end_index.x; ++col) {
            for (int row = begin_index.y; row <= end_index.y; ++row) {
                auto bin_center = binCenter(col, row);
                Point2<double> d = m - bin_center;

                bins[col][row] += sigmoidDensity(d.x, dim_m.x, dim_bin.x)
                                * sigmoidDensity(d.y, dim_m.y, dim_bin.y);
            }
        }
    }
    for (int col = 0; col < n_cols; ++col) {
        for (int row = 0; row < n_rows; ++row) {
            double exceed_density = max(0., bins[col][row] - 1);
            value_ += exceed_density * exceed_density;
        }
    }
    return value_ *= 0.5;
}

double diffSigmoid(double x) {
    double y = sigmoid(x);
    return 4 * y * (1 - y);
}

double diffSigmoidDensity(double x, double wm, double wb) {
    auto [w_min, w_max] = minmax(wm, wb);
    double left = diffSigmoid((x + 0.5 * w_max) / w_min);
    double right = diffSigmoid((x - 0.5 * w_max) / w_min);
    return 1 / wb * (left - right);
}

const std::vector<Point2<double>>& DensitySigmoid::Backward(const std::vector<Point2<double>>& input) {
    fill(grad_.begin(), grad_.end(), Point2<double>(0, 0));
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        Module& module = placement_.module(i);
        if (module.isFixed()) continue;

        Point2<double> m(module.centerX(), module.centerY());
        Point2<double> dim_m(module.width(), module.height());

        auto [begin_index, end_index] = getBeginEndBinIndex(module);
        for (int col = begin_index.x; col <= end_index.x; ++col) {
            for (int row = begin_index.y; row <= end_index.y; ++row) {
                auto bin_center = binCenter(col, row);
                Point2<double> d = m - bin_center;
                grad_[i] += {
                    max(0., bins[col][row] - 1)
                        * diffSigmoidDensity(d.x, dim_m.x, dim_bin.x)
                        * sigmoidDensity(d.y, dim_m.y, dim_bin.y),
                    max(0., bins[col][row] - 1)
                        * sigmoidDensity(d.x, dim_m.x, dim_bin.x)
                        * diffSigmoidDensity(d.y, dim_m.y, dim_bin.y)
                };
            }
        }
        grad_[i] *= dim_bin.x * dim_bin.y;
        // double density_norm = sqrt(grad_[i].x * grad_[i].x + grad_[i].y * grad_[i].y);
        // if (density_norm > 2.) grad_[i] *= 2. / density_norm;
    }
    return grad_;
}

double DensitySigmoid::overflowRatio() const {
    return accumulate(bins.begin(), bins.end(), 0., [&](const double a, const vector<double>& b) {
        return a + accumulate(b.begin(), b.end(), 0., [&](const double a, const double b) {
            return a + max(0., b - 1);
        });
    }) / bins.size() / bins[0].size();
}

ObjectiveFunction::ObjectiveFunction(Placement &placement) 
: BaseFunction(placement.numModules()), placement_(placement), wirelength_(placement, 1), density_(placement) {
    wirelength_();
    wirelength_.Backward();
    density_();
    density_.Backward();
    lambda = 1;
}

const double& ObjectiveFunction::operator()(const std::vector<Point2<double>>& input) {
    wirelength_(input);
    density_(input);
    value_ = 0;
    return value_;
}

const std::vector<Point2<double>>& ObjectiveFunction::Backward(const std::vector<Point2<double>>& input) {
    wirelength_.Backward(input);
    density_.Backward(input);
    Point2<double> chip(placement_.boundryRight() - placement_.boundryLeft(), 
        placement_.boundryTop() - placement_.boundryBottom());
    double step = sqrt(density_.binDim().x * density_.binDim().y) * 0.6;
    double max_density_grad_norm = -DBL_MAX;
    double max_wirelength_grad_norm = -DBL_MAX;
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        if (placement_.module(i).isFixed()) continue;
        auto gw = wirelength_.grad()[i];
        auto gd = density_.grad()[i];
        max_wirelength_grad_norm = max(max_wirelength_grad_norm, gw.x * gw.x + gw.y * gw.y);
        max_density_grad_norm = max(max_density_grad_norm, gd.x * gd.x + gd.y * gd.y);
    }
    max_wirelength_grad_norm = sqrt(max_wirelength_grad_norm);
    max_density_grad_norm = sqrt(max_density_grad_norm);
    double wirelength_scale = max_wirelength_grad_norm > 1e-12 ? 
        lambda / max_wirelength_grad_norm : 0;
    double density_scale = max_density_grad_norm > 1e-12 ? 
        (1 - lambda) / max_density_grad_norm : 0;
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        Module& module = placement_.module(i);
        if (module.isFixed()) {
            grad_[i] = {0, 0};
            continue;
        }
        Point2<double> out_of_bound_cost = {
            min(0., module.centerX() - placement_.boundryLeft() - 0.5 * density_.binDim().x) +
            max(0., module.centerX() - placement_.boundryRight() + 0.5 * density_.binDim().x),
            min(0., module.centerY() - placement_.boundryBottom() - 0.5 * density_.binDim().y) +
            max(0., module.centerY() - placement_.boundryTop() + 0.5 * density_.binDim().y)
        };
        grad_[i] = (
            wirelength_.grad()[i] * wirelength_scale +
            density_.grad()[i] * density_scale
        ) * step + 0.5 * out_of_bound_cost;
    }
    return grad_;
}
