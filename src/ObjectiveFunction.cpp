#include "ObjectiveFunction.h"

#include "cstdio"
#include "cfloat"
#include "numeric"

ExampleFunction::ExampleFunction(Placement &placement) : BaseFunction(1), placement_(placement)
{
    printf("Fetch the information you need from placement database.\n");
    printf("For example:\n");
    printf("    Placement boundary: (%.f,%.f)-(%.f,%.f)\n", placement_.boundryLeft(), placement_.boundryBottom(),
           placement_.boundryRight(), placement_.boundryTop());
}

const double &ExampleFunction::operator()(const std::vector<Point2<double>> &input)
{
    // Compute the value of the function
    value_ = 3. * input[0].x * input[0].x + 2. * input[0].x * input[0].y +
             2. * input[0].y * input[0].y + 7.;
    return value_;
}

const std::vector<Point2<double>> &ExampleFunction::Backward(const std::vector<Point2<double>> &input)
{
    // Compute the gradient of the function
    grad_[0].x = 6. * input[0].x + 2. * input[0].y;
    grad_[0].y = 2. * input[0].x + 4. * input[0].y;
    return grad_;
}

Wirelength::Wirelength(Placement& placement, double gamma) 
: BaseFunction(placement.numModules()), placement_(placement), gamma(gamma),
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

Density::Density(Placement& placement) 
: BaseFunction(placement.numModules()), placement_(placement) {

    double total_module_area = 0;
    for (unsigned i = 0; i < placement.numModules(); i++) {
        total_module_area += placement.module(i).area();
    }
    double chip_area = (placement.boundryRight() - placement.boundryLeft()) * 
        (placement.boundryTop() - placement.boundryBottom());
    target_density = total_module_area / chip_area;
    
    // creating the bin grid
    vector<unsigned> modules(placement.numModules());
    iota(modules.begin(), modules.end(), 0);
    sort(modules.begin(), modules.end(), [&](unsigned a, unsigned b) {
        return placement.module(a).width() < placement.module(b).width();
    });
    double width_medium = placement.numModules() > 0 ? 
        placement.module(placement.numModules() / 2).width() : 0;
    sort(modules.begin(), modules.end(), [&](unsigned a, unsigned b) {
        return placement.module(a).height() < placement.module(b).height();
    });
    double height_medium = placement.numModules() > 0 ? 
        placement.module(placement.numModules() / 2).height() : 0;
    double chip_width = placement.boundryRight() - placement.boundryLeft();
    double chip_height = placement.boundryTop() - placement.boundryBottom();
    unsigned column_count = max(1, (int)round(chip_width / (2 * width_medium)));
    unsigned row_count = max(1, (int)round(chip_height / (2 * height_medium)));
    bins.assign(column_count, vector<double>(row_count, 0));
    dim_bin.x = chip_width / column_count;
    dim_bin.y = chip_height / row_count;

    printf("total bins: %d * %d = %d\n", column_count, row_count, column_count * row_count);
    printf("target density: %f\n", target_density);
}

pair<Point2<int>, Point2<int>> Density::getBeginEndBinIndex(Module& module) {
    int num_cols = bins.size();
    int num_rows = num_cols ? bins[0].size() : 0;
    Point2<double> m(module.centerX(), module.centerY());
    Point2<double> dim_m(module.width(), module.height());
    Point2<double> boundry(placement_.boundryLeft(), placement_.boundryBottom());

    Point2<double> max_dist = (dim_m / 2.0) + 2.0 * dim_bin;
    Point2<double> begin = m - max_dist;
    Point2<double> end = m + max_dist;

    Point2<int> begin_index = Ceil((begin - boundry) / dim_bin - 0.5);
    Point2<int> end_index = Floor((end - boundry) / dim_bin - 0.5);

    begin_index = Max(begin_index, Point2<int>(0, 0));
    end_index = Min(end_index, Point2<int>(num_cols - 1, num_rows - 1));

    return make_pair(begin_index, end_index);
}

// return min(width(module), width(bin)) * p_x(bin, module)
// d = |centerX(module) - centerX(bin)|
// dim_m: width of module
// dim_bin: width of bin
double bell_shape(double d, double dim_m, double dim_bin) {
    double ret = dim_m * dim_bin / (2 / 3 * dim_m + 2 * dim_bin);
    double a = 4 / (dim_m + 2 * dim_bin) / (dim_m + 4 * dim_bin);
    double b = 2 / dim_bin / (dim_m + 4 * dim_bin);
    if (d < dim_m / 2 + dim_bin) {
        ret *= 1 - a * d * d;
    } else if (d < dim_m / 2 + 2 * dim_bin) {
        double t = (d - dim_m / 2 - 2 * dim_bin);
        ret *= b * t * t;
    } else {
        ret *= 0;
    }
    return ret;
}

// return ∂bell_shape / ∂d
// m: center of module
// bin: center of bin
// dim_m: width of module
// dim_m: width of bin
double grad_bell_shape(double m, double bin, double dim_m, double dim_bin) {
    double ret = dim_m * dim_bin / (2 / 3 * dim_m + 2 * dim_bin);
    double a = 4 / (dim_m + 2 * dim_bin) / (dim_m + 4 * dim_bin);
    double b = 2 / dim_bin / (dim_m + 4 * dim_bin);
    double sign = (bin < m) - (m < bin);
    double d = abs(m - bin);
    if (d < dim_m / 2 + dim_bin) {
        ret *= -2 * a * d * sign;
    } else if (d < dim_m / 2 + 2 * dim_bin) {
        ret *= 2 * b * (d - dim_m / 2 - 2 * dim_bin) * sign;
    } else {
        ret *= 0;
    }
    return ret;
}

const double& Density::operator()(const std::vector<Point2<double>>& input) {
    value_ = 0;
    for (auto& col : bins) fill(col.begin(), col.end(), 0);
    int num_cols = bins.size();
    int num_rows = num_cols ? bins[0].size() : 0;
    double target_overlap_area = target_density * dim_bin.x * dim_bin.y;
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        Module& module = placement_.module(i);

        Point2<double> m(module.centerX(), module.centerY());
        Point2<double> dim_m(module.width(), module.height());

        auto bin_index = getBeginEndBinIndex(module);

        // printf("module (%f, %f), dim: (%f, %f)", m.x, m.y, dim_m.x, dim_m.y);
        // printf(", bin index: (%d, %d) - (%d, %d)\n", bin_index.first.x, bin_index.first.y, bin_index.second.x, bin_index.second.y);

        for (int j = bin_index.first.x; j <= bin_index.second.x; ++j) {
            for (int k = bin_index.first.y; k <= bin_index.second.y; ++k) {
                auto bin_center = binCenter(j, k);
                Point2<double> d = Abs(bin_center - m);

                bins[j][k] += bell_shape(d.x, dim_m.x, dim_bin.x)
                            * bell_shape(d.y, dim_m.y, dim_bin.y);
            }
        }
    }
    for (int i = 0; i < num_cols; ++i) {
        for (int j = 0; j < num_rows; ++j) {
            double d_area = bins[i][j] - target_overlap_area;
            value_ += d_area * d_area;
        }
    }
    return value_;
}

const std::vector<Point2<double>>& Density::Backward(const std::vector<Point2<double>>& input) {
    fill(grad_.begin(), grad_.end(), Point2<double>(0, 0));
    double target_overlap_area = target_density * dim_bin.x * dim_bin.y;
    for (unsigned i = 0; i < placement_.numModules(); i++) {
        Module& module = placement_.module(i);
        Point2<double> m(module.centerX(), module.centerY());
        auto bin_index = getBeginEndBinIndex(module);
        for (int j = bin_index.first.x; j <= bin_index.second.x; ++j) {
            for (int k = bin_index.first.y; k <= bin_index.second.y; ++k) {
                auto bin_center = binCenter(j, k);
                Point2<double> d = Abs(bin_center - m);
                grad_[i].x += 2 * (bins[j][k] - target_overlap_area) 
                    * bell_shape(d.y, module.height(), dim_bin.y)
                    * grad_bell_shape(m.x, bin_center.x, module.width(), dim_bin.x);
                grad_[i].y += 2 * (bins[j][k] - target_overlap_area) 
                    * bell_shape(d.x, module.width(), dim_bin.x)
                    * grad_bell_shape(m.y, bin_center.y, module.height(), dim_bin.y);
            }
        }
    }
    return grad_;
}
