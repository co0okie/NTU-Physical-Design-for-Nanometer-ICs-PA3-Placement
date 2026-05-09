#include "ObjectiveFunction.h"

#include "cstdio"
#include "cfloat"

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
        double &max_x = max_of_net[i].x, &min_x = min_of_net[i].x;
        double &max_y = max_of_net[i].y, &min_y = min_of_net[i].y;
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            max_x = std::max(max_x, pin.x());
            min_x = std::min(min_x, pin.x());
            max_y = std::max(max_y, pin.y());
            min_y = std::min(min_y, pin.y());
        }
        double cache;
        double 
            numerator_x_max = 0, &denominator_x_max = denominator_max_of_net[i].x = 0,
            numerator_x_min = 0, &denominator_x_min = denominator_min_of_net[i].x = 0,
            numerator_y_max = 0, &denominator_y_max = denominator_max_of_net[i].y = 0,
            numerator_y_min = 0, &denominator_y_min = denominator_min_of_net[i].y = 0;
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            numerator_x_max += pin.x() * (cache = std::exp((pin.x() - max_x) / gamma));
            denominator_x_max += cache;
            numerator_x_min += pin.x() * (cache = std::exp((-pin.x() - -min_x) / gamma));
            denominator_x_min += cache;
            numerator_y_max += pin.y() * (cache = std::exp((pin.y() - max_y) / gamma));
            denominator_y_max += cache;
            numerator_y_min += pin.y() * (cache = std::exp((-pin.y() - -min_y) / gamma));
            denominator_y_min += cache;
        }
        max_wa_of_net[i].x = numerator_x_max / denominator_x_max;
        max_wa_of_net[i].y = numerator_y_max / denominator_y_max;
        min_wa_of_net[i].x = numerator_x_min / denominator_x_min;
        min_wa_of_net[i].y = numerator_y_min / denominator_y_min;
        value_ += max_wa_of_net[i].x - min_wa_of_net[i].x +
                  max_wa_of_net[i].y - min_wa_of_net[i].y;
    }
    return value_;
}

const std::vector<Point2<double>>& Wirelength::Backward(
const std::vector<Point2<double>>& input) {
    grad_.assign(placement_.numModules(), Point2<double>(0, 0));
    for (unsigned i = 0; i < placement_.numNets(); i++) {
        Net& net = placement_.net(i);
        for (unsigned j = 0; j < net.numPins(); j++) {
            Pin& pin = net.pin(j);
            grad_[pin.moduleId()].x += 
                std::exp((pin.x() - max_of_net[i].x) / gamma) / 
                denominator_max_of_net[i].x * 
                (1 + (pin.x() - max_wa_of_net[i].x) / gamma)
                -
                std::exp((-pin.x() - -min_of_net[i].x) / gamma) / 
                denominator_min_of_net[i].x * 
                (1 - (pin.x() - min_wa_of_net[i].x) / gamma);
            grad_[pin.moduleId()].y +=
                std::exp((pin.y() - max_of_net[i].y) / gamma) / 
                denominator_max_of_net[i].y * 
                (1 + (pin.y() - max_wa_of_net[i].y) / gamma)
                -
                std::exp((-pin.y() - -min_of_net[i].y) / gamma) / 
                denominator_min_of_net[i].y * 
                (1 - (pin.y() - min_wa_of_net[i].y) / gamma);
        }
    }
    return grad_;
}
