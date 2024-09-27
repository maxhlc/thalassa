#include <cthalassa/cthalassa.hpp>

int main() {
    // Declare THALASSA parameters
    cthalassa::Model model;
    cthalassa::Paths paths = {"./data/physical_constants.txt", "./data/earth_potential/GRIM5-S1.txt", "./data/kernels_to_load.furnsh", "./data/eop_data.txt"};
    cthalassa::Settings settings;
    cthalassa::Spacecraft spacecraft = {+8500.000000000000E+00, +13.00000000000000E+00, +13.00000000000000E+00, +2.200000000000000E+00, +1.500000000000000E+00};

    // Create THALASSA propagator
    cthalassa::Propagator propagator(model, paths, settings, spacecraft);

    // Declare initial state
    std::vector<double> stateIn = {7.000000000000000E+03, 0.000000000000000E+00, 0.000000000000000E+00,
                                   0.000000000000000E+00, 7.546053286792848E+00, 0.000000000000000E+00};

    // Set times
    const double tStart = 59043;
    const double tEnd = tStart + 7;
    const double tStep = 1;

    // Declare output vectors
    std::vector<double> timesOut;
    std::vector<std::vector<double>> statesOut;

    // Propagate with THALASSA
    propagator.propagate(tStart, tEnd, tStep, stateIn, timesOut, statesOut);

    // Return zero
    return 0;
}