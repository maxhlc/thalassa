#include <pythalassa/pythalassa.hpp>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

#include <Eigen/Core>

#include <cthalassa/cthalassa.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pythalassa, m) {
    // Module documentation
    m.doc() = "Python interface for the THALASSA orbit propagator";

    // Physical model enumerations
    py::enum_<cthalassa::ModelGravity>(m, "ModelGravity")
        .value("SPHERICAL", cthalassa::SPHERICAL)
        .value("NONSPHERICAL", cthalassa::NONSPHERICAL)
        .export_values();
    py::enum_<cthalassa::ModelSun>(m, "ModelSun")
        .value("SUN_DISABLED", cthalassa::SUN_DISABLED)
        .value("SUN_ENABLED", cthalassa::SUN_ENABLED)
        .export_values();
    py::enum_<cthalassa::ModelMoon>(m, "ModelMoon")
        .value("MOON_DISABLED", cthalassa::MOON_DISABLED)
        .value("MOON_ENABLED", cthalassa::MOON_ENABLED)
        .export_values();
    py::enum_<cthalassa::ModelDrag>(m, "ModelDrag")
        .value("DRAG_DISABLED", cthalassa::DRAG_DISABLED)
        .value("DRAG_WERTZ", cthalassa::DRAG_WERTZ)
        .value("DRAG_US76", cthalassa::DRAG_US76)
        .value("DRAG_J77", cthalassa::DRAG_J77)
        .value("DRAG_NRLMSISE00", cthalassa::DRAG_NRLMSISE00)
        .export_values();
    py::enum_<cthalassa::ModelFlux>(m, "ModelFlux")
        .value("FLUX_CONSTANT", cthalassa::FLUX_CONSTANT)
        .value("FLUX_VARIABLE", cthalassa::FLUX_VARIABLE)
        .export_values();
    py::enum_<cthalassa::ModelSRP>(m, "ModelSRP")
        .value("SRP_DISABLED", cthalassa::SRP_DISABLED)
        .value("SRP_ENABLED", cthalassa::SRP_ENABLED)
        .value("SRP_ENABLED_CONICAL", cthalassa::SRP_ENABLED_CONICAL)
        .export_values();
    py::enum_<cthalassa::ModelEphemerides>(m, "ModelEphemerides")
        .value("EPHEM_DE431", cthalassa::EPHEM_DE431)
        .value("EPHEM_SIMPLE", cthalassa::EPHEM_SIMPLE)
        .export_values();

    // Model structure
    py::class_<cthalassa::Model>(m, "Model")
        .def(py::init<>())
        .def_readwrite("insgrav", &cthalassa::Model::insgrav)
        .def_readwrite("isun", &cthalassa::Model::isun)
        .def_readwrite("imoon", &cthalassa::Model::imoon)
        .def_readwrite("idrag", &cthalassa::Model::idrag)
        .def_readwrite("iF107", &cthalassa::Model::iF107)
        .def_readwrite("iSRP", &cthalassa::Model::iSRP)
        .def_readwrite("iephem", &cthalassa::Model::iephem)
        .def_readwrite("gdeg", &cthalassa::Model::gdeg)
        .def_readwrite("gord", &cthalassa::Model::gord);

    // Paths structure
    py::class_<cthalassa::Paths>(m, "Paths")
        .def(py::init<>())
        .def_readwrite("phys_path", &cthalassa::Paths::phys_path)
        .def_readwrite("earth_path", &cthalassa::Paths::earth_path)
        .def_readwrite("kernel_path", &cthalassa::Paths::kernel_path)
        .def_readwrite("eop_path", &cthalassa::Paths::eop_path);

    // Equation enumeration
    py::enum_<cthalassa::Equations>(m, "Equations")
        .value("COWELL", cthalassa::COWELL)
        .value("EDROMO_T", cthalassa::EDROMO_T)
        .value("EDROMO_C", cthalassa::EDROMO_C)
        .value("EDROMO_L", cthalassa::EDROMO_L)
        .value("KS_T", cthalassa::KS_T)
        .value("KS_L", cthalassa::KS_L)
        .value("STISCHE_T", cthalassa::STISCHE_T)
        .value("STISCHE_L", cthalassa::STISCHE_L)
        .export_values();

    // Settings structure
    py::class_<cthalassa::Settings>(m, "Settings")
        .def(py::init<>())
        .def_readwrite("tol", &cthalassa::Settings::tol)
        .def_readwrite("mxstep", &cthalassa::Settings::mxstep)
        .def_readwrite("imcoll", &cthalassa::Settings::imcoll)
        .def_readwrite("eqs", &cthalassa::Settings::eqs)
        .def_readwrite("tspan", &cthalassa::Settings::tspan)
        .def_readwrite("tstep", &cthalassa::Settings::tstep);

    // Spacecraft structure
    py::class_<cthalassa::Spacecraft>(m, "Spacecraft")
        .def(py::init<>())
        .def_readwrite("mass", &cthalassa::Spacecraft::mass)
        .def_readwrite("area_drag", &cthalassa::Spacecraft::area_drag)
        .def_readwrite("area_srp", &cthalassa::Spacecraft::area_srp)
        .def_readwrite("cd", &cthalassa::Spacecraft::cd)
        .def_readwrite("cr", &cthalassa::Spacecraft::cr);

    // Propagator class
    py::class_<pythalassa::Propagator>(m, "Propagator")
        .def(py::init<const cthalassa::Model &, const cthalassa::Paths &, const cthalassa::Settings &, const cthalassa::Spacecraft &>())
        .def("propagate", py::overload_cast<const Eigen::VectorXd &, const Eigen::VectorXd &>(&pythalassa::Propagator::propagate, py::const_))
        .def_property_readonly("model", &pythalassa::Propagator::getModel)
        .def_property_readonly("paths", &pythalassa::Propagator::getPaths)
        .def_property("settings", &pythalassa::Propagator::getSettings, &pythalassa::Propagator::setSettings)
        .def_property("spacecraft", &pythalassa::Propagator::getSpacecraft, &pythalassa::Propagator::setSpacecraft);
}