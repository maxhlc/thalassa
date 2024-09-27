#ifndef CTHALASSA_HPP_
#define CTHALASSA_HPP_

#include <cstring>
#include <mutex>
#include <shared_mutex>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef CTHALASSA_USE_EIGEN
#include <Eigen/Core>
#endif

namespace cthalassa::internal {
    extern "C" {
    // clang-format off
        #include <cthalassa/cthalassa.h>
    // clang-format on
    }
} // namespace cthalassa::internal

namespace cthalassa {

    /// @brief Gravity models
    enum ModelGravity { SPHERICAL = 0, NONSPHERICAL = 1 };

    /// @brief Sun model
    enum ModelSun { SUN_DISABLED = 0, SUN_ENABLED = 1 };

    /// @brief Moon model
    enum ModelMoon { MOON_DISABLED = 0, MOON_ENABLED = 1 };

    /// @brief Drag model
    enum ModelDrag { DRAG_DISABLED = 0, DRAG_WERTZ = 1, DRAG_US76 = 2, DRAG_J77 = 3, DRAG_NRLMSISE00 = 4 };

    /// @brief Flux model
    enum ModelFlux { FLUX_CONSTANT = 0, FLUX_VARIABLE = 1 };

    /// @brief Solar Radiation Pressure (SRP) model
    enum ModelSRP { SRP_DISABLED = 0, SRP_ENABLED = 1, SRP_ENABLED_CONICAL = 2 };

    /// @brief Ephemerides
    enum ModelEphemerides { EPHEM_DE431 = 1, EPHEM_SIMPLE = 2 };

    /**
     * @brief Structure containing propagator model settings
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct Model {
        /// @brief Gravity model
        ModelGravity insgrav = NONSPHERICAL;

        /// @brief Sun model
        ModelSun isun = SUN_ENABLED;

        /// @brief Moon model
        ModelMoon imoon = MOON_ENABLED;

        /// @brief Drag model
        ModelDrag idrag = DRAG_WERTZ;

        /// @brief Flux model
        ModelFlux iF107 = FLUX_VARIABLE;

        /// @brief Solar Radiation Pressure (SRP) model
        ModelSRP iSRP = SRP_ENABLED_CONICAL;

        /// @brief Ephemerides
        ModelEphemerides iephem = EPHEM_SIMPLE;

        /// @brief Gravity model degree
        int gdeg = 2;

        /// @brief Gravity model order
        int gord = 2;

        /**
         * @brief Equality comparison with another Model structure
         *
         * @author Max Hallgarten La Casta
         *
         * @param rhs Other Model object
         * @return true Model match
         * @return false Model mismatch
         */
        bool operator==(const Model &rhs) const {
            return (insgrav == rhs.insgrav) && (isun == rhs.isun) && (imoon == rhs.imoon) && (idrag == rhs.idrag) && (iF107 == rhs.iF107) &&
                   (iSRP == rhs.iSRP) && (iephem == rhs.iephem) && (gdeg == rhs.gdeg) && (gord == rhs.gord);
        }

        /**
         * @brief Inequality comparison with another Model structure
         *
         * @author Max Hallgarten La Casta
         *
         * @param rhs Other Model structure
         * @return true Model mismatch
         * @return false Model match
         */
        bool operator!=(const Model &rhs) const { return !(*this == rhs); }

        /// @brief Method to implicitly cast the structure to its internal equivalent in CTHALASSA
        operator cthalassa::internal::THALASSAPhysicalModelStruct() const {
            return cthalassa::internal::THALASSAPhysicalModelStruct{insgrav, isun, imoon, idrag, iF107, iSRP, iephem, gdeg, gord};
        }
    } Model;

    /**
     * @brief Structure containing filepaths for THALASSA
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct Paths {
        /// @brief Filepath for THALASSA physical constants
        std::string phys_path;

        /// @brief Filepath for Earth geopotential model
        std::string earth_path;

        /// @brief Filepath for SPICE kernel
        std::string kernel_path;

        /// @brief Filepath for EOP
        std::string eop_path;

        /**
         * @brief Equality comparison with another Paths structure
         *
         * @author Max Hallgarten La Casta
         *
         * @param rhs Other Paths structure
         * @return true Paths match
         * @return false Paths mismatch
         */
        bool operator==(const Paths &rhs) const { return (phys_path == rhs.phys_path) && (earth_path == rhs.earth_path) && (kernel_path == rhs.kernel_path) && (eop_path == rhs.eop_path); }

        /**
         * @brief Inequality comparison with another Paths structure
         *
         * @author Max Hallgarten La Casta
         *
         * @param rhs Other Paths structure
         * @return true Paths mismatch
         * @return false Paths match
         */
        bool operator!=(const Paths &rhs) const { return !(*this == rhs); }

        /// @brief Method to implicitly cast the structure to its internal equivalent in CTHALASSA
        operator cthalassa::internal::THALASSAPathStruct() const {
            // Declare new path struct
            cthalassa::internal::THALASSAPathStruct paths;

            // Set physical constants path
            strcpy(paths.phys_path, phys_path.c_str());
            paths.phys_path_len = phys_path.size();

            // Set Earth model path
            strcpy(paths.earth_path, earth_path.c_str());
            paths.earth_path_len = earth_path.size();

            // Set SPICE kernel path
            strcpy(paths.kernel_path, kernel_path.c_str());
            paths.kernel_path_len = kernel_path.size();

            // Set EOP path
            strcpy(paths.eop_path, eop_path.c_str());
            paths.eop_path_len = eop_path.size();

            // Return path struct
            return paths;
        }
    } Paths;

    /// @brief Propagator formulations
    enum Equations { COWELL = 1, EDROMO_T = 2, EDROMO_C = 3, EDROMO_L = 4, KS_T = 5, KS_L = 6, STISCHE_T = 7, STISCHE_L = 8 };

    /**
     * @brief Structure containing propagator settings
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct Settings {
        /// @brief Integration tolerance
        double tol = 1.0E-13;

        /// @brief Maximum number of integration/output steps
        int mxstep = 1000000;

        /// @brief Moon collision check flag
        int imcoll = 0;

        /// @brief Propagator formulation
        Equations eqs = EDROMO_T;

        /// @brief Propagation time span [solar days]
        double tspan;

        /// @brief Step size [solar days]
        double tstep;

        /// @brief Method to implicitly cast the structure to its internal equivalent in CTHALASSA
        operator cthalassa::internal::THALASSAPropagatorStruct() const {
            return cthalassa::internal::THALASSAPropagatorStruct{tol, tspan, tstep, mxstep, imcoll, eqs};
        }
    } Settings;

    /**
     * @brief Structure containing spacecraft parameters
     *
     * @author Max Hallgarten La Casta
     */
    typedef struct Spacecraft {
        /// @brief Spacecraft mass [kg]
        double mass;

        /// @brief Drag area [m^2]
        double area_drag;

        /// @brief Solar Radiation Pressure (SRP) area [m^2]
        double area_srp;

        /// @brief Coefficient of drag [-]
        double cd;

        /// @brief Coefficient of reflectivity [-]
        double cr;

        /// @brief Method to implicitly cast the structure to its internal equivalent in CTHALASSA
        operator cthalassa::internal::THALASSAObjectStruct() const { return cthalassa::internal::THALASSAObjectStruct{mass, area_drag, area_srp, cd, cr}; }
    } Spacecraft;

    /**
     * @brief THALASSA propagator instances handler
     *
     * @author Max Hallgarten La Casta
     */
    class PropagatorInstances {

    protected:
        /// @brief Number of instantiated THALASSA propagator objects
        static size_t instances_;

        /// @brief Lock for making instance changes
        static std::mutex instancesMutex_;

        /// @brief Lock for propagations
        static std::mutex propagationMutex_;

        /// @brief Shared lock for propagations
        static std::shared_mutex propagationSharedMutex_;

        /// @brief Shared local lock for propagations
        mutable std::shared_mutex propagationLocalSharedMutex_;

        /// @brief Model settings
        static Model model_;

        /// @brief Filepaths
        static Paths paths_;

    public:
        /**
         * @brief Construct a new Propagator Instances object
         *
         * @author Max Hallgarten La Casta
         */
        PropagatorInstances() {
            // Increment the number of instances
            ++instances_;
        }

        /**
         * @brief Construct a new Propagator Instances object
         *
         * @author Max Hallgarten La Casta
         */
        PropagatorInstances(const PropagatorInstances &) {
            // Increment the number of instances
            ++instances_;
        }

        /**
         * @brief Destroy the Propagator Instances object
         *
         * @author Max Hallgarten La Casta
         */
        ~PropagatorInstances() {
            // Decrement the number of instances
            --instances_;
        }
    };

    /**
     * @brief THALASSA orbit propagator
     *
     * @author Max Hallgarten La Casta
     */
    class Propagator : public PropagatorInstances {

    private:
        /// @brief Propagator settings
        Settings settings_;

        /// @brief Spacecraft
        Spacecraft spacecraft_;

    public:
        /**
         * @brief Construct a new Propagator object
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] model Model settings
         * @param[in] paths Filepaths
         * @param[in] settings Propagator settings
         * @param[in] spacecraft Spacecraft
         */
        Propagator(const Model &model, const Paths &paths, const Settings &settings, const Spacecraft &spacecraft);

        /**
         * @brief Destroy the Propagator object
         *
         * @author Max Hallgarten La Casta
         */
        ~Propagator();

        /**
         * @brief Propagate a state
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] tStart Initial time
         * @param[in] tEnd Final time
         * @param[in] tStep Output time step
         * @param[in] stateIn Initial state
         * @param[out] timesOut Output times
         * @param[out] statesOut Output states
         */
        void propagate(const double &tStart, const double &tEnd, const double &tStep, const std::vector<double> &stateIn, std::vector<double> &timesOut,
                       std::vector<std::vector<double>> &statesOut) const;

        /**
         * @brief Propagate a state
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] times Times
         * @param[in] stateIn Initial state
         * @param[out] statesOut Output states
         */
        void propagate(const std::vector<double> &times, const std::vector<double> &stateIn, std::vector<std::vector<double>> &statesOut) const;

#ifdef CTHALASSA_USE_EIGEN

        /**
         * @brief Propagate a state
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] tStart Initial time
         * @param[in] tEnd Final time
         * @param[in] tStep Output time step
         * @param[in] stateIn Initial state
         * @param[out] timesOut Output times
         * @param[out] statesOut Output states
         */
        void propagate(const double &tStart, const double &tEnd, const double &tStep, const Eigen::VectorXd &stateIn, Eigen::VectorXd &timesOut,
                       Eigen::MatrixXd &statesOut) const;

        /**
         * @brief Propagate a state
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] times Times
         * @param[in] stateIn Initial state
         * @param[out] statesOut Output states
         */
        void propagate(const Eigen::VectorXd &times, const Eigen::VectorXd &stateIn, Eigen::MatrixXd &statesOut) const;

#endif

        /**
         * @brief Get the model settings
         *
         * @author Max Hallgarten La Casta
         *
         * @return Model Model settings
         */
        Model getModel() const;

        /**
         * @brief Get the filepaths
         *
         * @author Max Hallgarten La Casta
         *
         * @return Paths Filepaths
         */
        Paths getPaths() const;

        /**
         * @brief Set the propagator settings
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] settings Propagator settings
         */
        void setSettings(const Settings &settings);

        /**
         * @brief Get the propagator settings
         *
         * @author Max Hallgarten La Casta
         *
         * @return Settings Propagator settings
         */
        Settings getSettings();

        /**
         * @brief Set the spacecraft
         *
         * @author Max Hallgarten La Casta
         *
         * @param[in] spacecraft Spacecraft
         */
        void setSpacecraft(const Spacecraft &spacecraft);

        /**
         * @brief Get the spacecraft
         *
         * @author Max Hallgarten La Casta
         *
         * @return Spacecraft Spacecraft
         */
        Spacecraft getSpacecraft();
    };

} // namespace cthalassa

#endif