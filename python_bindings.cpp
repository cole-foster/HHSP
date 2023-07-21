#include <pybind11/pybind11.h>

namespace py = pybind11;

#include "HierarchicalHSP.hpp"


class HHSP_Bindings {
public:
    float multiplier_ = 0;

    HierarchicalHSP* alg_hhsp{nullptr};

    HHSP_Bindings(float multiplier): multiplier_(multiplier) {
        alg_hhsp = new HierarchicalHSP();
    };
    ~HHSP_Bindings() {
        if (alg_hhsp) delete alg_hhsp;
    };
    float multiply(float input) {
        return multiplier_ * input;
    }
    void setMultiplier(float multiplier) { multiplier_ = multiplier;}
    float getMultiplier() { return multiplier_;}
};


PYBIND11_MODULE(HHSP, handle) {
    handle.doc() = "Python Module for Hierarchical HSP";

    py::class_<HHSP_Bindings> (handle, "HHSP")
        .def(py::init<float>())
        .def_property("multiplier", &HHSP_Bindings::getMultiplier, &HHSP_Bindings::setMultiplier)
        .def("multiply", &HHSP_Bindings::multiply);
}