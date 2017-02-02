#ifndef SOLVERGENERIC_H
#define SOLVERGENERIC_H

#include <solver-generic_EXPORTS.hpp>
#include <Solver.hpp>

namespace mbsolve {

class SOLVER_GENERIC_EXPORT SolverGeneric : public ISolver
{
public:
    SolverGeneric(const Device& device, const Scenario& scenario);

    ~SolverGeneric();

    std::string getName() const;

    void run() const;
};

}

#endif
