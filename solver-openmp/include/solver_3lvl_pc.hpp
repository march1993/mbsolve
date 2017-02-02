#ifndef OPENMP_3LVL_PC_H
#define OPENMP_3LVL_PC_H

#include <solver-openmp_EXPORTS.hpp>
#include <Solver.hpp>
#include <DensityMatrixData.hpp>

namespace mbsolve {

static const unsigned int MaxRegions = 8;

struct sim_constants
{
    real M_CE;
    real M_CH;
    real M_CP;
    real sigma;

    /* transition frequencies */
    real w12;
    real w13;
    real w23;

    /* dipoles */
    real d23;

    /* anticrossing frequencies */
    real O12;

    /* scattering rates (inverse lifetimes) */
    real tau11;
    real tau12;
    real tau13;
    real tau21;
    real tau22;
    real tau23;
    real tau31;
    real tau32;
    real tau33;

    /* dephasing rates */
    real gamma12;
    real gamma13;
    real gamma23;

    unsigned int idx_start;
    unsigned int idx_end;

    real d_x;
    real d_t;
};

class CopyListEntry
{
protected:
    Result *m_res;
    unsigned int m_count;
    unsigned int m_interval;
    unsigned int m_position;

public:
    CopyListEntry(Result *result, unsigned int count,
		  unsigned int position, unsigned int interval) :
	m_res(result), m_count(count),
	m_position(position), m_interval(interval)
    {
    }

    ~CopyListEntry()
    {
    }

    virtual real *getSrc() const = 0;

    real *getDst(unsigned int idx) const {
	return m_res->data(idx/m_interval);
    }

    unsigned int getSize() const { return sizeof(real) * m_count; }

    unsigned int getCount() const { return m_count; }

    bool record(unsigned int idx) const { return (idx % m_interval) == 0; }
};

class CLEField : public CopyListEntry
{
private:
    real *m_address;
public:
    CLEField(real *address, Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_address(address)
    {
    }

    real *getSrc() const
    {
	return m_address + m_position;
    }
};

class CLEDensity : public CopyListEntry
{
private:
    DensityMatrixData *m_dm;
    unsigned int m_row;
    unsigned int m_col;

public:
    CLEDensity(DensityMatrixData *dm, unsigned int row, unsigned int col,
	     Result *result, unsigned int count,
	     unsigned int position, unsigned int interval) :
	CopyListEntry(result, count, position, interval), m_dm(dm), m_row(row),
	m_col(col)
    {
    }

    real *getSrc() const
    {
	return m_dm->oldDM(m_row, m_col) + m_position;
    }
};

class SOLVER_OPENMP_EXPORT SolverOMP_3lvl_pc : public ISolver
{
public:
    SolverOMP_3lvl_pc(const Device& device, const Scenario& scenario);

    ~SolverOMP_3lvl_pc();

    std::string getName() const;

    void run() const;

private:
    inline void estimate_step(int i, real src) const;

    DensityMatrixData *m_dm;

    real *m_h;
    real *m_e;
    real *m_e_est;

    std::vector<CopyListEntry *> m_copyListBlack;
    std::vector<CopyListEntry *> m_copyListRed;
};

}

#endif
