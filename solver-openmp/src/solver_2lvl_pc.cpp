#include <solver_2lvl_pc.hpp>

#include <boost/foreach.hpp>

namespace mbsolve{

static SolverFactory<SolverOMP_2lvl_pc> factory("openmp-2lvl-pc");

static struct sim_constants gsc[MaxRegions];

inline unsigned int get_region(unsigned int idx)
{
    for (unsigned int i = 0; i < MaxRegions; i++) {
	if (idx <= gsc[i].idx_end) {
	    return i;
	}
    }
    return 0;
}

inline void
SolverOMP_2lvl_pc::estimate_step(int i, real src) const
{
    int region = get_region(i);

    real rho11_o  = m_dm->oldDM(0, 0)[i];
    real rho11_e  = rho11_o;
    real rho22_o  = m_dm->oldDM(1, 1)[i];
    real rho22_e  = rho22_o;
    real rho12r_o = m_dm->oldDM(1, 0)[i];
    real rho12r_e = rho12r_o;
    real rho12i_o = m_dm->oldDM(0, 1)[i];
    real rho12i_e = rho12i_o;
    real field_o = m_e[i];
    real field_e = field_o;

    for (int pc_step = 0; pc_step < 4; pc_step++) {
	/* execute prediction - correction steps */

	real rho11  = 0.5 * (rho11_o + rho11_e);
	real rho22  = 0.5 * (rho22_o + rho22_e);
	real rho12r = 0.5 * (rho12r_o + rho12r_e);
	real rho12i = 0.5 * (rho12i_o + rho12i_e);
	real OmRabi = 0.5 * gsc[region].d12 * (field_o + field_e);

	rho11_e = rho11_o + gsc[region].d_t *
	    (- 2.0 * OmRabi * rho12i - gsc[region].tau1 * rho11);

	rho12i_e = rho12i_o + gsc[region].d_t *
	    (- gsc[region].w12 * rho12r
	     + OmRabi * (rho11 - rho22)
	     - gsc[region].gamma12 * rho12i);

	rho12r_e = rho12r_o + gsc[region].d_t *
	    (+ gsc[region].w12 * rho12i
	     - gsc[region].gamma12 * rho12r);

	rho22_e = rho22_o + gsc[region].d_t *
	    (+ 2.0 * OmRabi * rho12i
	     + gsc[region].tau1 * rho11);

	real j = 0;
	real p_t = gsc[region].M_CP * gsc[region].d12 *
	    (gsc[region].w12 * rho12i -
	     gsc[region].gamma12 * rho12r);

	field_e = field_o + gsc[region].M_CE *
	    (-j - p_t + (m_h[i + 1] - m_h[i])/gsc[region].d_x);

	if (i == 0) {
	    /* TODO rework soft source for pred-corr scheme */
	    //m_e_est[i] += src; /* soft source */
	    field_e = src; /* hard source */
	}
    }

    /* final update step */
    m_dm->newDM(0, 0)[i] = rho11_e;
    m_dm->newDM(0, 1)[i] = rho12i_e;
    m_dm->newDM(1, 0)[i] = rho12r_e;
    m_dm->newDM(1, 1)[i] = rho22_e;

    m_e[i] = field_e;
}

SolverOMP_2lvl_pc::SolverOMP_2lvl_pc(const Device& device,
				     const Scenario& scenario) :
    ISolver(device, scenario)
{
    /* total device length */
    Quantity length = device.XDim();

    /* minimum relative permittivity */
    Quantity minRelPermittivity = device.MinRelPermittivity();

    /* TODO: sanity check scenario? */
    if (m_scenario.NumGridPoints % 32 != 0) {
	throw std::invalid_argument("Number of grid points must be multiple"
				    " of 32");
    }

    /* determine grid point and time step size */
    real C = 0.5; /* courant number */
    real velocity = sqrt(MU0() * EPS0() * minRelPermittivity());
    m_scenario.GridPointSize = length()/(m_scenario.NumGridPoints - 1);
    real timestep  = C * m_scenario.GridPointSize * velocity;
    m_scenario.NumTimeSteps = ceil(m_scenario.SimEndTime/timestep) + 1;
    m_scenario.TimeStepSize = m_scenario.SimEndTime /
	(m_scenario.NumTimeSteps - 1);

    /* determine border indices and initialize region settings */
    if (device.Regions.size() > MaxRegions) {
	throw std::invalid_argument("Too many regions requested");
    }
    struct sim_constants sc[MaxRegions];

    unsigned int i = 0;
    BOOST_FOREACH(Region reg, device.Regions) {
	if (i > 0) {
	    gsc[i - 1].idx_end = round(reg.X0()/m_scenario.GridPointSize) - 1;
	}
	gsc[i].idx_start = round(reg.X0()/m_scenario.GridPointSize);
	gsc[i].M_CE = m_scenario.TimeStepSize/(EPS0() * reg.RelPermittivity());
	gsc[i].M_CH = m_scenario.TimeStepSize/(MU0() *
					      m_scenario.GridPointSize);
	/* TODO: overlap factor? */
	gsc[i].M_CP = -2.0 * reg.DopingDensity * HBAR();
	gsc[i].sigma = 2.0 * sqrt(EPS0 * reg.RelPermittivity/MU0) * reg.Losses;

	gsc[i].w12 = (reg.TransitionFrequencies.size() < 1) ? 0.0 :
	    reg.TransitionFrequencies[0]();
	/* TODO rename to rabi freqs or something */
	gsc[i].d12 = (reg.DipoleMoments.size() < 1) ? 0.0 :
	    reg.DipoleMoments[0]() * E0()/HBAR();
	gsc[i].tau1 = (reg.ScatteringRates.size() < 1) ? 0.0 :
	    reg.ScatteringRates[0]();
	gsc[i].gamma12 = (reg.DephasingRates.size() < 1) ? 0.0 :
	    reg.DephasingRates[0]();

	gsc[i].d_x = m_scenario.GridPointSize;
	gsc[i].d_t = m_scenario.TimeStepSize;

	if (reg.DopingDensity() < 1.0) {
	    gsc[i].dm11_init = 0.0;
	    gsc[i].dm22_init = 0.0;
	} else {
	    gsc[i].dm11_init = 0.0;
	    gsc[i].dm22_init = 1.0;
	}

	i++;
    }
    if (i > 0) {
	gsc[i - 1].idx_end = m_scenario.NumGridPoints - 1;
    }

    m_dm = new DensityMatrixData(m_scenario.NumGridPoints, 2, 5);

    m_h = new real[m_scenario.NumGridPoints + 1];
    m_e = new real[m_scenario.NumGridPoints];
    //m_e_est = new real[m_scenario.NumGridPoints];
    /*
#pragma omp parallel for
    for (int i = 0; i < m_scenario.NumGridPoints; i++) {
        int region = get_region(i);

        m_dm->oldDM(0, 0)[i] = gsc[region].dm11_init;
        m_dm->newDM(0, 0)[i] = gsc[region].dm11_init;
        m_dm->rhs(0, 0, 0)[i] = gsc[region].dm11_init;

        m_dm->oldDM(1, 1)[i] = gsc[region].dm22_init;
        m_dm->newDM(1, 1)[i] = gsc[region].dm22_init;
        m_dm->rhs(1, 1, 0)[i] = gsc[region].dm22_init;

        m_dm->oldDM(0, 1)[i] = 0.0;
        m_dm->newDM(0, 1)[i] = 0.0;
        m_dm->rhs(0, 1, 0)[i] = 0.0;

        m_dm->oldDM(1, 0)[i] = 0.0;
        m_dm->newDM(1, 0)[i] = 0.0;
        m_dm->rhs(1, 0, 0)[i] = 0.0;

        m_e[i] = 0.0;
        m_e_est[i] = 0.0;
    }

#pragma omp parallel for
    for (int i = 0; i < m_scenario.NumGridPoints + 1; i++) {
        m_h[i] = 0.0;
    }
    */
    /* set up results transfer data structures */
    BOOST_FOREACH(Record rec, m_scenario.Records) {
	unsigned int row_ct = m_scenario.SimEndTime/rec.Interval;
	unsigned int interval = ceil(1.0 * m_scenario.NumTimeSteps/row_ct);
	unsigned int position_idx;
	unsigned int col_ct;

	if (rec.Position() < 0.0) {
	    /* copy complete grid */
	    position_idx = 0;
	    col_ct = m_scenario.NumGridPoints;
	} else {
	    position_idx = round(rec.Position()/m_scenario.GridPointSize);
	    col_ct = 1;
	}

	/* allocate result memory */
	Result *res = new Result(rec.Name, col_ct, row_ct);
	m_results.push_back(res);

	/* create copy list entry */
	CopyListEntry *entry;
	if (rec.Type == EField) {
	    entry = new CLEField(m_e, res, col_ct, position_idx,
				 interval);
	    m_copyListRed.push_back(entry);
	} else if (rec.Type == HField) {
	    /* TODO: numGridPoints + 1 */
	    entry = new CLEField(m_h, res, col_ct, position_idx,
				 interval);
	    m_copyListBlack.push_back(entry);
	} else if (rec.Type == Density) {
	    if ((rec.I - 1 < 2) && (rec.J - 1 < 2)) {
		if (rec.I == rec.J) {
		    /* main diagonal entry */
		    entry = new CLEDensity(m_dm, rec.I - 1, rec.J - 1, res,
					   col_ct, position_idx,
					   interval);
		    //m_copyListBlack.push_back(entry);
		    m_copyListRed.push_back(entry);
		} else {
		    /* off-diagonal entry */
		    /* TODO */
		    /* if complex */
		    /* create two list entries */
		    /* create two Results, or one complex Result */

		    /* real part: GetSrcDensity(&dm, rec.I, rec.J); */
		    /* imag part: GetSrcDensity(&dm, rec.J, rec.I); */
		}
	    } else {
	    // throw exc
	    }
	} else {
	    // throw exc
	}
    }
}

SolverOMP_2lvl_pc::~SolverOMP_2lvl_pc()
{
    /* delete copy lists */
    BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
	delete entry;
    }
    BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
	delete entry;
    }

    delete m_dm;
    delete[] m_h;
    delete[] m_e;
    //delete[] m_e_est;
}

std::string
SolverOMP_2lvl_pc::getName() const
{
    return factory.getName();
}

void
SolverOMP_2lvl_pc::run() const
{
#pragma omp parallel
    {
	/* parallel initialization */
#pragma omp for schedule(static)
	for (int i = 0; i < m_scenario.NumGridPoints; i++) {
	    int region = get_region(i);

	    m_dm->oldDM(0, 0)[i] = gsc[region].dm11_init;
	    m_dm->newDM(0, 0)[i] = gsc[region].dm11_init;
	    m_dm->rhs(0, 0, 0)[i] = gsc[region].dm11_init;

	    m_dm->oldDM(1, 1)[i] = gsc[region].dm22_init;
	    m_dm->newDM(1, 1)[i] = gsc[region].dm22_init;
	    m_dm->rhs(1, 1, 0)[i] = gsc[region].dm22_init;

	    m_dm->oldDM(0, 1)[i] = 0.0;
	    m_dm->newDM(0, 1)[i] = 0.0;
	    m_dm->rhs(0, 1, 0)[i] = 0.0;

	    m_dm->oldDM(1, 0)[i] = 0.0;
	    m_dm->newDM(1, 0)[i] = 0.0;
	    m_dm->rhs(1, 0, 0)[i] = 0.0;

	    m_e[i] = 0.0;
	    //m_e_est[i] = 0.0;
	    m_h[i] = 0.0;
	    if (i == m_scenario.NumGridPoints - 1) {
		m_h[i + 1] = 0.0;
	    }
	}

	/* main loop */
	for (unsigned int n = 0; n < m_scenario.NumTimeSteps; n++) {
	    /* calculate source value */
	    /* TODO optimize divisions?`*/
	    real f_0 = 2e14;
	    real t = n * m_scenario.TimeStepSize;
	    real T_p = 20/f_0;
	    real gamma = 2 * t/T_p - 1;
	    real E_0 = 4.2186e9;
	    real src = E_0 * 1/std::cosh(10 * gamma) * sin(2 * M_PI * f_0 * t);
	    src = src/2; // pi pulse
	    /* TODO:  private(src) ? */
	    /* TODO: masteronly? */

	    /* update dm and e in parallel */
#pragma omp for schedule(static)
	    for (int i = 0; i < m_scenario.NumGridPoints; i++) {
		estimate_step(i, src);
	    }

	    /* update h in parallel */
#pragma omp for schedule(static)
	    for (int i = 0; i < m_scenario.NumGridPoints; i++) {
		int region = get_region(i);

		if (i != 0) {
		    m_h[i] += gsc[region].M_CH * (m_e[i] - m_e[i - 1]);
		}
	    }

#pragma omp master
	    {
		/* toggle density matrix data */
		m_dm->next();

		/* TODO parallel to update (or rather parallel copy) */
		/* gather e field and dm entries in copy stream */
		BOOST_FOREACH(CopyListEntry *entry, m_copyListRed) {
		    if (entry->record(n)) {
			std::copy(entry->getSrc(),
				  entry->getSrc() + entry->getCount(),
				  entry->getDst(n));
		    }
		}

		/* TODO parallel to update */
		/* gather h field in copy stream */
		BOOST_FOREACH(CopyListEntry *entry, m_copyListBlack) {
		    if (entry->record(n)) {
			std::copy(entry->getSrc(),
				  entry->getSrc() + entry->getCount(),
				  entry->getDst(n));
		    }
		}
	    }
#pragma omp barrier
	}
    }
}

}