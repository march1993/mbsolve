/*
 * mbsolve: Framework for solving the Maxwell-Bloch/-Lioville equations
 *
 * Copyright (c) 2016, Computational Photonics Group, Technical University of
 * Munich.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_MALLOC

#include <mpi.h>
#include <numeric>
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <unsupported/Eigen/MatrixFunctions>
#include <common_mpich.hpp>
#include <solver_mpich_clvl_os_red.hpp>

namespace mbsolve {

static solver_factory<solver_mpich_clvl_os_red<2> > f2("mpich-2lvl-os-red");
static solver_factory<solver_mpich_clvl_os_red<3> > f3("mpich-3lvl-os-red");
static solver_factory<solver_mpich_clvl_os_red<6> > f6("mpich-6lvl-os-red");

/* redundant calculation overlap */
const unsigned int OL = 32;
const unsigned int VEC = 4;

template<unsigned int num_lvl, unsigned int num_adj>
void
fill_rodr_coeff(const Eigen::Matrix<complex, num_adj, num_adj>& eigenvec,
                const Eigen::Matrix<complex, num_adj, 1>& eigenval,
                sim_constants_clvl_os<num_lvl>& sc)
{
    /* creating sorting order (descending eigenvalues) */
    std::vector<size_t> perm_idx(num_adj);
    std::iota(perm_idx.begin(), perm_idx.end(), 0);
    std::sort(perm_idx.begin(), perm_idx.end(),
              [&eigenval](size_t i1, size_t i2) {
                auto imag1 = eigenval(i1).imag(), imag2 = eigenval(i2).imag();
                imag1 = imag1 > 0 ? imag1 : -imag1;
                imag2 = imag2 > 0 ? imag2 : -imag2;
                  return imag1 > imag2;
              });

    /* sort eigenvectors */
    Eigen::Matrix<real, num_adj, num_adj> Q =
        Eigen::Matrix<real, num_adj, num_adj>::Zero();
    for (int i = 0; i < num_adj/2; i++) {
        unsigned int i1 = perm_idx[2 * i];
        unsigned int i2 = perm_idx[2 * i + 1];

        Q.col(2 * i) = 1.0/sqrt(2) *
            (eigenvec.col(i1) + eigenvec.col(i2)).real();
        Q.col(2 * i + 1) = 1.0/sqrt(2) *
            (-eigenvec.col(i1) + eigenvec.col(i2)).imag();
    }
    if (num_adj % 2 != 0) {
        Q(num_adj - 1, num_adj - 1) = 1.0;
    }

    /* TODO optimize
     * ignore eigenvalues = 0
     * group eigenvalues with multiplicity >= 2
     */
    for (int i = 0; i < num_adj/2; i++) {
        unsigned int i1 = perm_idx[2 * i];
        unsigned int i2 = perm_idx[2 * i + 1];
        Eigen::Matrix<real, num_adj, num_adj> b =
            Eigen::Matrix<real, num_adj, num_adj>::Zero();

        /* give warning if eigenvalues do not match */
        std::cout << "Eigenvalues pairwise: " <<
            eigenval(i1) << " and " << eigenval(i2) << std::endl;
        if (std::abs(eigenval(i1) - std::conj(eigenval(i2))) > 1e-5) {
            std::cout << "Warning: Above eigenvalues not pairwise" << std::endl;
        }
        sc.theta[i] = std::abs(eigenval(i1));

        b(2 * i, 2 * i + 1) = -1.0;
        b(2 * i + 1, 2 * i) = +1.0;

        sc.coeff_1[i] = Q * b * Q.transpose();
        sc.coeff_2[i] = Q * b * b * Q.transpose();

        std::cout << "theta: "<< std::endl << sc.theta[i] << std::endl;
        std::cout << "b = " << std::endl << b << std::endl;
    }
}


template<unsigned int num_lvl>
solver_mpich_clvl_os_red<num_lvl>::solver_mpich_clvl_os_red
(std::shared_ptr<const device> dev, std::shared_ptr<scenario> scen) :
    solver_int(dev, scen),
    m_name("mpich-" + std::to_string(num_lvl) + "lvl-os-red")
{
    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    Eigen::initParallel();
    Eigen::setNbThreads(1);

    if (dev->get_regions().size() == 0) {
        throw std::invalid_argument("No regions in device!");
    }

    /* determine simulation settings */
    init_fdtd_simulation(dev, scen, 0.5);

    setup_generators();

    /* set up simulaton constants */
    std::map<std::string, unsigned int> id_to_idx;
    unsigned int j = 0;

    for (const auto& mat_id : dev->get_used_materials()) {
        sim_constants_clvl_os<num_lvl> sc;

        auto mat = material::get_from_library(mat_id);

        /* factor for electric field update */
        sc.M_CE = scen->get_timestep_size()/
            (EPS0 * mat->get_rel_permittivity());

        /* factor for magnetic field update */
        sc.M_CH = scen->get_timestep_size()/
            (MU0 * mat->get_rel_permeability() * scen->get_gridpoint_size());

        /* convert loss term to conductivity */
        sc.sigma = sqrt(EPS0 * mat->get_rel_permittivity()/
                        (MU0 * mat->get_rel_permeability()))
            * mat->get_losses() * 2.0;

        /* 3-lvl quantum mechanical system */
        /* active region in 3-lvl description? */
        /* TODO: remove ugly dynamic cast to qm_desc_3lvl, allow other
         * representations? */
        std::shared_ptr<qm_desc_clvl<num_lvl> > qm =
            std::dynamic_pointer_cast<qm_desc_clvl<num_lvl> >(mat->get_qm());
        if (qm) {
            /* factor for macroscopic polarization */
            sc.M_CP = 0.5 * mat->get_overlap_factor() *
                qm->get_carrier_density();

            /* determine dipole operator as vector */
            sc.v = get_adj_op(qm->get_dipole_op());

            std::cout << "v: " << std::endl << sc.v << std::endl;

            /* time-independent hamiltionian in adjoint representation */
            Eigen::Matrix<real, num_adj, num_adj> M_0;
            M_0 = get_adj_liouvillian(qm->get_hamiltonian());

            /* determine lindblad term in adjoint representation */
            Eigen::Matrix<real, num_adj, num_adj> G;
            G = get_adj_sop(qm->get_lindblad_op());

            /* time-independent part */
            Eigen::Matrix<real, num_adj, num_adj> M = M_0 + G;
            std::cout << "M: " << std::endl << M << std::endl;

            /* determine equilibrium term */
            Eigen::Matrix<real, num_adj, 1> d_eq;
            d_eq = get_adj_deq(qm->get_lindblad_op());
            sc.d_eq = d_eq;

            /* determine inhomogeneous term */
            real eps = std::numeric_limits<real>::epsilon();
            if (d_eq.isZero(eps)) {
                sc.d_in = Eigen::Matrix<real, num_adj, 1>::Zero();
            } else {
                /* solve equation system M * d_in = d_eq */
                sc.d_in = M.fullPivLu().solve(d_eq);

                real err = (M * sc.d_in - d_eq).norm() / d_eq.norm();
                std::cout << "d_in solver error: " << err << std::endl;
                /* TODO: throw exception or only report warning? */
                if (err > 1e-3) {
                    throw std::invalid_argument("Time-indepent matrix not "
                                                "invertible!");
                }
            }
            std::cout << "d_in: " << std::endl << sc.d_in << std::endl;

            /* determine constant propagator */
            Eigen::Matrix<real, num_adj, num_adj> A_0;
            A_0 = (M * scen->get_timestep_size()/2).exp();

            /* determine dipole operator in adjoint representation */
            Eigen::Matrix<real, num_adj, num_adj> U;
            U = get_adj_liouvillian(-qm->get_dipole_op());
            std::cout << "U: " << std::endl << U << std::endl;

            /* diagonalize dipole operator */
            /* note: returned eigenvectors are normalized */
            Eigen::EigenSolver<Eigen::Matrix<real, num_adj, num_adj> > es(U);

            /* store propagators B1 and B2 */
            //sc.B_1 = A_0 * es.eigenvectors();
            //sc.B_2 = es.eigenvectors().adjoint() * A_0;

            sc.A_0 = A_0;
            sc.B = es.eigenvectors();

            sc.M = M;
            sc.U = U;

            /* for Rodrigues formula */
            sc.U2 = U * U;
            sc.theta_1 = sqrt(pow(U(0, 1), 2) + pow(U(0, 2), 2)
                              + pow(U(1, 2), 2));

            /* for general analytic approach */
            fill_rodr_coeff<num_lvl, num_adj>(es.eigenvectors(),
                                              es.eigenvalues(), sc);

            /* TODO refine check? */
            sc.has_qm = true;
            sc.has_dipole = true;

            /* store diagonal matrix containing the eigenvalues */
            sc.L = es.eigenvalues() * scen->get_timestep_size();

            sc.d_init = get_adj_op(qm->get_d_init());

            std::cout << "init: " << sc.d_init << std::endl;

            /* TODO remove?
            if (scen->get_dm_init_type() == scenario::lower_full) {

            } else if (scen->get_dm_init_type() == scenario::upper_full) {
                sc.d_init = Eigen::Matrix<real, num_adj, 1>::Zero();
                sc.d_init(num_adj - 1) = 1;
            } else {
            }
            */
        } else {
            /* set all qm-related factors to zero */
            sc.M_CP = 0.0;

            sc.has_qm = false;
            sc.has_dipole = false;

            sc.v = Eigen::Matrix<real, num_adj, 1>::Zero();

            //sc.B_1 = Eigen::Matrix<complex, num_adj, num_adj>::Zero();
            //sc.B_2 = Eigen::Matrix<complex, num_adj, num_adj>::Zero();
            sc.A_0 = Eigen::Matrix<real, num_adj, num_adj>::Zero();
            sc.B = Eigen::Matrix<complex, num_adj, num_adj>::Zero();

            sc.M = Eigen::Matrix<real, num_adj, num_adj>::Zero();
            sc.U = Eigen::Matrix<real, num_adj, num_adj>::Zero();

            sc.L = Eigen::Matrix<complex, num_adj, 1>::Zero();

            sc.d_in = Eigen::Matrix<real, num_adj, 1>::Zero();
            sc.d_eq = Eigen::Matrix<real, num_adj, 1>::Zero();
            sc.d_init = Eigen::Matrix<real, num_adj, 1>::Zero();
        }

        /* simulation settings */
        sc.d_x_inv = 1.0/scen->get_gridpoint_size();
        sc.d_t = scen->get_timestep_size();

        m_sim_consts.push_back(sc);
        id_to_idx[mat->get_id()] = j;
        j++;
    }

    /* set up indices array and initialize data arrays */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::cout << "Number of threads: " << world_size << std::endl;

    unsigned int *l_mat_indices = new unsigned int[scen->get_num_gridpoints()];

    for (unsigned int i = 0; i < scen->get_num_gridpoints(); i++) {
        unsigned int mat_idx = 0;
        real x = i * scen->get_gridpoint_size();

        for (const auto& reg : dev->get_regions()) {
            if ((x >= reg->get_start()) && (x <= reg->get_end())) {
                mat_idx = id_to_idx[reg->get_material()->get_id()];
                break;
            }
        }
        l_mat_indices[i] = mat_idx;
    }


    /* set up results and transfer data structures */
    uint64_t scratch_size = 0;
    for (const auto& rec : scen->get_records()) {
        /* create copy list entry */
        copy_list_entry entry(rec, scen, scratch_size);

        std::cout << "Rows: " << entry.get_rows() << " Cols: " << entry.get_cols() << std::endl;

        /* MASTER ONLY: add result to solver */
        if (world_rank == 0) {
            m_results.push_back(entry.get_result());
        }

        /* calculate scratch size */
        scratch_size += entry.get_size();

        /* take imaginary part into account */
        if (rec->is_complex()) {
            scratch_size += entry.get_size();
        }

        /* TODO check if result is available */
        /*
           throw std::invalid_argument("Requested result is not available!");
        */

        m_copy_list.push_back(entry);
    }

    /* MASTER ONLY */
    if (world_rank == 0) {
        /* allocate scratchpad result memory */
        m_result_scratch = (real *) mb_aligned_alloc(sizeof(real) * scratch_size);
        m_scratch_size = scratch_size;
    }

    /* create source data */
    m_source_data = new real[scen->get_num_timesteps() *
                             scen->get_sources().size()];
    unsigned int base_idx = 0;
    for (const auto& src : scen->get_sources()) {
        sim_source s;
        s.type = src->get_type();
        s.x_idx = src->get_position()/scen->get_gridpoint_size();
        s.data_base_idx = base_idx;
        m_sim_sources.push_back(s);

        /* calculate source values */
        for (unsigned int j = 0; j < scen->get_num_timesteps(); j++) {
            m_source_data[base_idx + j] =
                src->get_value(j * scen->get_timestep_size());
        }

        base_idx += scen->get_num_timesteps();
    }

    uint64_t num_gridpoints = m_scenario->get_num_gridpoints();
    uint64_t chunk_base = m_scenario->get_num_gridpoints() / world_size;
    uint64_t chunk_rem = m_scenario->get_num_gridpoints() % world_size;
    uint64_t num_timesteps = m_scenario->get_num_timesteps();


    l_copy_list = m_copy_list.data();
    l_sim_consts = m_sim_consts.data();
    l_sim_sources = m_sim_sources.data();


    unsigned int chunk = chunk_base;

    if (world_rank == world_size - 1) {
        chunk += chunk_rem;
    }

    t_sync_record_size = m_copy_list.size() * OL * (chunk_base + chunk_rem);
    t_sync_record = new sync_record[t_sync_record_size + 1];


    /* allocation */
    uint64_t size = chunk + 2 * OL;

    t_d = (Eigen::Matrix<real, num_adj, 1> *)
        mb_aligned_alloc(size *
                         sizeof(Eigen::Matrix<real, num_adj, 1>));
    t_h = (real *) mb_aligned_alloc(size * sizeof(real));
    t_e = (real *) mb_aligned_alloc(size * sizeof(real));
    t_p = (real *) mb_aligned_alloc(size * sizeof(real));
    t_mat_indices = (unsigned int *) mb_aligned_alloc(size * sizeof(unsigned int));

    __mb_assume_aligned(t_d);
    __mb_assume_aligned(t_e);
    __mb_assume_aligned(t_p);
    __mb_assume_aligned(t_h);
    __mb_assume_aligned(t_mat_indices);

    for (int i = 0; i < size; i++) {
        uint64_t global_idx = world_rank * chunk_base + (i - OL);
        if ((global_idx >= 0) && (global_idx < num_gridpoints)) {
            unsigned int mat_idx = l_mat_indices[global_idx];
            t_mat_indices[i] = mat_idx;

            t_d[i] = l_sim_consts[mat_idx].d_init;
        } else {
            t_mat_indices[i] = 0;

            t_d[i] = Eigen::Matrix<real, num_adj, 1>::Zero();
        }
        t_e[i] = 0.0;
        t_p[i] = 0.0;
        t_h[i] = 0.0;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    delete[] l_mat_indices;
}

template<unsigned int num_lvl>
solver_mpich_clvl_os_red<num_lvl>::~solver_mpich_clvl_os_red()
{
    unsigned int num_sources = m_sim_sources.size();
    unsigned int num_copy = m_copy_list.size();
    uint64_t num_gridpoints = m_scenario->get_num_gridpoints();
    uint64_t num_timesteps = m_scenario->get_num_timesteps();

    mb_aligned_free(t_h);
    mb_aligned_free(t_e);
    mb_aligned_free(t_p);
    mb_aligned_free(t_d);
    mb_aligned_free(t_mat_indices);

    /* MASTER ONLY */
    if (world_rank == 0) {
        mb_aligned_free(m_result_scratch);
    }
    delete[] m_source_data;
    delete[] t_sync_record;

    // Finalize the MPI environment.
    MPI_Finalize();
}

template<unsigned int num_lvl>
const std::string&
solver_mpich_clvl_os_red<num_lvl>::get_name() const
{
    return m_name;
}

template<unsigned int num_lvl, unsigned int num_adj>
void
update_fdtd(uint64_t size, unsigned int border, real *t_e, real *t_p,
            real *t_h, Eigen::Matrix<real, num_adj, 1>* t_d,
            unsigned int *t_mat_indices,
            sim_constants_clvl_os<num_lvl> *l_sim_consts)
{
    for (int i = border; i < size - border - 1; i++) {
        int mat_idx = t_mat_indices[i];

        real j = l_sim_consts[mat_idx].sigma * t_e[i];

        t_e[i] += l_sim_consts[mat_idx].M_CE *
            (-j - t_p[i] + (t_h[i + 1] - t_h[i]) *
             l_sim_consts[mat_idx].d_x_inv);
        /*
        if (i >= border + 1) {
            t_h[i] += l_sim_consts[mat_idx].M_CH * (t_e[i] - t_e[i - 1]);
        }
        */
    }
}

template<unsigned int num_lvl, unsigned int num_adj>
void
update_h(unsigned int size, unsigned int border, real *t_e, real *t_p,
            real *t_h, Eigen::Matrix<real, num_adj, 1>* t_d,
            unsigned int *t_mat_indices,
            sim_constants_clvl_os<num_lvl> *l_sim_consts)
{
    for (int i = border; i < size - border - 1; i++) {
        int mat_idx = t_mat_indices[i];

        if (i >= border + 1) {
            t_h[i] += l_sim_consts[mat_idx].M_CH * (t_e[i] - t_e[i - 1]);
        }
    }
}

void
apply_sources(real *t_e, real *source_data, unsigned int num_sources,
              sim_source *l_sim_sources, uint64_t time,
              unsigned int base_pos, uint64_t chunk)
{
    for (unsigned int k = 0; k < num_sources; k++) {
        int at = l_sim_sources[k].x_idx - base_pos + OL;
        if ((at > 0) && (at < chunk + 2 * OL)) {
            real src = source_data[l_sim_sources[k].data_base_idx + time];
            if (l_sim_sources[k].type == source::type::hard_source) {
                t_e[at] = src;
            } else if (l_sim_sources[k].type == source::type::soft_source) {
                /* TODO: fix source */
                t_e[at] += src;
            } else {
            }
        }
    }
}

complex mexp(const complex& arg)
{
    return std::exp(arg);
}

template<unsigned int num_lvl, unsigned int num_adj>
inline Eigen::Matrix<real, num_adj, num_adj>
mat_exp(const sim_constants_clvl_os<num_lvl>& s, real e)
{
    Eigen::Matrix<real, num_adj, num_adj> ret;

#if EXP_METHOD==1
    /* by diagonalization */
    Eigen::Matrix<complex, num_adj, 1> diag_exp = s.L * e;
    diag_exp = diag_exp.unaryExpr(&mexp);

    ret = (s.B * diag_exp.asDiagonal() * s.B.adjoint()).real();
#elif EXP_METHOD==2
    /* analytic solution */
    if (num_lvl == 2) {
        /* Rodrigues formula */
        ret = sin(s.theta_1 * e * s.d_t)/s.theta_1 * s.U +
            (1 - cos(s.theta_1 * e * s.d_t))/(s.theta_1 * s.theta_1) * s.U2 +
            Eigen::Matrix<real, num_adj, num_adj>::Identity();
    } else {
        ret = Eigen::Matrix<real, num_adj, num_adj>::Identity();
        for (int i = 0; i < num_adj/2; i++) {
            /* TODO nolias()? */
            ret += sin(s.theta[i] * e * s.d_t) * s.coeff_1[i] +
                (1 - cos(s.theta[i] * e * s.d_t)) * s.coeff_2[i];
        }
    }
#else
    /* Eigen matrix exponential */
    ret = (s.U * e * s.d_t).exp();
#endif

    return ret;
}

template<unsigned int num_lvl, unsigned int num_adj>
void
update_d(uint64_t size, unsigned int border, real *t_e, real *t_p,
         Eigen::Matrix<real, num_adj, 1>* t_d,
         unsigned int *t_mat_indices,
         sim_constants_clvl_os<num_lvl> *l_sim_consts)
{
    for (int i = border; i < size - border - 1; i++) {
        int mat_idx = t_mat_indices[i];

        if (l_sim_consts[mat_idx].has_qm) {
            /* update density matrix */
            Eigen::Matrix<real, num_adj, 1> d1, d2;

            /* time-indepedent half step */
            d1 = l_sim_consts[mat_idx].A_0 *
                (t_d[i] + l_sim_consts[mat_idx].d_in)
                - l_sim_consts[mat_idx].d_in;

            /* time-dependent full step */
            if (l_sim_consts[mat_idx].has_dipole) {
                /* determine time-dependent propagator */
                Eigen::Matrix<real, num_adj, num_adj> A_I =
                    mat_exp<num_lvl, num_adj>(l_sim_consts[mat_idx], t_e[i]);
                d2 = A_I * d1;
            } else {
                d2 = d1;
            }

            /* time-indepedent half step */
            t_d[i] = l_sim_consts[mat_idx].A_0 *
                (d2 + l_sim_consts[mat_idx].d_in)
                - l_sim_consts[mat_idx].d_in;

            /* update polarization */
            t_p[i] = l_sim_consts[mat_idx].M_CP *
                l_sim_consts[mat_idx].v.transpose() *
                (l_sim_consts[mat_idx].M * t_d[i] +
                 l_sim_consts[mat_idx].d_eq);
        } else {
            t_p[i] = 0.0;
        }
    }
}

template<unsigned int num_lvl>
inline void
solver_mpich_clvl_os_red<num_lvl>::write_sync_record(size_t offset, real value) const {
    size_t recorder_counter = t_sync_record[0].offset;

    if (recorder_counter >= t_sync_record_size) {
        throw std::invalid_argument("Too many results");
    } else {
        t_sync_record[1 + recorder_counter].offset = offset;
        t_sync_record[1 + recorder_counter].value = value;
        t_sync_record[0].offset += 1;
    }
};

template<unsigned int num_lvl>
void
solver_mpich_clvl_os_red<num_lvl>::run() const
{
    uint64_t num_gridpoints = m_scenario->get_num_gridpoints();
    uint64_t chunk_base = m_scenario->get_num_gridpoints() / world_size;
    uint64_t chunk_rem = m_scenario->get_num_gridpoints() % world_size;
    uint64_t num_timesteps = m_scenario->get_num_timesteps();
    unsigned int num_sources = m_sim_sources.size();
    unsigned int num_copy = m_copy_list.size();


    uint64_t chunk = chunk_base;
    if (world_rank == world_size - 1) {
        chunk += chunk_rem;
    }

    uint64_t size = chunk + 2 * OL;

    if (world_rank == 0) {
        __mb_assume_aligned(m_result_scratch);
    }

    /* main loop */
    for (uint64_t n = 0; n <= num_timesteps/OL; n++) {
        if (world_rank == 0) {
            std::cout << "[" << n << "/" << num_timesteps/OL << "]" << std::endl;
        }

        /* handle loop remainder */
        unsigned int subloop_ct = (n == num_timesteps/OL) ?
            num_timesteps % OL : OL;

        /* exchange data */
        /* 1. even send data to odd in the right */
        /* 2. even recv data from odd in the right */
        /* 3. even send data to odd in the left */
        /* 4. even recv data from odd in the left */

        sync_body<OL> body;
        /* 1. even send data to odd in the right */
        if (world_rank % 2 == 0) {
            if ((world_rank + 1) < world_size) {
                for (unsigned int i = 0; i < OL; i++) {
                    body.d[i] = t_d[OL + chunk_base + i];
                    body.e[i] = t_e[OL + chunk_base + i];
                    body.h[i] = t_h[OL + chunk_base + i];
                }
                MPI_Send(&body, sizeof(body), MPI_BYTE, world_rank + 1, 0, MPI_COMM_WORLD);
            }
        } else {
            if (world_rank > 0) {
                MPI_Recv(&body, sizeof(body), MPI_BYTE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[i] = body.d[i];
                    t_e[i] = body.e[i];
                    t_h[i] = body.h[i];
                }
            }
        }
        /* 2. even recv data from odd in the right */
        if (world_rank % 2 == 0) {
            if ((world_rank + 1) < world_size) {
                MPI_Recv(&body, sizeof(body), MPI_BYTE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[OL + chunk_base + i] = body.d[i];
                    t_e[OL + chunk_base + i] = body.e[i];
                    t_h[OL + chunk_base + i] = body.h[i];
                }
            }
        } else {
            if (world_rank > 0) {
                for (unsigned int i = 0; i < OL; i++) {
                    body.d[i] = t_d[i];
                    body.e[i] = t_e[i];
                    body.h[i] = t_h[i];
                }
                MPI_Send(&body, sizeof(body), MPI_BYTE, world_rank - 1, 0, MPI_COMM_WORLD);
            }
        }
        /* 3. even send data to odd in the left */
        if (world_rank % 2 == 0) {
            if (world_rank > 0) {
                for (unsigned int i = 0; i < OL; i++) {
                    body.d[i] = t_d[i];
                    body.e[i] = t_e[i];
                    body.h[i] = t_h[i];
                }
                MPI_Send(&body, sizeof(body), MPI_BYTE, world_rank - 1, 0, MPI_COMM_WORLD);
            }
        } else {
            if ((world_rank + 1) < world_size) {
                MPI_Recv(&body, sizeof(body), MPI_BYTE, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[OL + chunk_base + i] = body.d[i];
                    t_e[OL + chunk_base + i] = body.e[i];
                    t_h[OL + chunk_base + i] = body.h[i];
                }
            }
        }
        /* 4. even recv data from odd in the left */
        if (world_rank % 2 == 0) {
            if (world_rank > 0) {
                MPI_Recv(&body, sizeof(body), MPI_BYTE, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (unsigned int i = 0; i < OL; i++) {
                    t_d[i] = body.d[i];
                    t_e[i] = body.e[i];
                    t_h[i] = body.h[i];
                }
            }
        } else {
            if ((world_rank + 1) < world_size) {
                for (unsigned int i = 0; i < OL; i++) {
                    body.d[i] = t_d[OL + chunk_base + i];
                    body.e[i] = t_e[OL + chunk_base + i];
                    body.h[i] = t_h[OL + chunk_base + i];
                }
                MPI_Send(&body, sizeof(body), MPI_BYTE, world_rank + 1, 0, MPI_COMM_WORLD);
            }
        }

        /* sync after communication */
        MPI_Barrier(MPI_COMM_WORLD);

        /* sub-loop */
        /* clear sync records */
        t_sync_record[0].offset = 0;
        for (unsigned int m = 0; m < subloop_ct; m++) {
            /* align border to vector length */
            unsigned int border = m - (m % VEC);

            /* update d */
            update_d<num_lvl, num_adj>(size, border, t_e, t_p, t_d,
                                       t_mat_indices, l_sim_consts);

             /* update e + h with fdtd */
            update_fdtd<num_lvl, num_adj>(size, border, t_e, t_p, t_h,
                                          t_d, t_mat_indices,
                                          l_sim_consts);

            /* apply sources */
            apply_sources(t_e, m_source_data, num_sources,
                          l_sim_sources, n * OL + m, world_rank * chunk_base,
                          chunk);

            update_h<num_lvl, num_adj>(size, border, t_e, t_p, t_h,
                                       t_d, t_mat_indices,
                                       l_sim_consts);


            /* apply field boundary condition */
            if (world_rank == 0) {
                t_h[OL] = 0;
            }
            if (world_rank == world_size - 1) {
                t_h[OL + chunk] = 0;
            }

            /* save records to buffer */
            for (int k = 0; k < num_copy; k++) {
                if (l_copy_list[k].hasto_record(n * OL + m)) {
                    uint64_t pos = l_copy_list[k].get_position();
                    uint64_t cols = l_copy_list[k].get_cols();
                    uint64_t ridx = l_copy_list[k].get_row_idx();
                    uint64_t cidx = l_copy_list[k].get_col_idx();
                    record::type t = l_copy_list[k].get_type();

                    int64_t base_idx = world_rank * chunk_base - OL;
                    int64_t off_r = l_copy_list[k].get_offset_scratch_real
                        (n * OL + m, base_idx - pos);


                    for (uint64_t i = OL; i < chunk + OL; i++) {
                        int64_t idx = base_idx + i;
                        if ((idx >= pos) && (idx < pos + cols)) {
                            if (t == record::type::electric) {
                                write_sync_record(off_r + i, t_e[i]);
                            } else if (t == record::type::magnetic) {
                                write_sync_record(off_r + i, t_h[i]);
                            } else if (t == record::type::inversion) {
                                write_sync_record(off_r + i, t_d[i](num_lvl * (num_lvl - 1)));
                            } else if (t == record::type::density) {

                                /* right now only populations */
                                real temp = 1.0/num_lvl;
                                for (int l = num_lvl * (num_lvl - 1);
                                     l < num_adj; l++) {
                                    temp += 0.5 * t_d[i](l) *
                                        m_generators[l](ridx, cidx).real();
                                }
                                write_sync_record(off_r + i, temp);

                                /* TODO: coherences
                                 * remove 1/3
                                 * consider only two/one corresponding
                                 * entry */

                            } else {
                                /* TODO handle trouble */

                                /* TODO calculate regular populations */
                            }
                        }
                    }

                }
            }


        } /* end sub loop */

        /* sync records */
        if (world_rank == 0) {
            // save master's own records
            for (size_t i = 0; i < t_sync_record[0].offset; i++) {
                m_result_scratch[t_sync_record[1 + i].offset] = t_sync_record[1 + i].value;
            }

            // recv other process's records
            // std::cout << "receiving data..." << t_sync_record_size << std::endl;
            for (int i = 1; i < world_size; i++) {
                MPI_Recv(t_sync_record, sizeof(sync_record) * (1 + t_sync_record_size), MPI_BYTE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                // std::cout << "getting " << t_sync_record[0].offset << " records from process " << i << std::endl;
                for (size_t j = 0; j < t_sync_record[0].offset; j++) {
                    m_result_scratch[t_sync_record[1 + i].offset] = t_sync_record[1 + i].value;
                }
            }
        } else {
            // std::cout << "sending data..." << (1 + t_sync_record[0].offset) << std::endl;
            MPI_Send(t_sync_record, sizeof(sync_record) * (1 + t_sync_record[0].offset), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
        }

        /* sync after computation */
        MPI_Barrier(MPI_COMM_WORLD);
    } /* end main foor loop */


    /* bulk copy results into result classes */
    if (world_rank == 0) {
        for (const auto& cle : m_copy_list) {
            real *dr = m_result_scratch + cle.get_offset_scratch_real(0, 0);
            std::copy(dr, dr + cle.get_size(), cle.get_result_real(0, 0));
            if (cle.is_complex()) {
                real *di = m_result_scratch + cle.get_offset_scratch_imag(0, 0);
                std::copy(di, di + cle.get_size(), cle.get_result_imag(0, 0));
            }
        }
    }
}

}
