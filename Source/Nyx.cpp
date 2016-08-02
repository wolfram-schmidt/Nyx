#include <winstd.H>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include <unistd.h>

using std::cout;
using std::cerr;
using std::endl;
using std::istream;
using std::ostream;
using std::pair;
using std::string;

#include <Utility.H>
#include <CONSTANTS.H>
#include <Nyx.H>
#include <Nyx_F.H>
#include <Derive_F.H>
#include <VisMF.H>
#include <TagBox.H>
#include <Particles_F.H>

#if BL_USE_MPI
#include "MemInfo.H"
#endif

#ifdef GRAVITY
#include "Gravity.H"
#endif

#ifdef IN_SITU
#include <boxlib_in_situ_analysis.H>
#endif

static int sum_interval = -1;
static Real fixed_dt    = -1.0;
static Real initial_dt  = -1.0;
static Real dt_cutoff   =  0;

int Nyx::strict_subcycling = 0;

Real Nyx::old_a      = -1.0;
Real Nyx::new_a      = -1.0;
Real Nyx::old_a_time = -1.0;
Real Nyx::new_a_time = -1.0;

Array<Real> Nyx::plot_z_values;

bool Nyx::dump_old = false;
int Nyx::verbose      = 0;
int Nyx::show_timings = 0;

Real Nyx::cfl = 0.8;
Real Nyx::init_shrink = 1.0;
Real Nyx::change_max  = 1.1;

BCRec Nyx::phys_bc;
int Nyx::do_reflux = 1;
int Nyx::NUM_STATE = -1;
int Nyx::NUM_GROW = -1;

int Nyx::nsteps_from_plotfile = -1;

ErrorList Nyx::err_list;

int Nyx::Density = -1;
int Nyx::Eden = -1;
int Nyx::Eint = -1;
int Nyx::Xmom = -1;
int Nyx::Ymom = -1;
int Nyx::Zmom = -1;

int Nyx::Temp_comp = -1;
int Nyx::  Ne_comp = -1;

int Nyx::NumSpec  = 0;
int Nyx::NumAux   = 0;
int Nyx::NumAdv   = 0;

int Nyx::FirstSpec = -1;
int Nyx::FirstAux  = -1;
int Nyx::FirstAdv  = -1;

Real Nyx::small_dens = -1.e200;
Real Nyx::small_temp = -1.e200;
Real Nyx::gamma      =  0;

int Nyx::do_hydro = -1;
int Nyx::add_ext_src = 0;
int Nyx::heat_cool_type = 0;
int Nyx::strang_split = 0;

Real Nyx::average_gas_density = 0;
Real Nyx::average_dm_density = 0;
Real Nyx::average_neutr_density = 0;
Real Nyx::average_total_density = 0;

// Real Nyx::ave_lev_vorticity[10];
// Real Nyx::std_lev_vorticity[10];

#ifdef GRAVITY
Gravity* Nyx::gravity  =  0;
int Nyx::do_grav       = -1;
#else
int Nyx::do_grav       =  0;
#endif

int Nyx::allow_untagging    = 0;
int Nyx::use_const_species  = 0;
int Nyx::normalize_species  = 0;
int Nyx::ppm_type           = 1;
int Nyx::ppm_reference      = 1;
int Nyx::corner_coupling    = 1;

int Nyx::use_colglaz        = 0;
int Nyx::version_2          = 0;

int Nyx::use_flattening     = 1;
int Nyx::ppm_flatten_before_integrals = 0;

Real Nyx:: h_species        = 0.0;
Real Nyx::he_species        = 0.0;

int Nyx::use_exact_gravity  = 0;

#ifdef _OPENMP
#include <omp.h>
#endif

int Nyx::write_parameters_in_plotfile = true;
int Nyx::print_fortran_warnings       = true;

// Do we use separate SPH particles to initialize
//  the density and momentum on the grid?
int  Nyx::init_with_sph_particles = 0;

// Do we write the particles in single (IEEE32)
//  or doublue (NATIVE) precision?
#ifdef BL_SINGLE_PRECISION_PARTICLES
std::string Nyx::particle_plotfile_format = "IEEE32";
#else
std::string Nyx::particle_plotfile_format = "NATIVE";
#endif

// Note: Nyx::variableSetUp is in Nyx_setup.cpp
void
Nyx::variable_cleanup ()
{
#ifdef GRAVITY
    if (gravity != 0)
    {
        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "Deleting gravity in variable_cleanup...\n";
        delete gravity;
        gravity = 0;
    }
#endif

    desc_lst.clear();
}

void
Nyx::read_params ()
{
    static bool done = false;

    if (done) return;  // (caseywstark) when would this happen?

    done = true;  // ?

    ParmParse pp("nyx");

    pp.query("v", verbose);
    pp.query("show_timings", show_timings);
    //verbose = (verbose ? 1 : 0);
    pp.get("init_shrink", init_shrink);
    pp.get("cfl", cfl);
    pp.query("change_max", change_max);
    pp.query("fixed_dt", fixed_dt);
    pp.query("initial_dt", initial_dt);
    pp.query("sum_interval", sum_interval);
    pp.query("do_reflux", do_reflux);
    do_reflux = (do_reflux ? 1 : 0);
    pp.get("dt_cutoff", dt_cutoff);

    pp.query("dump_old", dump_old);

    pp.query("small_dens", small_dens);
    pp.query("small_temp", small_temp);
    pp.query("gamma", gamma);
    
    pp.query("strict_subcycling",strict_subcycling);

    // Get boundary conditions
    Array<int> lo_bc(BL_SPACEDIM), hi_bc(BL_SPACEDIM);
    pp.getarr("lo_bc", lo_bc, 0, BL_SPACEDIM);
    pp.getarr("hi_bc", hi_bc, 0, BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        phys_bc.setLo(i, lo_bc[i]);
        phys_bc.setHi(i, hi_bc[i]);
    }

    //
    // Check phys_bc against possible periodic geometry
    // if periodic, must have internal BC marked.
    //
    if (Geometry::isAnyPeriodic())
    {
        //
        // Do idiot check.  Periodic means interior in those directions.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (Geometry::isPeriodic(dir))
            {
                if (lo_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but low BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
                if (hi_bc[dir] != Interior)
                {
                    std::cerr << "Nyx::read_params:periodic in direction "
                              << dir
                              << " but high BC is not Interior" << std::endl;
                    BoxLib::Error();
                }
            }
        }
    }
    else
    {
        //
        // Do idiot check.  If not periodic, should be no interior.
        //
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            if (lo_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
            if (hi_bc[dir] == Interior)
            {
                std::cerr << "Nyx::read_params:interior bc in direction "
                          << dir
                          << " but not periodic" << std::endl;
                BoxLib::Error();
            }
        }
    }

    pp.get("do_hydro", do_hydro);
#ifdef NO_HYDRO
    if (do_hydro == 1) 
        BoxLib::Error("Cant have do_hydro == 1 when NO_HYDRO is true");
#endif

#ifdef NO_HYDRO
#ifndef GRAVITY
        BoxLib::Error("Dont know what to do with both hydro and gravity off");
#endif
#endif

    pp.query("add_ext_src", add_ext_src);
    pp.query("strang_split", strang_split);

    pp.query("heat_cool_type", heat_cool_type);

    pp.query("use_exact_gravity", use_exact_gravity);

#ifdef HEATCOOL
    if (heat_cool_type > 0 && add_ext_src == 0)
       BoxLib::Error("Nyx::must set add_ext_src to 1 if heat_cool_type > 0");
    if (heat_cool_type != 1 && heat_cool_type != 3)
       BoxLib::Error("Nyx:: nonzero heat_cool_type must equal 1 or 3"); 
    if (heat_cool_type == 0)
       BoxLib::Error("Nyx::contradiction -- HEATCOOL is defined but heat_cool_type == 0");
#else
    if (heat_cool_type > 0)
       BoxLib::Error("Nyx::you set heat_cool_type > 0 but forgot to set USE_HEATCOOL = TRUE");
#endif

    pp.query("allow_untagging", allow_untagging);
    pp.query("use_const_species", use_const_species);
    pp.query("normalize_species", normalize_species);
    pp.query("ppm_type", ppm_type);
    pp.query("ppm_reference", ppm_reference);
    pp.query("ppm_flatten_before_integrals", ppm_flatten_before_integrals);
    pp.query("use_flattening", use_flattening);
    pp.query("use_colglaz", use_colglaz);
    pp.query("version_2", version_2);
    pp.query("corner_coupling", corner_coupling);

    if (do_hydro == 1)
    {
        if (do_hydro == 1 && use_const_species == 1) 
        {
           pp.get("h_species" ,  h_species);
           pp.get("he_species", he_species);
           BL_FORT_PROC_CALL(SET_XHYDROGEN,set_xhydrogen)(h_species);
           if (ParallelDescriptor::IOProcessor())
           {
               std::cout << "Nyx::setting species concentrations to " 
                         << h_species << " and " << he_species
                         << " in the hydro and in the EOS " << std::endl;
           }
        }

        // 
        if (use_colglaz == 1)
        {
           if (ppm_type == 0 && ParallelDescriptor::IOProcessor())
               std::cout << "Nyx::setting use_colglaz = 1 with ppm_type = 0 \n";
           if (ppm_type != 0)
               BoxLib::Error("Nyx::ppm_type must be 0 with use_colglaz = 1");
        }

        // ppm_flatten_before_integrals is only done for ppm_type != 0
        if (ppm_type == 0 && ppm_flatten_before_integrals > 0)
        {
            std::cerr << "ppm_flatten_before_integrals > 0 not implemented for ppm_type != 0 \n";
            BoxLib::Error();
        }

        if (version_2 > 0 && ppm_type == 0)
           BoxLib::Error("Nyx::version_2 only defined for ppm_type = 1");

        if (version_2 !=0 && version_2 != 1 && version_2 != 2 && version_2 != 3)
           BoxLib::Error("Nyx:: don't know what to do with version_2 flag");

        // Make sure ppm_type is set correctly.
        if (ppm_type != 0 && ppm_type != 1 && ppm_type != 2) 
        {
           BoxLib::Error("Nyx::ppm_type must be 0, 1 or 2");
        }
    }
    
    // Make sure not to call refluxing if we're not actually doing any hydro.
    if (do_hydro == 0) do_reflux = 0;

#ifdef GRAVITY
    pp.get("do_grav", do_grav);
#endif

    read_particle_params();

    read_init_params();

    pp.query("write_parameter_file",write_parameters_in_plotfile);
    pp.query("print_fortran_warnings",print_fortran_warnings);

    read_comoving_params();

    if (pp.contains("plot_z_values"))
    {
      int num_z_values = pp.countval("plot_z_values");
      plot_z_values.resize(num_z_values);
      pp.queryarr("plot_z_values",plot_z_values,0,num_z_values);
    }
}

Nyx::Nyx ()
{
#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
    }
#endif
    fine_mask = 0;
}

Nyx::Nyx (Amr&            papa,
          int             lev,
          const Geometry& level_geom,
          const BoxArray& bl,
          Real            time)
    :
    AmrLevel(papa,lev,level_geom,bl,time)
{
    build_metrics();
    fine_mask = 0;

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        flux_reg = 0;
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    // Initialize to zero here in case we run with do_grav = false.
    MultiFab& new_grav_mf = get_new_data(Gravity_Type);
    new_grav_mf.setVal(0);

    if (do_grav)
    {
        // gravity is a static object, only alloc if not already there
        if (gravity == 0)
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);

        gravity->install_level(level, this);

        if (verbose && level == 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Setting the gravity type to "
                      << gravity->get_gravity_type() << '\n';
   }
#endif

    // Initialize the "a" variable
    if (level == 0 && time == 0.0 && old_a_time < 0.)
    {
       old_a_time = 0.0;
       new_a_time = 0.0;

       old_a = 1.0 / (1.0 + initial_z);
       new_a = old_a;
    }

     // Initialize "this_z" in the atomic_rates_module 
     if (heat_cool_type == 1 || heat_cool_type == 3)
         BL_FORT_PROC_CALL(INIT_THIS_Z, init_this_z)(&old_a);
    
    // Set grav_n_grow to 3 on init. It'll be reset in advance.
    grav_n_grow = 3;
}

Nyx::~Nyx ()
{
#ifndef NO_HYDRO
    if (do_hydro == 1)
        delete flux_reg;
#endif
    delete fine_mask;
}

void
Nyx::restart (Amr&     papa,
              istream& is,
              bool     b_read_special)
{
    AmrLevel::restart(papa, is, b_read_special);

    build_metrics();

#ifndef NO_HYDRO
    if (do_hydro == 1)
    {
        BL_ASSERT(flux_reg == 0);
        if (level > 0 && do_reflux)
            flux_reg = new FluxRegister(grids, crse_ratio, level, NUM_STATE);
    }
#endif

#ifdef GRAVITY
    if (do_grav && level == 0)
    {
        BL_ASSERT(gravity == 0);
        gravity = new Gravity(parent, parent->finestLevel(), &phys_bc, Density);
    }
#endif
}

void
Nyx::build_metrics ()
{
}

void
Nyx::setTimeLevel (Real time,
                   Real dt_old,
                   Real dt_new)
{
    if (ParallelDescriptor::IOProcessor()) {
       std::cout << "Setting the current time in the state data to " 
                 << parent->cumTime() << std::endl;
    }
    AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

void
Nyx::init (AmrLevel& old)
{
    Nyx* old_level = (Nyx*) &old;
    //
    // Create new grid data by fillpatching from old.
    //
    Real dt_new = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time = old_level->state[PhiGrav_Type].curTime();
    Real prev_time = old_level->state[PhiGrav_Type].prevTime();
#else
    Real cur_time = old_level->state[State_Type].curTime();
    Real prev_time = old_level->state[State_Type].prevTime();
#endif

    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

#ifndef NO_HYDRO
    if (do_hydro == 1) 
    {
        MultiFab& S_new = get_new_data(State_Type);
        MultiFab& D_new = get_new_data(DiagEOS_Type);

        for (FillPatchIterator
                 fpi(old, S_new, 0, cur_time,   State_Type, 0, NUM_STATE),
                dfpi(old, D_new, 0, cur_time, DiagEOS_Type, 0, 2);
                fpi.isValid() && dfpi.isValid();
                ++fpi,++dfpi)
        {
            FArrayBox&  tmp =  fpi();
            FArrayBox& dtmp = dfpi();
            S_new[fpi].copy(tmp);
            D_new[fpi].copy(dtmp);
        }
    }
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    for (FillPatchIterator fpi(old, Phi_new, 0, cur_time, PhiGrav_Type, 0, 1);
         fpi.isValid(); ++fpi)
    {
        Phi_new[fpi].copy(fpi());
    }
#endif

    // Set E in terms of e + kinetic energy
    // if (do_hydro)
    // enforce_consistent_e(S_new);
}

//
// This version inits the data on a new level that did not
// exist before regridding.
//
void
Nyx::init ()
{
    Real dt        = parent->dtLevel(level);
#ifdef NO_HYDRO
    Real cur_time  = get_level(level-1).state[PhiGrav_Type].curTime();
    Real prev_time = get_level(level-1).state[PhiGrav_Type].prevTime();
#else
    Real cur_time  = get_level(level-1).state[State_Type].curTime();
    Real prev_time = get_level(level-1).state[State_Type].prevTime();
#endif
    Real dt_old    = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level-1);

    setTimeLevel(cur_time, dt_old, dt);

#ifndef NO_HYDRO
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);
    FillCoarsePatch(S_new, 0, cur_time,   State_Type, 0, S_new.nComp());
    FillCoarsePatch(D_new, 0, cur_time, DiagEOS_Type, 0, D_new.nComp());
#endif

#ifdef GRAVITY
    MultiFab& Phi_new = get_new_data(PhiGrav_Type);
    FillCoarsePatch(Phi_new, 0, cur_time, PhiGrav_Type, 0, Phi_new.nComp());
#endif

    // We set dt to be large for this new level to avoid screwing up
    // computeNewDt.
    parent->setDtLevel(1.e100, level);

    // Set E in terms of e + kinetic energy
    // if (do_hydro)
    // enforce_consistent_e(S_new);
}

Real
Nyx::initial_time_step ()
{
    Real dummy_dt = 0;
    Real init_dt = 0;

    if (initial_dt > 0)
    {
        init_dt = initial_dt;
    }
    else
    {
        init_dt = init_shrink * est_time_step(dummy_dt);
    }

    bool dt_changed = false;
    if (level == 0 && plot_z_values.size() > 0)
        plot_z_est_time_step(init_dt,dt_changed);

    return init_dt;
}

Real
Nyx::est_time_step (Real dt_old)
{
    if (fixed_dt > 0)
        return fixed_dt;

    // This is just a dummy value to start with
    Real est_dt = 1.0e+200;

#ifndef NO_HYDRO
    const MultiFab& stateMF = get_new_data(State_Type);
#endif

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

#ifndef NO_HYDRO
    if (do_hydro)
    {
        Real a = get_comoving_a(cur_time);
        const Real* dx = geom.CellSize();

	Real dt = est_dt;
	  
#ifdef _OPENMP
#pragma omp parallel firstprivate(dt)
#endif
	{
	  for (MFIter mfi(stateMF,true); mfi.isValid(); ++mfi)
	    {
	      const Box& box = mfi.tilebox();

	      BL_FORT_PROC_CALL(FORT_ESTDT, fort_estdt)
                (BL_TO_FORTRAN(stateMF[mfi]), box.loVect(), box.hiVect(), dx,
                 &dt, &a);
	    }
#ifdef _OPENMP
#pragma omp critical (nyx_estdt)	      
#endif
	  {
	    est_dt = std::min(est_dt, dt);
	  }
	}

        // If in comoving coordinates, then scale dt (based on u and c) by a
        est_dt *= a;

        ParallelDescriptor::ReduceRealMin(est_dt);
        est_dt *= cfl;
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "...estdt from hydro at level "
                      << level << ": "
                      << est_dt << '\n';
    }
#endif

#ifdef GRAVITY
    particle_est_time_step(est_dt);
#endif

    if (level==0)
        comoving_est_time_step(cur_time,est_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Nyx::est_time_step at level "
                  << level
                  << ":  estdt = "
                  << est_dt << '\n';

    return est_dt;
}

void
Nyx::computeNewDt (int                   finest_level,
                   int                   sub_cycle,
                   Array<int>&           n_cycle,
                   const Array<IntVect>& ref_ratio,
                   Array<Real>&          dt_min,
                   Array<Real>&          dt_level,
                   Real                  stop_time,
                   int                   post_regrid_flag)
{
    //
    // We are at the start of a coarse grid timecycle.
    // Compute the timesteps for the next iteration.
    //
    if (level > 0)
        return;
    
    int i;

    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        Nyx& adv_level = get_level(i);
        dt_min[i] = adv_level.est_time_step(dt_level[i]);
    }

    if (fixed_dt <= 0.0)
    {
        if (post_regrid_flag == 1)
        {
            //
            // Limit dt's by pre-regrid dt
            //
            for (i = 0; i <= finest_level; i++)
            {
                dt_min[i] = std::min(dt_min[i], dt_level[i]);
            }
            //
            // Find the minimum over all levels
            //
            for (i = 0; i <= finest_level; i++)
            {
                n_factor *= n_cycle[i];
                dt_0 = std::min(dt_0, n_factor * dt_min[i]);
            }
        }
        else
        {
            bool sub_unchanged=true;
            if ((parent->maxLevel() > 0) && (level == 0) &&
                (parent->subcyclingMode() == "Optimal") && 
                (parent->okToRegrid(level) || parent->levelSteps(0) == 0) )
            {
                int new_cycle[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                    new_cycle[i] = n_cycle[i];
                // The max allowable dt
                Real dt_max[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    dt_max[i] = dt_min[i];
                }
                // find the maximum number of cycles allowed.
                int cycle_max[finest_level+1];
                cycle_max[0] = 1;
                for (i = 1; i <= finest_level; i++)
                {
                    cycle_max[i] = parent->MaxRefRatio(i-1);
                }
                // estimate the amout of work to advance each level.
                Real est_work[finest_level+1];
                for (i = 0; i <= finest_level; i++)
                {
                    est_work[i] = parent->getLevel(i).estimateWork();
                }
                // this value will be used only if the subcycling pattern is changed.
                dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
                for (i = 0; i <= finest_level; i++)
                {
                    if (n_cycle[i] != new_cycle[i])
                    {
                        sub_unchanged = false;
                        n_cycle[i] = new_cycle[i];
                    }
                }
                
            }
            
            if (sub_unchanged)
            //
            // Limit dt's by change_max * old dt
            //
            {
                for (i = 0; i <= finest_level; i++)
                {
                    if (verbose && ParallelDescriptor::IOProcessor())
                    {
                        if (dt_min[i] > change_max*dt_level[i])
                        {
                            std::cout << "Nyx::compute_new_dt : limiting dt at level "
                                      << i << '\n';
                            std::cout << " ... new dt computed: " << dt_min[i]
                                      << '\n';
                            std::cout << " ... but limiting to: "
                                      << change_max * dt_level[i] << " = " << change_max
                                      << " * " << dt_level[i] << '\n';
                        }
                    }

                    dt_min[i] = std::min(dt_min[i], change_max * dt_level[i]);
                }
                //
                // Find the minimum over all levels
                //
                for (i = 0; i <= finest_level; i++)
                {
                    n_factor *= n_cycle[i];
                    dt_0 = std::min(dt_0, n_factor * dt_min[i]);
                }
            }
            else
            {
                if (verbose && ParallelDescriptor::IOProcessor())
                {
                   std::cout << "Nyx: Changing subcycling pattern. New pattern:\n";
                   for (i = 1; i <= finest_level; i++)
                    std::cout << "   Lev / n_cycle: " << i << " " << n_cycle[i] << '\n';
                }
            }
        }
    }
    else
    {
        dt_0 = fixed_dt;
    }

    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0.0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    // Shrink the time step if necessary in order to hit the next plot_z_value
    if (level == 0 && plot_z_values.size() > 0)
    {
        bool dt_changed = false;
        plot_z_est_time_step(dt_0,dt_changed);

        // Update the value of a if we didn't change dt in the call to plot_z_est_time_step.
        // If we didn't change dt there, then we have already done the integration.
        // If we did    change dt there, then we need to re-integrate here.
        if (dt_changed)
            integrate_comoving_a(cur_time,dt_0);
    }
    else
    {
        integrate_comoving_a(cur_time,dt_0);
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

void
Nyx::computeInitialDt (int                   finest_level,
                       int                   sub_cycle,
                       Array<int>&           n_cycle,
                       const Array<IntVect>& ref_ratio,
                       Array<Real>&          dt_level,
                       Real                  stop_time)
{
    //
    // Grids have been constructed, compute dt for all levels.
    //
    if (level > 0)
        return;

    int i;
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    if (parent->subcyclingMode() == "Optimal")
    {
        int new_cycle[finest_level+1];
        for (i = 0; i <= finest_level; i++)
            new_cycle[i] = n_cycle[i];
        Real dt_max[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            dt_max[i] = get_level(i).initial_time_step();
        }
        // Find the maximum number of cycles allowed
        int cycle_max[finest_level+1];
        cycle_max[0] = 1;
        for (i = 1; i <= finest_level; i++)
        {
            cycle_max[i] = parent->MaxRefRatio(i-1);
        }
        // estimate the amout of work to advance each level.
        Real est_work[finest_level+1];
        for (i = 0; i <= finest_level; i++)
        {
            est_work[i] = parent->getLevel(i).estimateWork();
        }
        dt_0 = parent->computeOptimalSubcycling(finest_level+1, new_cycle, dt_max, est_work, cycle_max);
        for (i = 0; i <= finest_level; i++)
        {
            n_cycle[i] = new_cycle[i];
        }
        if (verbose && ParallelDescriptor::IOProcessor() && finest_level > 0)
        {
           std::cout << "Nyx: Initial subcycling pattern:\n";
           for (i = 0; i <= finest_level; i++)
               std::cout << "Level " << i << ": " << n_cycle[i] << '\n';
        }
    }
    else
    {
        for (i = 0; i <= finest_level; i++)
        {
            dt_level[i] = get_level(i).initial_time_step();
            n_factor *= n_cycle[i];
            dt_0 = std::min(dt_0, n_factor * dt_level[i]);
        }
    }
    //
    // Limit dt's by the value of stop_time.
    //
    const Real eps = 0.001 * dt_0;
#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif
    if (stop_time >= 0)
    {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (i = 0; i <= finest_level; i++)
    {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }

    integrate_comoving_a(cur_time,dt_0);
}

bool
Nyx::writePlotNow ()
{
    if (level > 0)
        BoxLib::Error("Should only call writePlotNow at level 0!");

    bool found_one = false;

    if (plot_z_values.size() > 0)
    {

#ifdef NO_HYDRO
        Real prev_time = state[PhiGrav_Type].prevTime();
        Real  cur_time = state[PhiGrav_Type].curTime();
#else
        Real prev_time = state[State_Type].prevTime();
        Real  cur_time = state[State_Type].curTime();
#endif
        Real a_old = get_comoving_a(prev_time);
        Real z_old = (1. / a_old) - 1.;

        Real a_new = get_comoving_a( cur_time);
        Real z_new = (1. / a_new) - 1.;

        for (int i = 0; i < plot_z_values.size(); i++)
        {
            if (std::abs(z_new - plot_z_values[i]) < (0.01 * (z_old - z_new)) ) 
                found_one = true;
        }
    }

    if (found_one) {
        return true;
    } else {
        return false;
    }
}

void
Nyx::do_energy_diagnostics ()
{
    // nothing to see here, folks
}

void
Nyx::post_timestep (int iteration)
{
    //
    // Integration cycle on fine level grids is complete
    // do post_timestep stuff here.
    //
    int finest_level = parent->finestLevel();
    const int ncycle = parent->nCycle(level);
    
    //
    // Remove virtual particles at this level if we have any.
    //
    remove_virtual_particles();

    //
    // Remove Ghost particles on the final iteration 
    //
    if (iteration == ncycle)
        remove_ghost_particles();

    //
    // Redistribute if it is not the last subiteration
    //
    if (iteration < ncycle || level == 0)
    {
         for (int i = 0; i < theActiveParticles().size(); i++)
         {
             theActiveParticles()[i]->Redistribute(false, true, level, grav_n_grow);    
         }
    }

#ifndef NO_HYDRO
    if (do_reflux && level < finest_level)
    {
        MultiFab& S_new_crse = get_new_data(State_Type);
#ifdef GRAVITY
        MultiFab drho_and_drhoU;
#ifdef CGRAV	  
        if (do_grav &&
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav"))
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            // Define the update to rho and rhoU due to refluxing.
            drho_and_drhoU.define(grids, BL_SPACEDIM + 1, 0, Fab_allocate);
            MultiFab::Copy(drho_and_drhoU, S_new_crse, Density, 0,
                           BL_SPACEDIM + 1, 0);
            drho_and_drhoU.mult(-1.0);
        }
#endif // GRAVITY

        //We must reflux if the next finer level is subcycled relative to this level;
        //   otherwise the reflux was done as part of the multilevel advance
        if (parent->nCycle(level+1) != 1)
          reflux();

        // We need to do this before anything else because refluxing changes the
        // values of coarse cells underneath fine grids with the assumption
        // they'll be over-written by averaging down
        if (level < finest_level)
            average_down();

        // This needs to be done after any changes to the state from refluxing.
        enforce_nonnegative_species(S_new_crse);

#ifdef GRAVITY
#ifdef CGRAV
        if (do_grav && 
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
            && gravity->get_no_sync() == 0)
#else
        if (do_grav && gravity->get_gravity_type() == "PoissonGrav" && gravity->get_no_sync() == 0)
#endif
        {
            MultiFab::Add(drho_and_drhoU, S_new_crse, Density, 0, BL_SPACEDIM+1, 0);

            MultiFab dphi(grids, 1, 0);
            dphi.setVal(0);

            gravity->reflux_phi(level, dphi);

            // Compute (cross-level) gravity sync based on drho, dphi
            PArray<MultiFab> grad_delta_phi_cc(finest_level - level + 1,
                                               PArrayManage);
            for (int lev = level; lev <= finest_level; lev++)
            {
                grad_delta_phi_cc.set(lev - level,
                                      new MultiFab(get_level(lev).boxArray(),BL_SPACEDIM, 0));
                grad_delta_phi_cc[lev-level].setVal(0);
            }

            gravity->gravity_sync(level,finest_level,iteration,ncycle,drho_and_drhoU,dphi,grad_delta_phi_cc);
            dphi.clear();

            for (int lev = level; lev <= finest_level; lev++)
            {
                Real dt_lev = parent->dtLevel(lev);
                MultiFab&  S_new_lev = get_level(lev).get_new_data(State_Type);
                Real cur_time = state[State_Type].curTime();
                Real a_new = get_comoving_a(cur_time);

                const BoxArray& ba = get_level(lev).boxArray();
                MultiFab grad_phi_cc(ba, BL_SPACEDIM, 0);
                gravity->get_new_grav_vector(lev, grad_phi_cc, cur_time);

#ifdef _OPENMP
#pragma omp parallel	      
#endif
		{
		  FArrayBox sync_src;
		  FArrayBox dstate;

		  for (MFIter mfi(S_new_lev,true); mfi.isValid(); ++mfi)
                  {
                    const Box& bx = mfi.tilebox();
                    dstate.resize(bx, BL_SPACEDIM + 1);
                    if (lev == level)
                    {
		      dstate.copy(drho_and_drhoU[mfi]);
                    }
                    else
                    {
		      dstate.setVal(0);
                    }

                    // Compute sync source
                    sync_src.resize(bx, BL_SPACEDIM+1);
                    int i = mfi.index();
                    BL_FORT_PROC_CALL(FORT_SYNCGSRC,fort_syncgsrc)
                        (bx.loVect(), bx.hiVect(), BL_TO_FORTRAN(grad_phi_cc[i]),
                         BL_TO_FORTRAN(grad_delta_phi_cc[lev-level][i]),
                         BL_TO_FORTRAN(S_new_lev[i]), BL_TO_FORTRAN(dstate),
                         BL_TO_FORTRAN(sync_src), &a_new, dt_lev);

                    sync_src.mult(0.5 * dt_lev);
                    S_new_lev[mfi].plus(sync_src, 0, Xmom, BL_SPACEDIM);
                    S_new_lev[mfi].plus(sync_src, BL_SPACEDIM, Eden, 1);
		  }
		}
	    }
        }
#endif
    }
#endif // end ifndef NO_HYDRO 

    if (level < finest_level)
        average_down();

    if (level == 0)
    {
        int nstep = parent->levelSteps(0);
#ifndef NO_HYDRO
        if ( (do_hydro == 1) && (sum_interval > 0) && (nstep % sum_interval == 0))
        {
            sum_integrated_quantities();
        }
#endif
        write_info();

#if BL_USE_MPI
        // Memory monitoring:
        MemInfo* mInfo = MemInfo::GetInstance();
        char info[32];
        snprintf(info, sizeof(info), "Step %4d", nstep);
        mInfo->LogSummary(info);
#endif
    }

#ifndef NO_HYDRO
    if (do_hydro)
    {
       // Re-compute temperature after all the other updates.
       compute_new_temp();
    }
#endif
}

void
Nyx::post_restart ()
{
    if (level == 0)
        particle_post_restart(parent->theRestartFile());

    if (level == 0)
        comoving_a_post_restart(parent->theRestartFile());

#ifdef NO_HYDRO
    Real cur_time = state[PhiGrav_Type].curTime();
#else
    Real cur_time = state[State_Type].curTime();
#endif

    // Update the value of a only if restarting from chk00000
    //   (special case for which computeNewDt is *not* called from Amr::coarseTimeStep)
    if (level == 0 && cur_time == 0.0)
        integrate_comoving_a(cur_time,parent->dtLevel(0));

#ifdef TISF
     int blub = parent->finestLevel();
     BL_FORT_PROC_CALL(FORT_SET_FINEST_LEVEL, fort_set_finest_level)(&blub);
#endif

#ifdef GRAVITY

    if (do_grav)
    {
        if (level == 0)
        {
            for (int lev = 0; lev <= parent->finestLevel(); lev++)
            {
                AmrLevel& this_level = get_level(lev);
                gravity->install_level(lev, &this_level);
            }

            gravity->set_mass_offset(cur_time);

            if (
#ifdef CGRAV	  
            (gravity->get_gravity_type() == "PoissonGrav"||gravity->get_gravity_type() == "CompositeGrav")
#else
	    gravity->get_gravity_type() == "PoissonGrav"
#endif
)
            {
                // Do multilevel solve here.  We now store phi in the checkpoint file so we can use it
                //  at restart.
                int use_previous_phi_as_guess = 1;
                gravity->multilevel_solve_for_phi(0,parent->finestLevel(),use_previous_phi_as_guess);

#ifndef AGN
                if (do_dm_particles)
#endif
                {
                    for (int k = 0; k <= parent->finestLevel(); k++)
                    {
                        const BoxArray& ba = get_level(k).boxArray();
                        MultiFab grav_vec_new(ba, BL_SPACEDIM, 0);
                        gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
                    }
                }
            }
        }
    }
#endif

#ifndef NO_HYDRO
    if (level == 0)
    {
       // Need to compute this *before* regridding in case this is needed
       compute_average_density();
       set_small_values();
    }
#endif
}

#ifndef NO_HYDRO
void 
Nyx::set_small_values ()
{
       if (do_hydro == 0)
          return;

       Real small_pres;

       const Real cur_time = state[State_Type].curTime();
       Real a = get_comoving_a(cur_time);

       Real average_temperature;
       compute_average_temperature(average_temperature);
       //
       // Get the number of species from the network model.
       //
       BL_FORT_PROC_CALL(GET_NUM_SPEC, get_num_spec)(&NumSpec);
       BL_FORT_PROC_CALL(GET_NUM_AUX , get_num_aux )(&NumAux);

       BL_FORT_PROC_CALL(SET_SMALL_VALUES, set_small_values) 
            (&average_gas_density, &average_temperature,
             &a,  &small_dens, &small_temp, &small_pres);

       if (verbose && ParallelDescriptor::IOProcessor())
       {
          std::cout << "... setting small_dens to " << small_dens << '\n';
          std::cout << "... setting small_temp to " << small_temp << '\n';
          std::cout << "... setting small_pres to " << small_pres << '\n';
       }
}
#endif

void
Nyx::postCoarseTimeStep (Real cumtime)
{
#ifdef IN_SITU
   const Real cur_time = state[State_Type].curTime();
#if 1
   runInSituAnalysis(get_data(State_Type, cur_time), Geom(), nStep());
#else
   std::vector<Halo> halos = findHalos(*this, cur_time, State_Type);
   if (ParallelDescriptor::IOProcessor())
   {
       std::cout << "Found " << halos.size() << " halos." << std::endl;
       for (std::vector<Halo>::iterator it = halos.begin(); it != halos.end(); ++it)
           std::cout << it->pos << " " << it->mass << std::endl;
   }
#endif
#endif
    //
    // postCoarseTimeStep() is only called by level 0.
    //
    if (Nyx::theDMPC() && particle_move_type == "Random")
        particle_move_random();
}

void
Nyx::post_regrid (int lbase,
                  int new_finest)
{
#ifndef NO_HYDRO
#ifdef TISF
     BL_FORT_PROC_CALL(FORT_SET_FINEST_LEVEL, fort_set_finest_level)(&new_finest);
#endif
#endif

    if (level == lbase)
        particle_redistribute(lbase);

    int which_level_being_advanced = parent->level_being_advanced();
 
#ifdef GRAVITY
    bool do_grav_solve_here;
    if (which_level_being_advanced >= 0)
    {
        do_grav_solve_here = (level == which_level_being_advanced) && (lbase == which_level_being_advanced);
    } else {
        do_grav_solve_here = (level == lbase);
    }

    // Only do solve here if we will be using it in the timestep right after without re-solving,
    //      or if this is called from somewhere other than Amr::timeStep
    const Real cur_time = state[PhiGrav_Type].curTime();
    if (do_grav && (cur_time > 0) && do_grav_solve_here)
    {
#ifdef CGRAV	  
        if (gravity->get_gravity_type() == "PoissonGrav" || gravity->get_gravity_type() == "CompositeGrav")
#else
        if (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            int use_previous_phi_as_guess = 1;
            gravity->multilevel_solve_for_phi(level, new_finest, use_previous_phi_as_guess);
        }
    }
#endif
    delete fine_mask;
    fine_mask = 0;
}

void
Nyx::post_init (Real stop_time)
{
    if (level > 0)
        return;

    // If we restarted from a plotfile, we need to reset the level_steps counter
    if (!parent->theRestartPlotFile().empty())
        parent->setLevelSteps(0,nsteps_from_plotfile);

    //
    // Average data down from finer levels
    // so that conserved data is consistent between levels.
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--)
        get_level(k).average_down();

#ifdef GRAVITY
    if (do_grav)
    {
        const Real cur_time = state[PhiGrav_Type].curTime();
        if 
#ifdef CGRAV	  
            (gravity->get_gravity_type() == "PoissonGrav" || 
             gravity->get_gravity_type() == "CompositeGrav")
#else
	    (gravity->get_gravity_type() == "PoissonGrav")
#endif
        {
            //
            // Calculate offset before first multilevel solve.
            //
            gravity->set_mass_offset(cur_time);

            //
            // Solve on full multilevel hierarchy
            //
            gravity->multilevel_solve_for_phi(0, finest_level);
        }

        // Make this call just to fill the initial state data.
        for (int k = 0; k <= finest_level; k++)
        {
            const BoxArray& ba = get_level(k).boxArray();
            MultiFab grav_vec_new(ba, BL_SPACEDIM, 0);
            gravity->get_new_grav_vector(k, grav_vec_new, cur_time);
        }
    }
#endif

#ifndef NO_HYDRO
    if ( (do_hydro == 1) && (sum_interval > 0) && (parent->levelSteps(0) % sum_interval == 0) )
    {
        sum_integrated_quantities();
    }
    else
    {
        // Even if we don't call `sum_integrated_quantities` we need to compute
        // average_density before regridding
        compute_average_density();
    }

    if (do_hydro == 1)
    {
        set_small_values();
    }
#endif

    write_info();
}

int
Nyx::okToContinue ()
{
    if (level > 0)
        return 1;

    int test = 1;
    if (parent->dtLevel(0) < dt_cutoff)
        test = 0;

    if ((test == 1) && (final_a > 0))
    {
#ifdef NO_HYDRO
        Real cur_time = state[PhiGrav_Type].curTime();
#else
        Real cur_time = state[State_Type].curTime();
#endif
        Real a = get_comoving_a(cur_time);
        if (a >= final_a) test = 0;
        if (verbose && ParallelDescriptor::IOProcessor())
        {
            if (test == 0)
                std::cout << "...a " << a
                          << " is greater than or equal to final_a " << final_a
                          << '\n';
        }
    }
    return test;
}

#ifdef AUX_UPDATE
void
Nyx::advance_aux (Real time,
                  Real dt)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... special update for auxiliary variables \n";

    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_old,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        FArrayBox& old_fab = S_old[mfi];
        FArrayBox& new_fab = S_new[mfi];
        BL_FORT_PROC_CALL(FORT_AUXUPDATE, fort_auxupdate)
            (BL_TO_FORTRAN(old_fab), BL_TO_FORTRAN(new_fab), box.loVect(),
             box.hiVect(), &dt);
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::reflux ()
{
    BL_ASSERT(level<parent->finestLevel());

    const Real strt = ParallelDescriptor::second();

    get_flux_reg(level+1).Reflux(get_new_data(State_Type), 1.0, 0, 0, NUM_STATE,
                                 geom);

    if (show_timings)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        Real end = ParallelDescriptor::second() - strt;
        ParallelDescriptor::ReduceRealMax(end, IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "Nyx::reflux() at level "
                      << level
                      << " : time = "
                      << end << '\n';
    }
}
#endif // NO_HYDRO

void
Nyx::average_down ()
{
    if (level == parent->finestLevel()) return;

#ifndef NO_HYDRO
    // With State_Type we do DiagEOS_Type
    average_down(State_Type);
#endif

#ifdef GRAVITY
    average_down(PhiGrav_Type);
    average_down(Gravity_Type);
#endif
}

#ifndef NO_HYDRO
void
Nyx::enforce_nonnegative_species (MultiFab& S_new)
{
#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        BL_FORT_PROC_CALL(ENFORCE_NONNEGATIVE_SPECIES,
			  enforce_nonnegative_species)
	  (BL_TO_FORTRAN(S_new[mfi]), bx.loVect(), bx.hiVect(),
	   &print_fortran_warnings);
    }
}

void
Nyx::enforce_consistent_e (MultiFab& S)
{
#ifdef _OPENMP
#pragma omp parallel
#endif    
  for (MFIter mfi(S,true); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.tilebox();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();
        BL_FORT_PROC_CALL(ENFORCE_CONSISTENT_E, 
			  enforce_consistent_e)
	  (lo, hi, BL_TO_FORTRAN(S[mfi]));
    }
}
#endif

void
Nyx::average_down (int state_index)
{
#ifndef NO_HYDRO
    // We average DiagEOS_Type when average_down is called with State_Type
    if (state_index == DiagEOS_Type) return;
#endif

    if (level == parent->finestLevel()) return;

    Nyx& fine_level = get_level(level + 1);

    MultiFab& S_crse = get_new_data(state_index);
    MultiFab& S_fine = fine_level.get_new_data(state_index);
    const int num_comps = S_fine.nComp();

#ifndef NO_HYDRO
    if (state_index == State_Type)
    {
        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA(S_fine.boxArray().size());
 
        for (int i = 0; i < S_fine.boxArray().size(); ++i)
        {
            crse_S_fine_BA.set(i, BoxLib::coarsen(S_fine.boxArray()[i],
                                                  fine_ratio));
        }
        MultiFab crse_S_fine(crse_S_fine_BA, num_comps, 0);

        MultiFab& D_crse =            get_new_data(DiagEOS_Type);
        MultiFab& D_fine = fine_level.get_new_data(DiagEOS_Type);
        const int num_D_comps = D_fine.nComp();
        MultiFab crse_D_fine(crse_S_fine_BA, num_D_comps, 0);

#ifdef _OPENMP
#pragma omp parallel
#endif    
        for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
        {
 	    const Box&       overlap  = mfi.tilebox();

            const FArrayBox& fine_S_fab = S_fine[mfi];
            FArrayBox&       crse_S_fab = crse_S_fine[mfi];

            const FArrayBox& fine_D_fab = D_fine[mfi];
            FArrayBox&       crse_D_fab = crse_D_fine[mfi];

            BL_FORT_PROC_CALL(FORT_AVGDOWN, fort_avgdown)
                (BL_TO_FORTRAN(crse_S_fab), num_comps, 
		 BL_TO_FORTRAN(fine_S_fab),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());

            BL_FORT_PROC_CALL(FORT_AVGDOWN, fort_avgdown)
                (BL_TO_FORTRAN(crse_D_fab), num_D_comps, 
		 BL_TO_FORTRAN(fine_D_fab),
                 overlap.loVect(), overlap.hiVect(), fine_ratio.getVect());
        }
        D_crse.copy(crse_D_fine);
        S_crse.copy(crse_S_fine);
    }
    else
#endif
    {
      const Geometry& fine_geom = parent->Geom(level+1);
      const Geometry& crse_geom = parent->Geom(level  );
      BoxLib::average_down(S_fine,S_crse,fine_geom,crse_geom,0,num_comps,fine_ratio);
    }
}

void
Nyx::errorEst (TagBoxArray& tags,
               int          clearval,
               int          tagval,
               Real         time,
               int          n_error_buf,
               int          ngrow)
{
    const int*  domain_lo = geom.Domain().loVect();
    const int*  domain_hi = geom.Domain().hiVect();
    const Real* dx        = geom.CellSize();
    const Real* prob_lo   = geom.ProbLo();

    Real avg;

    for (int j = 0; j < err_list.size(); j++)
    {
        MultiFab* mf = derive(err_list[j].name(), time, err_list[j].nGrow());

        BL_ASSERT(!(mf == 0));

        if (err_list[j].errType() == ErrorRec::UseAverage)
        {
            if (err_list[j].name() == "density")
            {
                avg = average_gas_density;
            } 
            else if (err_list[j].name() == "particle_mass_density")
            {
                avg = average_dm_density;
#ifdef NEUTRINO_PARTICLES
                avg += average_neutr_density;
#endif
            }
            else if (err_list[j].name() == "total_density")
            {
                avg = average_total_density;
            }
#if 0
            else if (err_list[j].name() == "magvort")
            {
                avg = std::fabs(ave_lev_vorticity[level]);
                stddev = std_lev_vorticity[level];
                thresh = avg + std::max(stddev,avg);
                //std::cout << "errorEst, level " << level << ": magvort avg " << avg << ", stddev " << stddev
                //        << ", max " << std::max(stddev,avg) << ", thresh " << thresh << std::endl;
                thresh = std::max(thresh, 500.0);
                //std::cout << "errorEst, level " << level << ": thresh cut " << thresh << std::endl;
                avg = thresh;
            }
#endif
            else
            {
                if (ParallelDescriptor::IOProcessor())
                    std::cout << "Dont know the average of this variable "
                              << err_list[j].name() << '\n';
                avg = 0;
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            Array<int> itags;

            for (MFIter mfi(*mf,true); mfi.isValid(); ++mfi)
            {
		// FABs
		FArrayBox&  datfab  = (*mf)[mfi];
		TagBox&     tagfab  = tags[mfi];

		// Box in physical space
		int         idx     = mfi.index();
		RealBox     gridloc = RealBox(grids[idx],geom.CellSize(),geom.ProbLo());

		// tile box
		const Box&  tilebx  = mfi.tilebox();

		//fab box
		const Box&  datbox  = datfab.box();

		// We cannot pass tagfab to Fortran becuase it is BaseFab<char>.
		// So we are going to get a temporary integer array.
		tagfab.get_itags(itags, tilebx);

		// data pointer and index space
		int*        tptr    = itags.dataPtr();
		const int*  tlo     = tilebx.loVect();
		const int*  thi     = tilebx.hiVect();
		//
		const int*  lo      = tlo;
		const int*  hi      = thi;
		//
		const Real* xlo     = gridloc.lo();
		//
		Real*       dat     = datfab.dataPtr();
		const int*  dlo     = datbox.loVect();
		const int*  dhi     = datbox.hiVect();
		const int   ncomp   = datfab.nComp();
		
                if (err_list[j].errType() == ErrorRec::Standard)
                {
                    err_list[j].errFunc()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, xlo, prob_lo, &time, &level);
                }
                else if (err_list[j].errType() == ErrorRec::UseAverage)
                {
                   err_list[j].errFunc2()(tptr, ARLIM(tlo), ARLIM(thi), &tagval,
                                          &clearval, dat, ARLIM(dlo), ARLIM(dhi),
                                          lo, hi, &ncomp, domain_lo, domain_hi,
                                          dx, &level, &avg);
                }

                //
                // Don't forget to set the tags in the TagBox.
                //
                if (allow_untagging == 1)
                {
                    tagfab.tags_and_untags(itags, tilebx);
                }
                else
                {
                   tagfab.tags(itags, tilebx);
                }
            }
        }

        delete mf;
    }
}

MultiFab*
Nyx::derive (const std::string& name,
             Real               time,
             int                ngrow)
{
    if (name == "Rank")
    {
        MultiFab* derive_dat = new MultiFab(grids, 1, 0);
        for (MFIter mfi(*derive_dat); mfi.isValid(); ++mfi)
        {
           (*derive_dat)[mfi].setVal(ParallelDescriptor::MyProc());
        }
        return derive_dat;
    } else 
        return particle_derive(name, time, ngrow);
}

void
Nyx::derive (const std::string& name,
             Real               time,
             MultiFab&          mf,
             int                dcomp)
{
    AmrLevel::derive(name, time, mf, dcomp);
}

void
Nyx::network_init ()
{
    BL_FORT_PROC_CALL(NYX_NETWORK_INIT, nyx_network_init)();
}

#ifndef NO_HYDRO
void
Nyx::reset_internal_energy (MultiFab& S_new, MultiFab& D_new)
{
    // Synchronize (rho e) and (rho E) so they are consistent with each other

    const Real  cur_time = state[State_Type].curTime();
    Real        a        = get_comoving_a(cur_time);
    const Real* dx       = geom.CellSize();
    const Real  vol      = D_TERM(dx[0],*dx[1],*dx[2]);

    Real sum_energy_added = 0;
    Real sum_energy_total = 0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum_energy_added,sum_energy_total)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        Real s  = 0;
        Real se = 0;
        BL_FORT_PROC_CALL(RESET_INTERNAL_E, reset_internal_e)
            (BL_TO_FORTRAN(S_new[mfi]), BL_TO_FORTRAN(D_new[mfi]), 
             bx.loVect(), bx.hiVect(), 
             &print_fortran_warnings, &a, &s, &se);
        sum_energy_added += s;
        sum_energy_total += se;
    }

    if (verbose > 1)
    {
        Real sums[2] = {sum_energy_added,sum_energy_total};

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealSum(sums,2,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            sum_energy_added = vol*sums[0];
            sum_energy_total = vol*sums[1];

            if (sum_energy_added > (1.e-12)*sum_energy_total)
            {
                std::cout << "Adding to (rho E) "
                          << sum_energy_added
                          << " out of total (rho E) "
                          << sum_energy_total << '\n';
            }
        }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_new_temp ()
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real cur_time   = state[State_Type].curTime();

    reset_internal_energy(S_new,D_new);

    Real a = get_comoving_a(cur_time);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        BL_FORT_PROC_CALL(COMPUTE_TEMP, compute_temp)
            (bx.loVect(), bx.hiVect(), 
            BL_TO_FORTRAN(S_new[mfi]), 
            BL_TO_FORTRAN(D_new[mfi]), &a,
             &print_fortran_warnings);
    }

    // Compute the maximum temperature
    Real max_temp = D_new.norm0(Temp_comp);

    int imax = -1;
    int jmax = -1;
    int kmax = -1;
 
    Real den_maxt;

    // Find the cell which has the maximum temp -- but only if not the first
    // time step because in the first time step too many points have the same
    // value.
    Real prev_time   = state[State_Type].prevTime();
    if (prev_time > 0.0 && verbose > 0)
    {
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            BL_FORT_PROC_CALL(COMPUTE_MAX_TEMP_LOC, compute_max_temp_loc)
                (bx.loVect(), bx.hiVect(), 
                BL_TO_FORTRAN(S_new[mfi]), 
                BL_TO_FORTRAN(D_new[mfi]), 
                &max_temp,&den_maxt,&imax,&jmax,&kmax);
        }

        if (verbose > 1 && ParallelDescriptor::IOProcessor())
            if (imax > -1 && jmax > -1 && kmax > -1)
            {
              std::cout << "Maximum temp. at level " << level << " is " << max_temp 
                        << " at density " << den_maxt 
                        << " at (i,j,k) " << imax << " " << jmax << " " << kmax << std::endl;
            }
    }
}
#endif

#ifndef NO_HYDRO
void
Nyx::compute_rho_temp (Real& rho_T_avg, Real& T_avg, Real& T_meanrho)
{
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& D_new = get_new_data(DiagEOS_Type);

    Real rho_T_sum=0.0,   T_sum=0.0, T_meanrho_sum=0.0;
    Real   rho_sum=0.0, vol_sum=0.0,    vol_mn_sum=0.0;

#ifdef _OPENMP
#pragma omp parallel reduction(+:rho_T_sum, rho_sum, T_sum, T_meanrho_sum, vol_sum, vol_mn_sum)
#endif
    for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        BL_FORT_PROC_CALL(COMPUTE_RHO_TEMP, compute_rho_temp)
            (bx.loVect(), bx.hiVect(), geom.CellSize(), 
             BL_TO_FORTRAN(S_new[mfi]),
             BL_TO_FORTRAN(D_new[mfi]), &average_gas_density,  
             &rho_T_sum, &T_sum, &T_meanrho_sum, &rho_sum, &vol_sum, &vol_mn_sum);
    }
    Real sums[6] = {rho_T_sum, rho_sum, T_sum, T_meanrho_sum, vol_sum, vol_mn_sum};
    ParallelDescriptor::ReduceRealSum(sums,6);

    rho_T_avg = sums[0] / sums[1];  // density weighted T
        T_avg = sums[2] / sums[4];  // volume weighted T
    if (sums[5] > 0) { 
       T_meanrho = sums[3] / sums[5];  // T at mean density
       T_meanrho = pow(10.0, T_meanrho);
    }
}
#endif