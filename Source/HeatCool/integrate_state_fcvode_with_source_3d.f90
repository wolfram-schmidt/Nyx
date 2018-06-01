subroutine integrate_state_with_source_fcvode(lo, hi, &
                                state   , s_l1, s_l2, s_l3, s_h1, s_h2, s_h3, &
                                state_n ,sn_l1,sn_l2,sn_l3,sn_h1,sn_h2,sn_h3, &
                                diag_eos, d_l1, d_l2, d_l3, d_h1, d_h2, d_h3, &
                                hydro_src, src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                                reset_src,srcr_l1,srcr_l2,srcr_l3,srcr_h1,srcr_h2,srcr_h3, &
                                I_R, ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3, &
                                a, delta_time, min_iter, max_iter) &
                                bind(C, name="integrate_state_fcvode_with_source")
!
!   Calculates the sources to be added later on.
!
!   Parameters
!   ----------
!   lo : double array (3)
!       The low corner of the current box.
!   hi : double array (3)
!       The high corner of the current box.
!   state_* : double arrays
!       The state vars
!   diag_eos_* : double arrays
!       Temp and Ne
!   hydro_src_* : doubles arrays
!       The source terms to be added to state (iterative approx.)
!   reset_src_* : doubles arrays
!       The source terms based on the reset correction
!   double array (3)
!       The low corner of the entire domain
!   a : double
!       The current a
!   delta_time : double
!       time step size, in Mpc km^-1 s ~ 10^12 yr.
!
!   Returns
!   -------
!   state : double array (dims) @todo
!       The state vars
!
    use amrex_fort_module, only : rt => amrex_real
    use amrex_error_module, only : amrex_abort
    use meth_params_module, only : NVAR, URHO, UEDEN, UEINT, &
                                   NDIAG, TEMP_COMP, NE_COMP, ZHI_COMP, gamma_minus_1
    use bl_constants_module, only: M_PI, ONE, HALF
    use eos_params_module
    use network
    use eos_module, only: nyx_eos_T_given_Re, nyx_eos_given_RT
    use fundamental_constants_module
    use comoving_module, only: comoving_h, comoving_OmB
    use comoving_nd_module, only: fort_integrate_comoving_a
    use atomic_rates_module, only: YHELIUM
    use vode_aux_module    , only: JH_vode, JHe_vode, z_vode, i_vode, j_vode, k_vode
    use reion_aux_module   , only: zhi_flash, zheii_flash, flash_h, flash_he, &
                                   T_zhi, T_zheii, inhomogeneous_on
    use cvode_interface
    use fnvector_serial
    use fcvode_extras_src
    use, intrinsic :: iso_c_binding

    implicit none

    integer         , intent(in) :: lo(3), hi(3)
    integer         , intent(in) :: s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer         , intent(in) :: sn_l1, sn_l2, sn_l3, sn_h1, sn_h2, sn_h3
    integer         , intent(in) :: d_l1, d_l2, d_l3, d_h1, d_h2, d_h3
    integer         , intent(in) :: src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
    integer         , intent(in) :: srcr_l1, srcr_l2, srcr_l3, srcr_h1, srcr_h2, srcr_h3
    integer         , intent(in) :: ir_l1, ir_l2, ir_l3, ir_h1, ir_h2, ir_h3
    real(rt), intent(in   ) ::    state(s_l1:s_h1, s_l2:s_h2,s_l3:s_h3, NVAR)
    real(rt), intent(inout) ::  state_n(sn_l1:sn_h1, sn_l2:sn_h2,sn_l3:sn_h3, NVAR)
    real(rt), intent(inout) :: diag_eos(d_l1:d_h1, d_l2:d_h2,d_l3:d_h3, NDIAG)
    real(rt), intent(in   ) :: hydro_src(src_l1:src_h1, src_l2:src_h2,src_l3:src_h3, NVAR)
    real(rt), intent(in   ) :: reset_src(srcr_l1:srcr_h1, srcr_l2:srcr_h2,srcr_l3:srcr_h3, 1)
    real(rt), intent(inout) :: I_R(ir_l1:ir_h1, ir_l2:ir_h2,ir_l3:ir_h3)
    real(rt), intent(in)    :: a, delta_time
    integer         , intent(inout) :: max_iter, min_iter

    integer :: i, j, k
    real(rt) :: asq,aendsq,ahalf,ahalf_inv,delta_rho,delta_e,delta_rhoe
    real(rt) :: z, z_end, a_end, rho, H_reion_z, He_reion_z
    real(rt) :: rho_orig, T_orig, ne_orig, e_orig
    real(rt) :: rho_out, T_out, ne_out, e_out
    real(rt) :: rho_src, rhoe_src, e_src
    real(rt) :: mu, mean_rhob, T_H, T_He
    real(rt) :: species(5)

    integer(c_int) :: ierr       ! error flag from C functions
    real(c_double) :: tstart     ! initial time

    real(c_double), pointer, dimension(:) :: atol
    real(c_double) :: rtol
    type(c_ptr) :: sunvec_y      ! sundials vector
    type(c_ptr) :: CVmem         ! CVODE memory
    integer(c_long), parameter :: neq = 2
    real(c_double), pointer :: yvec(:)
    type(c_ptr) :: sunvec_atol

    allocate(yvec(neq))
    allocate(atol(neq))

    z = 1.d0/a - 1.d0
    call fort_integrate_comoving_a(a, a_end, delta_time)
    z_end = 1.0d0/a_end - 1.0d0
 
    asq = a*a
    aendsq = a_end*a_end
    ahalf     = HALF * (a + a_end)
    ahalf_inv  = ONE / ahalf

    mean_rhob = comoving_OmB * 3.d0*(comoving_h*100.d0)**2 / (8.d0*M_PI*Gconst)

    ! Flash reionization?
    if ((flash_h .eqv. .true.) .and. (z .gt. zhi_flash)) then
       JH_vode = 0
    else
       JH_vode = 1
    endif
    if ((flash_he .eqv. .true.) .and. (z .gt. zheii_flash)) then
       JHe_vode = 0
    else
       JHe_vode = 1
    endif

    if (flash_h ) H_reion_z  = zhi_flash
    if (flash_he) He_reion_z = zheii_flash

    ! Note that (lo,hi) define the region of the box containing the grow cells
    ! Do *not* assume this is just the valid region
    ! apply heating-cooling to UEDEN and UEINT
    sunvec_y = N_VMake_Serial(NEQ, yvec)
    if (.not. c_associated(sunvec_y)) then
        call amrex_abort('integrate_state_fcvode: sunvec = NULL')
    end if

    CVmem = FCVodeCreate(CV_BDF, CV_NEWTON)
    if (.not. c_associated(CVmem)) then
        call amrex_abort('integrate_state_fcvode: CVmem = NULL')
    end if

    tstart = 0.0
    ! CVodeMalloc allocates variables and initialize the solver. We can initialize the solver with junk because once we enter the
    ! (i,j,k) loop we will immediately call fcvreinit which reuses the same memory allocated from CVodeMalloc but sets up new
    ! initial conditions.
    ierr = FCVodeInit(CVmem, c_funloc(RhsFn_src), tstart, sunvec_y)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVodeInit() failed')
    end if

    ! Set dummy tolerances. These will be overwritten as soon as we enter the loop and reinitialize the solver.
    rtol = 1.0d-5
    atol(1) = 1.0d-10
    atol(2) = 1.0d-10
    sunvec_atol = N_VMake_Serial(neq, atol)

    ierr = FCVodeSVtolerances(CVmem, rtol, sunvec_atol)
    if (ierr /= 0) then
      call amrex_abort('integrate_state_fcvode: FCVodeSStolerances() failed')
    end if

    ierr = FCVDense(CVmem, neq)
    if (ierr /= 0) then
       call amrex_abort('integrate_state_fcvode: FCVDense() failed')
    end if

    do k = lo(3),hi(3)
        do j = lo(2),hi(2)
            do i = lo(1),hi(1)

                ! Original values
                rho_orig  = state(i,j,k,URHO)
                e_orig    = state(i,j,k,UEINT) / rho_orig
                T_orig    = diag_eos(i,j,k,TEMP_COMP)
                ne_orig   = diag_eos(i,j,k,  NE_COMP)

                if (inhomogeneous_on) then
                   H_reion_z = diag_eos(i,j,k,ZHI_COMP)
                   if (z .gt. H_reion_z) then
                      JH_vode = 0
                   else
                      JH_vode = 1
                   endif
                endif

                if (e_orig .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e entering strang integration ', z, i,j,k, rho_orig/mean_rhob, e_orig
                    call bl_abort('bad e in strang')
                    !$OMP END CRITICAL
                end if

                rho_src  = hydro_src(i,j,k,URHO)
                rhoe_src = hydro_src(i,j,k,UEINT)

                !                             
                ! This term satisfies the equation anewsq * (rho_new e_new) = 
                !                                  aoldsq * (rho_old e_old) + dt * H_{rho e} + anewsq * (reset_src)
                !  where e_new = e_old + dt * e_src                           
                e_src = ( ((asq*state(i,j,k,UEINT) + delta_time * rhoe_src ) / aendsq )/ &
                        (      state(i,j,k,URHO ) + delta_time * rho_src  ) - e_orig) / delta_time
                e_src = ( ((asq*state(i,j,k,UEINT) + delta_time * rhoe_src ) / aendsq + reset_src(i,j,k,1))/ &
                        (      state(i,j,k,URHO ) + delta_time * rho_src  ) - e_orig) / delta_time

                i_vode = i
                j_vode = j
                k_vode = k

                call fcvode_wrapper_with_source(delta_time,rho_orig,T_orig,ne_orig,e_orig, neq, CVmem, &
                                                        sunvec_y, yvec, rho_out ,T_out ,ne_out ,e_out,rho_src,e_src)
                !                             
                ! I_R satisfies the equation anewsq * (rho_out  e_out ) = 
                !                            aoldsq * (rho_orig e_orig) + dt * a_half * I_R + dt * H_{rho e} + anewsq * reset_src
                I_R(i,j,k) = ( aendsq * rho_out *e_out - ( (asq*rho_orig* e_orig + delta_time*rhoe_src) ) ) / &
                              (delta_time * ahalf) - aendsq * reset_src(i,j,k,1) / (delta_time * ahalf)

                if (e_out .lt. 0.d0) then
                    !$OMP CRITICAL
                    print *,'negative e exiting strang integration ', z, i,j,k, rho_orig/mean_rhob, e_out, &
                             state_n(i,j,k,UEINT)/state_n(i,j,k,URHO)
                    call flush(6)
                    !$OMP END CRITICAL
                    T_out  = 10.0
                    ne_out = 0.0
                    mu     = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)
                    e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                    !call bl_abort('bad e out of strang')
                end if

                ! Update T and ne (do not use stuff computed in f_rhs, per vode manual)
                call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho_out, e_out, a, species)

                ! Instanteneous heating from reionization:
                T_H = 0.0d0
                if (inhomogeneous_on .or. flash_h) then
                   if ((H_reion_z  .lt. z) .and. (H_reion_z  .ge. z_end)) T_H  = (1.0d0 - species(2))*max((T_zhi-T_out), 0.0d0)
                endif

                T_He = 0.0d0
                if (flash_he) then
                   if ((He_reion_z .lt. z) .and. (He_reion_z .ge. z_end)) T_He = (1.0d0 - species(5))*max((T_zheii-T_out), 0.0d0)
                endif

                if ((T_H .gt. 0.0d0) .or. (T_He .gt. 0.0d0)) then
                   T_out = T_out + T_H + T_He                            ! For simplicity, we assume
                   ne_out = 1.0d0 + YHELIUM                              !    completely ionized medium at
                   if (T_He .gt. 0.0d0) ne_out = ne_out + YHELIUM        !    this point.  It's a very minor
                   mu = (1.0d0+4.0d0*YHELIUM) / (1.0d0+YHELIUM+ne_out)   !    detail compared to the overall approximation.
                   e_out  = T_out / (gamma_minus_1 * mp_over_kB * mu)
                   call nyx_eos_T_given_Re(JH_vode, JHe_vode, T_out, ne_out, rho, e_out, a, species)
                endif

                ! Update (rho e) and (rho E)
                ! Note that we add to state_n because those already have hydro_source in them
                state_n(i,j,k,UEINT) = state_n(i,j,k,UEINT) + delta_time * ahalf * I_R(i,j,k) / aendsq
                state_n(i,j,k,UEDEN) = state_n(i,j,k,UEDEN) + delta_time * ahalf * I_R(i,j,k) / aendsq
                 
                ! Update T and ne
                diag_eos(i,j,k,TEMP_COMP) = T_out
                diag_eos(i,j,k,  NE_COMP) = ne_out

            end do ! i
        end do ! j
    end do ! k

    call N_VDestroy_Serial(sunvec_y)
    call N_VDestroy_Serial(sunvec_atol)
    call FCVodeFree(cvmem)

    deallocate(yvec)
    deallocate(atol)

end subroutine integrate_state_with_source_fcvode
