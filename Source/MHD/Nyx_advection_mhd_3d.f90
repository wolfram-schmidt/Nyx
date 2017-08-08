
! :::
! ::: ----------------------------------------------------------------
! :::

      subroutine fort_advance_mhd(time,lo,hi,&
           uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
           uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
		   bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		   byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		   bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
		   bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
		   byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
		   bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
           ugdnvx,ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3, &
           ugdnvy,ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3, &
           ugdnvz,ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3, &
           src ,src_l1,src_l2,src_l3,src_h1,src_h2,src_h3, &
           grav,gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3, &
           delta,dt, &
           flux1,flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3, &
           flux2,flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3, &
           flux3,flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3, &
           courno,a_old,a_new,e_added,ke_added,print_fortran_warnings,do_grav) &
           bind(C, name="fort_advance_mhd")

!--------------------- Dependencies ------------------------------------------------
      use amrex_fort_module, only : rt => amrex_real
      use mempool_module, only : bl_allocate, bl_deallocate
	  use ct_upwind, only : corner_transport
	  use mhd_plm_module, only : plm
      use meth_params_module!, only : QVAR, NTHERM, NHYP, normalize_species, NVAR, URHO, UEDEN
      use enforce_module, only : enforce_nonnegative_species
      use bl_constants_module

      implicit none

!-------------------- Variables -----------------------------------------------------

      integer lo(3),hi(3),print_fortran_warnings,do_grav
      integer uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3
      integer uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
	  integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
	  integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
	  integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
	  integer bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
	  integer byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
	  integer bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
      integer ugdnvx_l1,ugdnvx_l2,ugdnvx_l3,ugdnvx_h1,ugdnvx_h2,ugdnvx_h3
      integer ugdnvy_l1,ugdnvy_l2,ugdnvy_l3,ugdnvy_h1,ugdnvy_h2,ugdnvy_h3
      integer ugdnvz_l1,ugdnvz_l2,ugdnvz_l3,ugdnvz_h1,ugdnvz_h2,ugdnvz_h3
      integer flux1_l1,flux1_l2,flux1_l3,flux1_h1,flux1_h2,flux1_h3
	  integer flux2_l1,flux2_l2,flux2_l3,flux2_h1,flux2_h2,flux2_h3
	  integer flux3_l1,flux3_l2,flux3_l3,flux3_h1,flux3_h2,flux3_h3
      integer src_l1,src_l2,src_l3,src_h1,src_h2,src_h3
      integer gv_l1,gv_l2,gv_l3,gv_h1,gv_h2,gv_h3
      real(rt)  uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3,  NVAR)
      real(rt)  uout(uout_l1:uout_h1, uout_l2:uout_h2, uout_l3:uout_h3, NVAR)
      real(rt)  bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt)  bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
      real(rt)  byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt)  byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
      real(rt)  bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
      real(rt)  bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)
      real(rt)  src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, NTHERM)
	  real(rt) ugdnvx(ugdnvx_l1:ugdnvx_h1,ugdnvx_l2:ugdnvx_h2,ugdnvx_l3:ugdnvx_h3)
      real(rt) ugdnvy(ugdnvy_l1:ugdnvy_h1,ugdnvy_l2:ugdnvy_h2,ugdnvy_l3:ugdnvy_h3)
      real(rt) ugdnvz(ugdnvz_l1:ugdnvz_h1,ugdnvz_l2:ugdnvz_h2,ugdnvz_l3:ugdnvz_h3)
      real(rt)  grav( gv_l1:gv_h1, gv_l2:gv_h2, gv_l3:gv_h3, 3)
      real(rt)  flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,NVAR)
	  real(rt)  flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,NVAR)
	  real(rt)  flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,NVAR)
      real(rt)  delta(3),dt,time,courno
      real(rt)  a_old, a_new
      real(rt)  e_added,ke_added

      ! Automatic arrays for workspace
      real(rt), pointer :: q(:,:,:,:)
      real(rt), pointer :: c(:,:,:)
      real(rt), pointer :: csml(:,:,:)
	  real(rt), pointer :: flatn(:,:,:)
 !     real(rt), pointer :: div(:,:,:)
 !     real(rt), pointer :: pdivu(:,:,:)
      real(rt), pointer :: srcQ(:,:,:,:)
	  real(rt), allocatable :: E(:,:,:,:,:)
	  real(rt), allocatable :: flx(:,:,:,:,:)
	  real(rt), allocatable :: qp(:,:,:,:,:)
	  real(rt), allocatable :: qm(:,:,:,:,:)

      real(rt) dx,dy,dz
      integer ngq,ngf
      integer q_l1, q_l2, q_l3, q_h1, q_h2, q_h3
      integer srcq_l1, srcq_l2, srcq_l3, srcq_h1, srcq_h2, srcq_h3
      integer 	:: i

      ngq = NHYP
      ngf = 1

      q_l1 = lo(1)-NHYP
      q_l2 = lo(2)-NHYP
      q_l3 = lo(3)-NHYP
      q_h1 = hi(1)+NHYP
      q_h2 = hi(2)+NHYP
      q_h3 = hi(3)+NHYP
      srcq_l1 = lo(1)-1
      srcq_l2 = lo(2)-1
      srcq_l3 = lo(3)-1
      srcq_h1 = hi(1)+1
      srcq_h2 = hi(2)+1
      srcq_h3 = hi(3)+1
	  write(*,*) "lo = ", lo
	  write(*,*) "hi = ", hi
	  write(*,*) "loq = ", q_l1, q_l2, q_l3
	  write(*,*) "hiq = ", q_h1, q_h2, q_h3

      call bl_allocate(     q, lo-NHYP, hi+NHYP, QVAR)
      call bl_allocate( flatn, lo-NHYP, hi+NHYP      )
      call bl_allocate(     c, lo-NHYP, hi+NHYP      )
      call bl_allocate(  csml, lo-NHYP, hi+NHYP      )
	  call bl_allocate(  srcQ, lo-1, hi+1, QVAR)
	  allocate(	E(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3, 3, 4))
	  allocate( flx(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR , 3))
	  allocate( qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR , 3))
	  allocate( qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR , 3))
	  q = 0.d0

      
      dx = delta(1)
      dy = delta(2)
      dz = delta(3)



!Step One, Calculate Primitives based on conservatives
	  call ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
					bxin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
					byin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
					bzin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
                    q , c , csml, flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                    src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                    srcQ, srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                    grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                    courno,dx,dy,dz,dt,ngq,ngf,a_old,a_new)

!Step Two, Interpolate Cell centered values to faces
	  call plm(q, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3,&	
				 bxin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
				 byin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
				 bzin, lo(1), lo(2), lo(3), hi(1), hi(2), hi(3), &
                 qp, qm, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, dx, dy, dz, dt ,a_old)
!write(*,*) "Bounds", q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
!qm(:,:,:,:,1) = q
!qm(:,:,:,:,2) = q
!qm(:,:,:,:,3) = q
!qp(:,:,:,:,1) = q
!qp(:,:,:,:,2) = q
!qp(:,:,:,:,3) = q
flx = 0.d0
!Step Three, Corner Couple and find the correct fluxes + electric fields
	  call corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &	
							flx, E, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, dx , dy, dz, dt)
!Step Four, Conservative update
      call consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  src ,src_l1 ,src_l2 ,src_l3 ,src_h1 ,src_h2 ,src_h3, &
                  flx , q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, &
                  lo ,hi ,dx ,dy ,dz ,dt ,a_old ,a_new)
!Step Five Magnetic Update
     call magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		 byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		 bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
		 bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
		 byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
		 bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
		 src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
		 E,q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3,lo, hi, dx, dy, dz, dt, a_old, a_new)

	  flux1(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,URHO:UEDEN) = flx(flux1_l1:flux1_h1,flux1_l2:flux1_h2, flux1_l3:flux1_h3,URHO:UEDEN,1)
	  flux2(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,URHO:UEDEN) = flx(flux2_l1:flux2_h1,flux2_l2:flux2_h2, flux2_l3:flux2_h3,URHO:UEDEN,2)
	  flux3(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,URHO:UEDEN) = flx(flux3_l1:flux3_h1,flux3_l2:flux3_h2, flux3_l3:flux3_h3,URHO:UEDEN,3)

      ! We are done with these here so can go ahead and free up the space.
      call bl_deallocate(q)
      call bl_deallocate(flatn)
      call bl_deallocate(c)
      call bl_deallocate(csml)
!      call bl_deallocate(div)
      call bl_deallocate(srcQ)
!     call bl_deallocate(pdivu)
	  deallocate(qm)
	  deallocate(qp)
	  deallocate(flx)
	  deallocate(E)


      ! Enforce the density >= small_dens.  Make sure we do this immediately after consup.
      call enforce_minimum_density(uin, uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3, &
                                        uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                                        lo,hi,print_fortran_warnings)
      
      if (do_grav .gt. 0)  then
          call add_grav_source(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                               uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                               grav, gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                               lo,hi,dx,dy,dz,dt,a_old,a_new,e_added,ke_added)
      endif
      ! Enforce species >= 0
 !     call enforce_nonnegative_species(uout,uout_l1,uout_l2,uout_l3, &
  !                                     uout_h1,uout_h2,uout_h3,lo,hi,0)

      ! Re-normalize the species
   !   if (normalize_species .eq. 1) then
    !     call normalize_new_species(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
     !                               lo,hi)
     ! end if

end subroutine fort_advance_mhd



! :::
! ::: ------------------------------------------------------------------
! :::

      subroutine ctoprim(lo,hi,uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3,&
						 bx, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		 				 by, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
						 bz, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
                         q,c,csml,flatn,  q_l1,  q_l2,  q_l3,  q_h1,  q_h2,  q_h3, &
                         src,  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3, &
                         srcQ,srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3, &
                         grav,gv_l1, gv_l2, gv_l3, gv_h1, gv_h2, gv_h3, &
                         courno,dx,dy,dz,dt,ngp,ngf,a_old,a_new)
      !
      !     Will give primitive variables on lo-ngp:hi+ngp, and flatn on lo-ngf:hi+ngf
      !     if use_flattening=1.  Declared dimensions of q,c,csml,flatn are given
      !     by DIMS(q).  This declared region is assumed to encompass lo-ngp:hi+ngp.
      !     Also, uflaten call assumes ngp>=ngf+3 (ie, primitve data is used by the
      !     routine that computes flatn).
      !
      use amrex_fort_module, only : rt => amrex_real
      use network, only : nspec, naux
      use eos_params_module
      use eos_module
!      use flatten_module
      use bl_constants_module
      use meth_params_module, only : NTHERM, URHO, UMX, UMY, UMZ, &
                                     UEDEN, UEINT, UFA, UFS, &
                                     QVAR, QRHO, QU, QV, QW, &
                                     QREINT, QPRES, QFA, QFS, &
									 QMAGX,  QMAGY, QMAGZ, &
                                     nadv, small_dens, small_pres, &
                                     gamma_const, gamma_minus_1, use_flattening

      implicit none

      real(rt), parameter:: small = 1.d-8

      integer lo(3), hi(3)
      integer  uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
	  integer bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
	  integer byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
	  integer bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
      integer    q_l1,   q_l2,   q_l3,   q_h1,   q_h2,   q_h3
      integer   gv_l1,  gv_l2,  gv_l3,  gv_h1,  gv_h2,  gv_h3
      integer  src_l1, src_l2, src_l3, src_h1, src_h2, src_h3
      integer srcq_l1,srcq_l2,srcq_l3,srcq_h1,srcq_h2,srcq_h3

      real(rt) :: uin(uin_l1:uin_h1,uin_l2:uin_h2,uin_l3:uin_h3,NTHERM)
	  real(rt) :: bx(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
      real(rt) :: by(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
      real(rt) :: bz(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
      real(rt) :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) !Contains Cell Centered Mag Field
      real(rt) :: c(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: csml(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: flatn(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3)
      real(rt) :: src( src_l1: src_h1, src_l2: src_h2, src_l3: src_h3,NTHERM)
      real(rt) :: srcQ(srcq_l1:srcq_h1,srcq_l2:srcq_h2,srcq_l3:srcq_h3,QVAR)
      real(rt) :: grav( gv_l1: gv_h1, gv_l2: gv_h2, gv_l3: gv_h3,3)
      real(rt) :: dx, dy, dz, dt, courno, a_old, a_new
      real(rt) :: dpdr, dpde

      integer          :: i, j, k
      integer          :: ngp, ngf, loq(3), hiq(3)
      integer          :: n, nq
      integer          :: iadv, ispec
	  integer		   :: btile(3)
      real(rt) :: courx, coury, courz, courmx, courmy, courmz
      real(rt) :: a_half, a_dot, rhoInv
      real(rt) :: dtdxaold, dtdyaold, dtdzaold, small_pres_over_dens

      do i=1,3
         loq(i) = lo(i)-ngp
         hiq(i) = hi(i)+ngp
      enddo
      !
      ! Make q (all but p), except put e in slot for rho.e, fix after eos call.
      ! The temperature is used as an initial guess for the eos call and will be overwritten.
      !
      do k = loq(3),hiq(3)
         do j = loq(2),hiq(2)
            do i = loq(1),hiq(1)
               if (uin(i,j,k,URHO) .le. ZERO) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                  print *,'>>> ... negative density ',uin(i,j,k,URHO)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
               end if

               rhoInv = ONE/uin(i,j,k,URHO)

               q(i,j,k,QRHO) = uin(i,j,k,URHO)
               q(i,j,k,QU)   = uin(i,j,k,UMX)*rhoInv
               q(i,j,k,QV)   = uin(i,j,k,UMY)*rhoInv
               q(i,j,k,QW)   = uin(i,j,k,UMZ)*rhoInv

               ! Convert "rho e" to "e"
               q(i,j,k,QREINT ) = uin(i,j,k,UEINT)*rhoInv

            enddo
         enddo
      enddo

      ! Load advected quatities, c, into q, assuming they arrived in uin as rho.c
      do iadv = 1, nadv
         n = UFA + iadv - 1
         nq = QFA + iadv - 1
         do k = loq(3),hiq(3)
            do j = loq(2),hiq(2)
               do i = loq(1),hiq(1)
                  q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
               enddo
            enddo
         enddo
      enddo

      ! Load chemical species and aux. variables, c, into q, assuming they arrived in uin as rho.c
!      if (UFS .gt. 0) then
!         do ispec = 1, nspec+naux
!            n  = UFS + ispec - 1
!            nq = QFS + ispec - 1
!            do k = loq(3),hiq(3)
!               do j = loq(2),hiq(2)
!                  do i = loq(1),hiq(1)
!                     q(i,j,k,nq) = uin(i,j,k,n)/q(i,j,k,QRHO)
!                  enddo
!               enddo
!            enddo
!         enddo
!      end if ! UFS > 0

      small_pres_over_dens = small_pres / small_dens
	  !Calculate Cell Centered Magnetic Field x
      do k = bxin_l3, bxin_h3 -1
         do j = bxin_l2, bxin_h2 -1
            do i = bxin_l1, bxin_h1 -1
					q(i,j,k,QMAGX) = 0.5d0*(bx(i+1,j,k) + bx(i,j,k))
				end do
			end do
		end do

      do k = byin_l3, byin_h3 -1
         do j = byin_l2, byin_h2 -1
            do i = byin_l1, byin_h1 -1
					q(i,j,k,QMAGY) = 0.5d0*(by(i,j+1,k) + by(i,j,k))
				end do
			end do
		end do
      do k = bzin_l3, bzin_h3 -1
         do j = bzin_l2, bzin_h2 -1
            do i = bzin_l1, bzin_h1 -1
					q(i,j,k,QMAGZ) = 0.5d0*(bz(i,j,k+1) + bz(i,j,k))
				end do
			end do
		end do
	   q(bxin_l1:bxin_h1,bxin_l2:bxin_h2, bxin_h3, QMAGX) = bx(:,:,bxin_h3)
	   q(bxin_l1:bxin_h1,bxin_h2, bxin_l3:bxin_h3, QMAGX) = bx(:,bxin_h2,:)	
	   q(bxin_h1,bxin_l2:bxin_h2, bxin_l3:bxin_h3, QMAGX) = bx(bxin_h1,:,:)
	   q(byin_l1:byin_h1,byin_l2:byin_h2, byin_h3, QMAGY) = by(:,:,byin_h3)
	   q(byin_l1:byin_h1,byin_h2, byin_l3:byin_h3, QMAGY) = by(:,byin_h2,:)	
	   q(byin_h1,byin_l2:byin_h2, byin_l3:byin_h3, QMAGY) = by(byin_h1,:,:)
	   q(bzin_l1:bzin_h1,bzin_l2:bzin_h2, bzin_h3, QMAGZ) = bz(:,:,bzin_h3)
	   q(bzin_l1:bzin_h1,bzin_h2, bzin_l3:bzin_h3, QMAGZ) = bz(:,bzin_h2,:)	
	   q(bzin_h1,bzin_l2:bzin_h2, bzin_l3:bzin_h3, QMAGZ) = bz(bzin_h1,:,:)
      ! Get p, T, c, csml using q state
      do k = loq(3), hiq(3)
         do j = loq(2), hiq(2)
            do i = loq(1), hiq(1)

               ! If necessary, reset the energy using small_temp
               if (q(i,j,k,QREINT) .lt. ZERO) then

!                 HACK HACK HACK 
!                 call nyx_eos_given_RT(q(i,j,k,QREINT),q(i,j,k,QPRES),q(i,j,k,QRHO), &
!                                       small_temp,diag_eos(i,j,k,NE_COMP),a_old)

                  if (q(i,j,k,QREINT) .lt. ZERO) then
                     !
                     ! A critical region since we usually can't write from threads.
                     !
                     print *,'   '
                     print *,'>>> Error: Nyx_advection_3d::ctoprim ',i,j,k
                     print *,'>>> ... new e from eos_given_RT call is negative ',q(i,j,k,QREINT)
                     print *,'    '
                     call bl_error("Error:: Nyx_advection_3d.f90 :: ctoprim")
                  end if
               end if

               ! Define the soundspeed from the EOS
               call nyx_eos_soundspeed(c(i,j,k), q(i,j,k,QRHO), q(i,j,k,QREINT), &
					   q(i,j,k,QMAGX), q(i,j,k,QMAGY), q(i,j,k,QMAGZ))

               ! Set csmal based on small_pres and small_dens
               csml(i,j,k) = sqrt(gamma_const * small_pres_over_dens)

               ! Convert "e" back to "rho e"
               q(i,j,k,QREINT) = q(i,j,k,QREINT)*q(i,j,k,QRHO)

               ! Pressure = (gamma - 1) * rho * e + 0.5 B dot B
               q(i,j,k,QPRES) = gamma_minus_1 * q(i,j,k,QREINT) &
				+ 0.5d0*(q(i,j,k,QMAGX)**2 + q(i,j,k,QMAGY)**2 + q(i,j,k,QMAGZ)**2)
            end do
         end do
      end do

      a_half = HALF * (a_old + a_new)
      a_dot   = (a_new - a_old) / dt

      ! Make sure these are initialized to zero.
      srcQ = ZERO

      ! NOTE - WE ASSUME HERE THAT src(i,j,k,URHO) = 0. --
      !        IF NOT THEN THE FORMULAE BELOW ARE INCOMPLETE.

      ! compute srcQ terms
  !    do k = lo(3)-1, hi(3)+1
  !       do j = lo(2)-1, hi(2)+1
  !          do i = lo(1)-1, hi(1)+1

  !             rhoInv = ONE/q(i,j,k,QRHO)

  !            srcQ(i,j,k,QRHO  ) = src(i,j,k,URHO)
  !             srcQ(i,j,k,QU    ) = src(i,j,k,UMX) * rhoInv - a_dot * q(i,j,k,QU) + &
  !                                  grav(i,j,k,1)
  !             srcQ(i,j,k,QV    ) = src(i,j,k,UMY) * rhoInv - a_dot * q(i,j,k,QV) + &
  !                                  grav(i,j,k,2)
  !             srcQ(i,j,k,QW    ) = src(i,j,k,UMZ) * rhoInv - a_dot * q(i,j,k,QW) + &
  !                                  grav(i,j,k,3)
  !             srcQ(i,j,k,QREINT) = src(i,j,k,UEDEN) - q(i,j,k,QU)*src(i,j,k,UMX) - &
  !                                                     q(i,j,k,QV)*src(i,j,k,UMY) - &
  !                                                     q(i,j,k,QW)*src(i,j,k,UMZ) - &
  !                                                     a_dot * THREE * gamma_minus_1 * q(i,j,k,QREINT)

   !            dpde = gamma_minus_1 * q(i,j,k,QRHO)
   !            dpdr = gamma_minus_1 * q(i,j,k,QREINT)/q(i,j,k,QRHO)
   !            srcQ(i,j,k,QPRES ) = dpde * srcQ(i,j,k,QREINT) * rhoInv &
   !                               + dpdr * srcQ(i,j,k,QRHO)

   !            if (UFS .gt. 0) then
   !               do ispec = 1,nspec+naux
   !                  srcQ(i,j,k,QFS+ispec-1) = src(i,j,k,UFS+ispec-1)*rhoInv
   !               enddo
   !            end if ! UFS > 0

   !            do iadv = 1,nadv
   !               srcQ(i,j,k,QFA+iadv-1) = src(i,j,k,UFA+iadv-1)*rhoInv
   !            enddo

   !         enddo
   !      enddo
   !   enddo

      ! Compute running max of Courant number over grids
      courmx = courno
      courmy = courno
      courmz = courno

      dtdxaold = dt / dx / a_old
      dtdyaold = dt / dy / a_old
      dtdzaold = dt / dz / a_old

      do k = lo(3),hi(3)
         do j = lo(2),hi(2)
            do i = lo(1),hi(1)

               courx = ( c(i,j,k)+abs(q(i,j,k,QU)) ) * dtdxaold
               coury = ( c(i,j,k)+abs(q(i,j,k,QV)) ) * dtdyaold
               courz = ( c(i,j,k)+abs(q(i,j,k,QW)) ) * dtdzaold

               courmx = max( courmx, courx )
               courmy = max( courmy, coury )
               courmz = max( courmz, courz )

               if (courx .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (u+c) * a * dt / dx > 1 ', courx
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... u, c                ',q(i,j,k,QU), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in x-dir in ctoprim")
               end if

               if (coury .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (v+c) * a * dt / dx > 1 ', coury
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... v, c                ',q(i,j,k,QV), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in y-dir in ctoprim")
               end if

               if (courz .gt. ONE) then
                  !
                  ! A critical region since we usually can't write from threads.
                  !
                  print *,'   '
                  print *,'>>> ... (w+c) * a * dt / dx > 1 ', courz
                  print *,'>>> ... at cell (i,j,k)   : ',i,j,k
                  print *,'>>> ... w, c                ',q(i,j,k,QW), c(i,j,k)
                  print *,'>>> ... density             ',q(i,j,k,QRHO)
                  call bl_error("Error:: Nyx_advection_3d.f90 :: CFL violation in z-dir in ctoprim")
               end if

            enddo
         enddo
      enddo
      courno = max( courmx, courmy, courmz )
      end subroutine ctoprim
! :::
! ::: ========================== Conservative Update ===============================================================
! ::: 

	subroutine consup(uin,uin_l1,uin_l2,uin_l3,uin_h1,uin_h2,uin_h3, &
                  uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, &
                  src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
                  flux,flux_l1,flux_l2,flux_l3,flux_h1,flux_h2,flux_h3, &
                  lo,hi,dx,dy,dz,dt,a_old,a_new)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : QVAR, QRHO, UMX,UMY,UMZ, QPRES, NVAR, URHO, UEDEN, UEINT

	implicit none

 	  integer,  intent(in)  :: uin_l1, uin_l2, uin_l3, uin_h1, uin_h2, uin_h3
	  integer,  intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
      integer,  intent(in)  :: flux_l1,flux_l2,flux_l3,flux_h1,flux_h2,flux_h3
	  integer,  intent(in)  :: src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
	  integer, intent(in) 	:: lo(3), hi(3)

	  real(rt), intent(in)  :: uin(uin_l1:uin_h1, uin_l2:uin_h2, uin_l3:uin_h3, NVAR)
	  real(rt), intent(in)  :: src(src_l1:src_h1,src_l2:src_h2,src_l3:src_h3, NVAR)
	  real(rt), intent(in)  :: flux(flux_l1:flux_h1,flux_l2:flux_h2,flux_l3:flux_h3,QVAR,3)
	  real(rt), intent(in) 	:: dx,dy,dz,dt,a_old, a_new 
	  real(rt), intent(out) :: uout(uout_l1:uout_h1,uout_l2:uout_h2, uout_l3:uout_h3,NVAR)


	  integer 				:: i, j, k		
	 uout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
	  !****TO DO ******* SOURCES
		do k = lo(3), hi(3)
			do j = lo(2), hi(2)
				do i = lo(1), hi(1)
					uout(i,j,k,URHO:UEDEN) = uin(i,j,k,URHO:UEDEN) - dt/dx*(flux(i+1,j,k,URHO:UEDEN,1) - flux(i,j,k,URHO:UEDEN,1)) &
											 -dt/dy*(flux(i,j+1,k,URHO:UEDEN,2) - flux(i,j,k,URHO:UEDEN,2)) &
											 -dt/dz*(flux(i,j,k+1,URHO:UEDEN,3) - flux(i,j,k,URHO:UEDEN,3)) !Add source terms later
					uout(i,j,k,UEINT) = uout(i,j,k,UEDEN) - 0.5*uout(i,j,k,URHO)*dot_product(uout(i,j,k,UMX:UMZ)/uout(i,j,k,URHO),uout(i,j,k,UMX:UMZ)/uout(i,j,k,URHO))
					if(uout(i,j,k,UEDEN).lt.0.d0) then
						write(*,*) "Negative Energy at", i, j, k, ' NRG in = ', uin(i,j,k,UEDEN)
						write(*,*) "flux x= ", - dt/dx*(flux(i+1,j,k,UEDEN,1) - flux(i,j,k,UEDEN,1))
						write(*,*) "flux y= ", - dt/dy*(flux(i,j+1,k,UEDEN,2) - flux(i,j,k,UEDEN,2))
						write(*,*) "flux z= ", - dt/dz*(flux(i,j,k+1,UEDEN,3) - flux(i,j,k,UEDEN,3))
						write(*,*) "internal NRG = ", uout(i,j,k,UEINT)
						pause
					endif

					!write(*,*) "Density = ",i,j,k, uout(i,j,k,URHO)
				enddo
			enddo
		enddo

	end subroutine consup

! :::
! ::: ========================== Magnetic Update ===============================================================
! ::: 

	subroutine magup(bxin, bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3, &
		 byin, byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3, &
		 bzin, bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3, &
		 bxout, bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3, &
		 byout, byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3, &
		 bzout, bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3, &
		 src ,  src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3, &
		 E,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,lo, hi, dx, dy, dz, dt, a_old, a_new)

     use amrex_fort_module, only : rt => amrex_real
     use meth_params_module, only : QVAR

	implicit none
	
		integer, intent(in)   :: bxin_l1, bxin_l2, bxin_l3, bxin_h1, bxin_h2, bxin_h3
		integer, intent(in)   :: byin_l1, byin_l2, byin_l3, byin_h1, byin_h2, byin_h3
		integer, intent(in)   :: bzin_l1, bzin_l2, bzin_l3, bzin_h1, bzin_h2, bzin_h3
		integer, intent(in)   :: bxout_l1, bxout_l2, bxout_l3, bxout_h1, bxout_h2, bxout_h3
		integer, intent(in)   :: byout_l1, byout_l2, byout_l3, byout_h1, byout_h2, byout_h3
		integer, intent(in)   :: bzout_l1, bzout_l2, bzout_l3, bzout_h1, bzout_h2, bzout_h3
		integer, intent(in)   :: src_l1,  src_l2,  src_l3,  src_h1,  src_h2,  src_h3
		integer, intent(in)   :: q_l1, q_l2,q_l3 ,q_h1,q_h2, q_h3
		integer, intent(in)   :: lo(3), hi(3)

		real(rt), intent(in)  :: bxin(bxin_l1:bxin_h1, bxin_l2:bxin_h2, bxin_l3:bxin_h3)
		real(rt), intent(in)  :: byin(byin_l1:byin_h1, byin_l2:byin_h2, byin_l3:byin_h3)
		real(rt), intent(in)  :: bzin(bzin_l1:bzin_h1, bzin_l2:bzin_h2, bzin_l3:bzin_h3)
		real(rt), intent(in)  :: src(src_l1:src_h1, src_l2:src_h2, src_l3:src_h3, QVAR)
		real(rt), intent(in)  :: E(q_l1:q_h1, q_l2:q_h2, q_l3:q_h3, 3, 4) 
		real(rt), intent(in)  :: dx, dy, dz, dt, a_old, a_new

		real(rt), intent(out) :: bxout(bxout_l1:bxout_h1, bxout_l2:bxout_h2, bxout_l3:bxout_h3)
		real(rt), intent(out) :: byout(byout_l1:byout_h1, byout_l2:byout_h2, byout_l3:byout_h3)
		real(rt), intent(out) :: bzout(bzout_l1:bzout_h1, bzout_l2:bzout_h2, bzout_l3:bzout_h3)

		integer				  :: i, j, k
		
		
		!***** TO DO ***** SOURCES
		bxout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
		byout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
		bzout(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
		!-------------------------------- bx --------------------------------------------------
			do k = lo(3), hi(3)
				do j = lo(2), hi(2)
					do i = lo(1), hi(1)
						bxout(i,j,k) = bxin(i,j,k) - dt/dx*(E(i,j,k,2,3) - E(i,j,k,2,1) - (E(i,j,k,3,2) - E(i,j,k,3,1)))
						print *, "Bx = ", bxout(i,j,k)
					enddo
				enddo
			enddo

		!------------------------------- by --------------------------------------------------
			do k = lo(3), hi(3)
				do j = lo(2), hi(2)
					do i = lo(1), hi(1)
						byout(i,j,k) = byin(i,j,k) - dt/dy*(E(i,j,k,3,2) - E(i,j,k,3,3) - (E(i,j,k,1,4) - E(i,j,k,1,1)))
						print *, "By = ", byout(i,j,k)
					enddo
				enddo
			enddo

		!------------------------------- bz --------------------------------------------------
			do k = lo(3), hi(3)
				do j = lo(2), hi(2)
					do i = lo(1), hi(1)
						bzout(i,j,k) = bzin(i,j,k) - dt/dz*(E(i,j,k,1,4) - E(i,j,k,1,3) - (E(i,j,k,2,3) - E(i,j,k,2,4)))
						print *, "Bz = ", bzout(i,j,k)
					enddo
				enddo
			enddo
		!--------------------------------------------------------------------------------------
	end subroutine magup
