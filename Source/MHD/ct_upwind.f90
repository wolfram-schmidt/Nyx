module ct_upwind

 use amrex_fort_module, only : rt => amrex_real
 use hlld_solver, only : hlld
 use meth_mhd_params_module

implicit none 

private interpolate, primtocons
public corner_transport

contains

subroutine corner_transport( !Cons
			     !left and right Prim
			     !Fluxes
				)

 use amrex_fort_module, only : rt => amrex_real
 use meth_mhd_params_module
implicit none

	integer, intent(in)   :: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3

	real(rt), intent(in)  :: qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3)
	real(rt), intent(in)  :: qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)	
	real(rt), intent(out) :: elec(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)
 
	real(rt)			  :: um(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3) !PtoC Vars
	real(rt)			  :: up(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3)
	real(rt)			  :: cons_temp_L(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3,2) !2D Temporary Conservative Vars
	real(rt)			  :: cons_temp_R(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3,2) !2D Temporary Conservative Vars
	real(rt)			  :: cons_half_L(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3) !Flux Corrected Conservative Vars
	real(rt)			  :: cons_half_R(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR,3)
	real(rt)			  :: flx1D(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3) !Flux1d for all directions
	real(rt)			  :: flx2D(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3, 2) !Flux2d for all directions 2 perpendicular directions
	real(rt) 			  :: Etemp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,3) !Temporary Electric Field

	
!Prim to Cons
	call PrimToCons(qm, um, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3)
	call PrimToCons(qp, up, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3)
	
!Calculate Flux 1D
	!x-dir
	call hlld(qm(:,:,:,:,1),qp(:,:,:,:,1),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx1D(:,:,:,:,1),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	!y-dir	
	call hlld(qm(:,:,:,:,2),qp(:,:,:,:,2),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx1D(:,:,:,:,2),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	!z-dir
	call hlld(qm(:,:,:,:,3),qp(:,:,:,:,3),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx1D(:,:,:,:,3),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)

!Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields

	call elec_interp()

!Corner Couple

	call corner_couple() !Correct Conservative vars using Transverse Fluxes
	call corner_couple_mag()

!Calculate Flux 2D
do i = 1,2
	!x-dir
	call hlld(cons_temp_L(:,:,:,:,1,i),cons_temp_R(:,:,:,:,1,i),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx2D(:,:,:,:,1,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	!y-dir	
	call hlld(cons_temp_L(:,:,:,:,2,i),cons_temp_R(:,:,:,:,2,i),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx2D(:,:,:,:,2,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	!z-dir
	call hlld(cons_temp_L(:,:,:,:,3,i),cons_temp_R(:,:,:,:,3,i),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx2D(:,:,:,:,3,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)
enddo
	
!Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse flx1D
	flx1D(:,:,:,:,:) = 0.5d0*(flx2D(:,:,:,:,:,1) + flx2D(:,:,:,:,:,2))
	
	call elec_interp()

!Half Step conservative vars
	call half_step()
	call half_step_mag()

!Final Fluxes
	!x-dir
	call hlld(cons_half_L(:,:,:,:,1),cons_half_R(:,:,:,:,1),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx(:,:,:,:,1),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	!y-dir	
	call hlld(cons_half_L(:,:,:,:,2),cons_half_R(:,:,:,:,2),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx(:,:,:,:,2),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	!z-dir
	call hlld(cons_temp_L(:,:,:,:,3),cons_temp_R(:,:,:,:,3),qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,flx(:,:,:,:,3),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)

!Primitive update

	call prim_half()

!Final Electric Field Update
	call elec_interp()

end subroutine corner_transport

subroutine PrimToCons(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_mhd_params_module

implicit none

	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(in)	::q(q_l1,q_l2,q_l3,q_h1,q_h2, q_h3)
	real(rt), intent(out)	::u(q_l1,q_l2,q_l3,q_h1,q_h2, q_h3)

	u(QRHO)  = q(QRHO)
	u(QU)    = q(QRHO)*q(QU)
	u(QV)    = q(QRHO)*q(QV)
	u(QW)    = q(QRHO)*q(QW)
	u(QPRES) = q(QPRES)*gamma_minus_1*q(QRHO)
	u(QMAGX:QMAGZ) = q(QMAGX:QMAGZ)

end subroutine PrimToCons


end module ct_updwind
