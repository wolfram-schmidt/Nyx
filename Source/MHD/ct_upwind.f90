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

	real(rt), intent(in)  :: q(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR) !prim vars at time t^n
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

	call elec_interp(Etemp, q, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
			flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)

!Corner Couple

	call corner_couple(cons_temp_L, cons_temp_R, um, up, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,&
					   flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3) !Correct Conservative vars using Transverse Fluxes
	call corner_couple_mag(cons_temp_L, cons_temp_R, um, up, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,&
						Etemp) !Correct Magnetic Vars using Etemp

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
	real(rt), intent(in)	::q(q_l1,q_l2,q_l3,q_h1,q_h2, q_h3,QVAR)
	real(rt), intent(out)	::u(q_l1,q_l2,q_l3,q_h1,q_h2, q_h3,QVAR)

 do k = q_l3,q_h3
 	 do j = q_l2,q_h2
 		 do i = q_l1, q_l2
 			u(i,j,k,QRHO)  = q(i,j,k,QRHO)
 			u(i,j,k,QU)    = q(i,j,k,QRHO)*q(i,j,k,QU)
 			u(i,j,k,QV)    = q(i,j,k,QRHO)*q(i,j,k,QV)
 			u(i,j,k,QW)    = q(i,j,k,QRHO)*q(i,j,k,QW)
			u(i,j,k,QPRES) = q(i,j,k,QPRES)*gamma_minus_1*q(i,j,k,QRHO) &
 							 + 0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ)*q(i,j,k,QMAGX:QMAGZ))) !Energy
			u(i,j,k,QMAGX:QMAGZ) = q(i,j,k,QMAGX:QMAGZ)
		 enddo
	 enddo
 enddo
end subroutine PrimToCons


subroutine corner_couple(uL, uR, um, up, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3,&
					   flx, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)
use amrex_fort_module, only : rt => amrex_real
use meth_mhd_params_module

implicit none
	
	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	
	real(rt), intent(in)	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)

	real(rt), intent(out)	::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
	real(rt), intent(out)	::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)

	integer					:: i ,j ,k

	do k = q_l3,q_h3
		do j = q_l2,q_h2
			do i = q_l1,q_h1
	!Left Corrected States
				uL(i,j,k,QRHO:QPRES,1,1) = um(i,j,k,QRHO:QPRES,1) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,2) - flx(i,j-1,k,QRHO:QPRES,2))
				uL(i,j,k,QRHO:QPRES,1,2) = um(i,j,k,QRHO:QPRES,1) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,3) - flx(i,j-1,k,QRHO:QPRES,3))
				uL(i,j,k,QRHO:QPRES,2,1) = um(i,j,k,QRHO:QPRES,2) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,1) - flx(i,j-1,k,QRHO:QPRES,1))
				uL(i,j,k,QRHO:QPRES,2,2) = um(i,j,k,QRHO:QPRES,2) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,3) - flx(i,j-1,k,QRHO:QPRES,3))
				uL(i,j,k,QRHO:QPRES,3,1) = um(i,j,k,QRHO:QPRES,3) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,1) - flx(i,j-1,k,QRHO:QPRES,1))
				uL(i,j,k,QRHO:QPRES,3,2) = um(i,j,k,QRHO:QPRES,3) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,2) - flx(i,j-1,k,QRHO:QPRES,2))
	!Right Corrected States
				uR(i,j,k,QRHO:QPRES,1,1) = up(i,j,k,QRHO:QPRES,1) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,2) - flx(i,j-1,k,QRHO:QPRES,2))
				uR(i,j,k,QRHO:QPRES,1,2) = up(i,j,k,QRHO:QPRES,1) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,3) - flx(i,j-1,k,QRHO:QPRES,3))
				uR(i,j,k,QRHO:QPRES,2,1) = up(i,j,k,QRHO:QPRES,2) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,1) - flx(i,j-1,k,QRHO:QPRES,1))
				uR(i,j,k,QRHO:QPRES,2,2) = up(i,j,k,QRHO:QPRES,2) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,3) - flx(i,j-1,k,QRHO:QPRES,3))
				uR(i,j,k,QRHO:QPRES,3,1) = up(i,j,k,QRHO:QPRES,3) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,1) - flx(i,j-1,k,QRHO:QPRES,1))
				uR(i,j,k,QRHO:QPRES,3,2) = up(i,j,k,QRHO:QPRES,3) - dt/(3.d0*dx)*(flx(i,j,k,QRHO:QPRES,2) - flx(i,j-1,k,QRHO:QPRES,2))
			enddo
		enddo
	enddo
end subroutine corner_couple

subroutine corner_couple_mag(uL, uR, um, up,&
                            q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
       						Etemp)
use amrex_fort_module, only : rt => amrex_real
use meth_mhd_params_module

!Correction using Faraday's Law
implicit none

	integer, intent(in)			::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(inout)		::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(inout)     ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)    	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)    	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)	
	real(rt), intent(in)		::E(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,3,2,2)
	
	integer					:: i ,j ,k

   	do k = q_l3,q_h3
		do j = q_l2,q_h2
			do i = q_l1,q_h1	
		!Left State
				!X-direction 
				!-> Affected by Y flux
				uL(i,j,k,QMAGX,1) = um(i,j,k,QMAGX,1) + deltat/(3.d0*deltax)*&
										(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
				uL(i,j,k,QMAGZ,1) = um(i,j,k,QMAGZ,1) + deltat/(6.d0*deltax)*&
										((E(i,j,k,1,2,2) - E(i,j,k,1,2,1)) + &
										(E(i,j,k,1,1,2) - E(i,j,k,1,1,1))
				!-> Affected by Z flux
				uL(i,j,k,QMAGX,1) = uL(i,j,k,QMAGX,1) - deltat/(3.d0*deltax)*&
									 	(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				uL(i,j,k,QMAGY,1) = um(i,j,k,QMAGY,1) - deltat/(6.d0*deltax)*&
									 	((E(i,j,k,1,2,2) - E(i,j,k,1,2,1)) + &
										(E(i,j,k,1,1,2) - E(i,j,k,1,1,1))
				!Y-direction
				!-> Affected by X flux
				uL(i,j,k,QMAGY,2) = um(i,j,k,QMAGY,2) - deltat/(3.d0*deltay)*&
										(E(i,j,k,3,2,2) - E(i,j,k,3,2,1))
				uL(i,j,k,QMAGZ,2) = um(i,j,k,QMAGZ,2) - deltat/(6.d0*deltay)*&
									((E(i,j,k,2,2,2) - E(i,j,k,2,2,1)) + &
									(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				!-> Affected by Z flux
				uL(i,j,k,QMAGY,2) = uL(i,j,k,QMAGY,2) + deltat/(3.d0*deltay)*&
										(E(i,j,k,1,2,2) - E(i,j,k,1,2,1))
				uL(i,j,k,QMAGX,2) = um(i,j,k,QMAGX,2) + deltat/(6.d0*deltay)*&
									((E(i,j,k,2,2,2) - E(i,j,k,2,2,1)) + &
									(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))

				!Z-Direction
				!-> Affected by X flux
				uL(i,j,k,QMAGZ,3) = um(i,j,k,QMAGZ,3) - deltat/(3.d0*deltaz)*&
										(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				uL(i,j,k,QMAGY,3) = um(i,j,k,QMAGZ,3) - deltat/(6.d0*deltaz)*&
									((E(i,j,k,3,2,2) - E(i,j,k,3,2,1)) + &
									(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
				!-> Affected by Y flux
				uL(i,j,k,QMAGZ,3) = uL(i,j,k,QMAGZ,3) + deltat/(3.d0*deltaz)*&
										(E(i,j,k,1,2,2) - E(i,j,k,1,2,1))
				UL(i,j,k,QMAGX,3) = um(i,j,k,QMAGX,3) + deltat/(6.d0*deltaz)*&
									((E(i,j,k,3,2,2) - E(i,j,k,3,2,1)) + &
									(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
		!Right State
				!X-direction 
				!-> Affected by Y flux
				uR(i,j,k,QMAGX,1) = up(i,j,k,QMAGX,1) + deltat/(3.d0*deltax)*&
										(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
				uR(i,j,k,QMAGZ,1) = up(i,j,k,QMAGZ,1) + deltat/(6.d0*deltax)*&
										((E(i,j,k,1,2,2) - E(i,j,k,1,2,1)) + &
										(E(i,j,k,1,1,2) - E(i,j,k,1,1,1))
				!-> Affected by Z flux
				uR(i,j,k,QMAGX,1) = uR(i,j,k,QMAGX,1) - deltat/(3.d0*deltax)*&
									 	(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				uR(i,j,k,QMAGY,1) = up(i,j,k,QMAGY,1) - deltat/(6.d0*deltax)*&
									 	((E(i,j,k,1,2,2) - E(i,j,k,1,2,1)) + &
										(E(i,j,k,1,1,2) - E(i,j,k,1,1,1))
				!Y-direction
				!-> Affected by X flux
				uR(i,j,k,QMAGY,2) = up(i,j,k,QMAGY,2) - deltat/(3.d0*deltay)*&
										(E(i,j,k,3,2,2) - E(i,j,k,3,2,1))
				uR(i,j,k,QMAGZ,2) = up(i,j,k,QMAGZ,2) - deltat/(6.d0*deltay)*&
									((E(i,j,k,2,2,2) - E(i,j,k,2,2,1)) + &
									(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				!-> Affected by Z flux
				uR(i,j,k,QMAGY,2) = uR(i,j,k,QMAGY,2) + deltat/(3.d0*deltay)*&
										(E(i,j,k,1,2,2) - E(i,j,k,1,2,1))
				uR(i,j,k,QMAGX,2) = up(i,j,k,QMAGX,2) + deltat/(6.d0*deltay)*&
									((E(i,j,k,2,2,2) - E(i,j,k,2,2,1)) + &
									(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))

				!Z-Direction
				!-> Affected by X flux
				uR(i,j,k,QMAGZ,3) = up(i,j,k,QMAGZ,3) - deltat/(3.d0*deltaz)*&
										(E(i,j,k,2,1,2) - E(i,j,k,2,1,1))
				uR(i,j,k,QMAGY,3) = up(i,j,k,QMAGZ,3) - deltat/(6.d0*deltaz)*&
									((E(i,j,k,3,2,2) - E(i,j,k,3,2,1)) + &
									(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
				!-> Affected by Y flux
				uR(i,j,k,QMAGZ,3) = uR(i,j,k,QMAGZ,3) + deltat/(3.d0*deltaz)*&
										(E(i,j,k,1,2,2) - E(i,j,k,1,2,1))
				UR(i,j,k,QMAGX,3) = up(i,j,k,QMAGX,3) + deltat/(6.d0*deltaz)*&
									((E(i,j,k,3,2,2) - E(i,j,k,3,2,1)) + &
									(E(i,j,k,3,1,2) - E(i,j,k,3,1,1))
			enddo
		enddo
	enddo

end subroutine corner_couple_mag


subroutine electric_interp(Etemp, q, qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
			flx, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)

use amrex_fort_module, only : rt => amrex_real
use meth_mhd_params_module

implicit none
	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)

	real(rt), intent(out)	::Etemp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,3,2,2) !Three coordinates 4 edges
	
	real(rt)				::Ecen(3)
	real(rt)				::dxE_face(2), dxE_edge(2,2,2)
	real(rt)				::dyE_face(2), dyE_edge(2,2,2)
	real(rt)				::dzE_face(2), dzE_edge(2,2,2)
	
	integer					::i,j,k	

!Interpolate Electric Fields to edge
	do k = q_l3,q_h3
		do j = q_l2, q_h2
			do i = q_l1, q_h1
				!X-direction
				!1/4 interpolation
				call electric(q(i,j,k,:),Ecen)
				dxE_face(1) = 2.d0*(-flx(i,j,k,QMAGY,3) - Ecen(1))
				dyE_face(1) = 2.d0*(-flx(i,j,k,QMAGZ,1) - Ecen(2))
 				dzE_face(1) = 2.d0*(-flx(i,j,k,QMAGX,2) - Ecen(3))
				!3/4 interp
				call electric(q(i+1,j,k,:),Ecen)
				dxE_face(2) = 2.d0*(Ecen(1) + flx(i,j,k,QMAGY,3))
				call electric(q(i,j+1,k,:),Ecen)
				dyE_face(2) = 2.d0*(Ecen(2) + flx(i,j,k,QMAGZ,1))
				call electric(q(i,j,k+1,:),Ecen)
				dzE_face(2) = 2.d0*(Ecen(3) + flx(i,j,k,QMAGX,2))
				if(q(i,j,k,QU).gt. 0.d0) then
					dxE_edge(1, 1, 1) = dxE_face(1)
				elseif(q(i,j,k,QU).lt. 0.d0) then 
					dxE_edge(1, 1, 1) = dxE_face(2)
				else
					dxE_edge(1, 1, 1) = 0.5d0*(dxE_face(1) + dxE_face(2))
				endif
			enddo
		enddo
	enddo

	

end subroutine electric_interp	

end module ct_updwind
