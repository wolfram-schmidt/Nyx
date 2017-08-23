module ct_upwind

 use amrex_fort_module, only : rt => amrex_real
 use hlld_solver, only : hlld
 use meth_params_module

 implicit none 

 private primtocons
 public corner_transport
 
 
interface checkisnan
       module procedure checkisnanmult
       module procedure checkisnans
end interface checkisnan

 contains

subroutine corner_transport( q, qm, qp, q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3, &	
							flx, elec, flx_l1 , flx_l2 , flx_l3 , flx_h1 , flx_h2 , flx_h3, dx,dy,dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : QVAR
 use electric_field
implicit none

	integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3

	real(rt), intent(in)  :: q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR) !prim vars at time t^n
	real(rt), intent(in)  :: qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)  :: qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)	!Half Step Fluxes
!	real(rt), intent(out) :: elec(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,3,4) !Half Step Electric Fields
	real(rt), intent(out) :: elec(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,3) !Half Step Electric Fields
 
	real(rt)			  :: um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !PtoC Vars
	real(rt)			  :: up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt)			  :: cons_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Conservative Vars
	real(rt)			  :: cons_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Conservative Vars
	real(rt)			  :: cons_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Conservative Vars
	real(rt)			  :: cons_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt)			  :: q_temp_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Primitive Vars
	real(rt)			  :: q_temp_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2) !2D Temporary Primitive Vars
	real(rt)			  :: q_half_M(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Primitive Vars
	real(rt)			  :: q_half_P(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3) !Flux Corrected Primitive Vars
	real(rt)			  :: flx1D(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3) !Flux1d for all directions
	real(rt)			  :: flx2D(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3, 2) !Flux2d for all directions 2 perpendicular directions
!	real(rt) 			  :: Etemp(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,3,4) !Temporary Electric Field
	real(rt) 			  :: Etemp(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,3) !Temporary Electric Field
	real(rt)			  :: q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt) 			  :: dx, dy, dz, dt
	integer				  :: i, j , k


um = 0.d0
up = 0.d0
cons_temp_M = 0.d0
cons_temp_P = 0.d0
q_temp_M = 0.d0
q_temp_P = 0.d0
cons_half_M = 0.d0
cons_half_P = 0.d0
q_half_M = 0.d0
q_half_P = 0.d0
q2D = 0.d0
flx1D = 0.d0
flx2D = 0.d0


!Calculate Flux 1D
write(*,*) "Do Flux 1D"
	!x-dir
	call hlld(qm(:,:,:,:,1),qp(:,:,:,:,1),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx1D(:,:,:,:,1),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	!y-dir	
	call hlld(qm(:,:,:,:,2),qp(:,:,:,:,2),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx1D(:,:,:,:,2),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	!z-dir
	call hlld(qm(:,:,:,:,3),qp(:,:,:,:,3),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx1D(:,:,:,:,3),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)
!Prim to Cons
do i = 1,3
	call PrimToCons(qm(:,:,:,:,i), um(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
	call PrimToCons(qp(:,:,:,:,i), up(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo
!Use "1D" fluxes To interpolate Temporary Edge Centered Electric Fields
write(*,*) "Do Electric Field 1D"
	call electric_edge(Etemp, q, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
			flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)

write(*,*) "Corner Couple Cons"
!Corner Couple
	call corner_couple(cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
					   flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, dx, dy, dz, dt) !Correct Conservative vars using Transverse Fluxes
write(*,*) "Corner Couple Mag"
	call corner_couple_mag(cons_temp_M, cons_temp_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
						Etemp, dx, dy, dz, dt) !Correct Magnetic Vars using Etemp


!Cons To Prim
do i = 1,3
	call ConsToPrim(q_temp_M(:,:,:,:,i,1), cons_temp_M(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
	call ConsToPrim(q_temp_P(:,:,:,:,i,1), cons_temp_P(:,:,:,:,i,1), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
	call ConsToPrim(q_temp_M(:,:,:,:,i,2), cons_temp_M(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
	call ConsToPrim(q_temp_P(:,:,:,:,i,2), cons_temp_P(:,:,:,:,i,2), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

print *, "Flux 2D"
!Calculate Flux 2D
do i = 1,2
	print *, "x flux", i
	!x-dir
	call hlld(q_temp_M(:,:,:,:,1,i),q_temp_P(:,:,:,:,1,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx2D(:,:,:,:,1,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	print *, "y flux", i
	!y-dir	
	call hlld(q_temp_M(:,:,:,:,2,i),q_temp_P(:,:,:,:,2,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx2D(:,:,:,:,2,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	print *, "z flux", i
	!z-dir
	call hlld(q_temp_M(:,:,:,:,3,i),q_temp_P(:,:,:,:,3,i),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx2D(:,:,:,:,3,i),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)
enddo

!Use Averaged 2D fluxes to interpolate temporary Edge Centered Electric Fields, reuse "flx1D"
	flx1D(:,:,:,:,:) = 0.5d0*(flx2D(:,:,:,:,:,1) + flx2D(:,:,:,:,:,2))
	call electric_edge(Etemp, q, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
			flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)


!Half Step conservative vars
	call half_step(cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
			flx2D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, dx, dy, dz, dt)
	call half_step_mag(cons_half_M, cons_half_P, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
			   Etemp, dx, dy, dz, dt)
do i = 1,3
	call ConsToPrim(q_half_M(:,:,:,:,i), cons_half_M(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
	call ConsToPrim(q_half_P(:,:,:,:,i), cons_half_P(:,:,:,:,i), q_l1 , q_l2 , q_l3 , q_h1 , q_h2 , q_h3)
enddo

!Final Fluxes
!flx = 0.d0
	!x-dir
print *, "Flux Half Step"
	call hlld(q_half_M(:,:,:,:,1),q_half_P(:,:,:,:,1),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx(:,:,:,:,1),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 1)
	!y-dir	
	call hlld(q_half_M(:,:,:,:,2),q_half_P(:,:,:,:,2),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx(:,:,:,:,2),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 2)
	!z-dir
	call hlld(q_half_M(:,:,:,:,3),q_half_P(:,:,:,:,3),q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx(:,:,:,:,3),&
			  flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, 3)

	
!Primitive update
call prim_half(q2D,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx1D,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, dx, dy, dz, dt)
!Final Electric Field Update
call electric_edge(Etemp, q2D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
			flx1D, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)


end subroutine corner_transport

!================================================= Calculate the Conservative Variables ===============================================

subroutine PrimToCons(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt), intent(out)	::u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	integer					:: i ,j ,k

 do k = q_l3,q_h3
 	 do j = q_l2,q_h2
 		 do i = q_l1, q_h1
			u(i,j,k,URHO)  = q(i,j,k,QRHO)
 			u(i,j,k,UMX)    = q(i,j,k,QRHO)*q(i,j,k,QU)
 			u(i,j,k,UMY)    = q(i,j,k,QRHO)*q(i,j,k,QV)
 			u(i,j,k,UMZ)    = q(i,j,k,QRHO)*q(i,j,k,QW)
			u(i,j,k,UEDEN) = (q(i,j,k,QPRES)-0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ))))/(gamma_minus_1) &
							  + 0.5d0*q(i,j,k,QRHO)*dot_product(q(i,j,k,QU:QW),q(i,j,k,QU:QW)) &
							  + 0.5d0*(dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ)))
			u(i,j,k,UEINT) = (q(i,j,k,QPRES) - 0.5d0*dot_product(q(i,j,k,QMAGX:QMAGZ),q(i,j,k,QMAGX:QMAGZ)))/(gamma_minus_1)
			u(i,j,k,QMAGX:QMAGZ) = q(i,j,k,QMAGX:QMAGZ)
		 enddo
	 enddo
 enddo
end subroutine PrimToCons

!================================================= Calculate the Primitve Variables ===============================================

subroutine ConsToPrim(q, u, q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(in)	::u(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt), intent(out)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	integer					:: i ,j ,k

	q = u
 do k = q_l3,q_h3
 	 do j = q_l2,q_h2
 		 do i = q_l1, q_h1
 			q(i,j,k,QRHO)  = u(i,j,k,URHO)
 			q(i,j,k,QU)    = u(i,j,k,UMX)/q(i,j,k,QRHO)
 			q(i,j,k,QV)    = u(i,j,k,UMY)/q(i,j,k,QRHO)
 			q(i,j,k,QW)    = u(i,j,k,UMZ)/q(i,j,k,QRHO)
			q(i,j,k,QPRES) = u(i,j,k,UEINT)*gamma_minus_1 + 0.5*dot_product(u(i,j,k,QMAGX:QMAGZ),u(i,j,k,QMAGX:QMAGZ))
			q(i,j,k,QREINT) = (q(i,j,k,QPRES) - 0.5d0**dot_product(u(i,j,k,QMAGX:QMAGZ),u(i,j,k,QMAGX:QMAGZ)))/gamma_minus_1
			q(i,j,k,QMAGX:QMAGZ) = u(i,j,k,QMAGX:QMAGZ)
		 enddo
	 enddo
 enddo
end subroutine ConsToPrim

!======================================= Update the Temporary Conservative Variables with Transverse 1D Fluxes ========================
subroutine corner_couple(uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
					   flx, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module

implicit none
	
	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	
	real(rt), intent(in)	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)

	real(rt), intent(out)	::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
	real(rt), intent(out)	::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)

	real(rt)				:: dx, dy, dz, dt, u, v, w
	integer					:: i ,j ,k
	

	uL(:,:,:,:,:,1) = um
	uL(:,:,:,:,:,2) = um
	uR(:,:,:,:,:,1) = up
	uR(:,:,:,:,:,2) = up
	do k = flx_l3,flx_h3 -1
		do j = flx_l2,flx_h2 -1
			do i = flx_l1,flx_h1 -1
	!Left Corrected States
				uL(i,j,k,URHO:UEDEN,1,1) = um(i,j,k,URHO:UEDEN,1) - dt/(3.d0*dx)*(flx(i,j+1,k,URHO:UEDEN,2) - flx(i,j,k,URHO:UEDEN,2))! y corrected x
				uL(i,j,k,URHO:UEDEN,1,2) = um(i,j,k,URHO:UEDEN,1) - dt/(3.d0*dx)*(flx(i,j,k+1,URHO:UEDEN,3) - flx(i,j,k,URHO:UEDEN,3))! z corrected x
				uL(i,j,k,URHO:UEDEN,2,1) = um(i,j,k,URHO:UEDEN,2) - dt/(3.d0*dy)*(flx(i+1,j,k,URHO:UEDEN,1) - flx(i,j,k,URHO:UEDEN,1))! x corrected y
				uL(i,j,k,URHO:UEDEN,2,2) = um(i,j,k,URHO:UEDEN,2) - dt/(3.d0*dy)*(flx(i,j,k+1,URHO:UEDEN,3) - flx(i,j,k,URHO:UEDEN,3))! z corrected y
				uL(i,j,k,URHO:UEDEN,3,1) = um(i,j,k,URHO:UEDEN,3) - dt/(3.d0*dz)*(flx(i+1,j,k,URHO:UEDEN,1) - flx(i,j,k,URHO:UEDEN,1))! x corrected z
				uL(i,j,k,URHO:UEDEN,3,2) = um(i,j,k,URHO:UEDEN,3) - dt/(3.d0*dz)*(flx(i,j+1,k,URHO:UEDEN,2) - flx(i,j,k,URHO:UEDEN,2))! y corrected z

				u = uL(i,j,k,UMX,1,1)/uL(i,j,k,URHO,1,1)
				v = uL(i,j,k,UMY,1,1)/uL(i,j,k,URHO,1,1)
				w = uL(i,j,k,UMZ,1,1)/uL(i,j,k,URHO,1,1)
				uL(i,j,k,UEINT,1,1) = uL(i,j,k,UEDEN,1,1) - 0.5d0*uL(i,j,k,URHO,1,1)*(u**2 + v**2 + w**2)

				u = uL(i,j,k,UMX,1,2)/uL(i,j,k,URHO,1,2)
				v = uL(i,j,k,UMY,1,2)/uL(i,j,k,URHO,1,2)
				w = uL(i,j,k,UMZ,1,2)/uL(i,j,k,URHO,1,2)
				uL(i,j,k,UEINT,1,2) = uL(i,j,k,UEDEN,1,2) - 0.5d0*uL(i,j,k,URHO,1,2)*(u**2 + v**2 + w**2)
			
				u = uL(i,j,k,UMX,2,1)/uL(i,j,k,URHO,2,1)
				v = uL(i,j,k,UMY,2,1)/uL(i,j,k,URHO,2,1)
				w = uL(i,j,k,UMZ,2,1)/uL(i,j,k,URHO,2,1)
				uL(i,j,k,UEINT,2,1) = uL(i,j,k,UEDEN,2,1) - 0.5d0*uL(i,j,k,URHO,2,1)*(u**2 + v**2 + w**2)

				u = uL(i,j,k,UMX,2,2)/uL(i,j,k,URHO,2,2)
				v = uL(i,j,k,UMY,2,2)/uL(i,j,k,URHO,2,2)
				w = uL(i,j,k,UMZ,2,2)/uL(i,j,k,URHO,2,2)
				uL(i,j,k,UEINT,2,2) = uL(i,j,k,UEDEN,2,2) - 0.5d0*uL(i,j,k,URHO,2,2)*(u**2 + v**2 + w**2)

				u = uL(i,j,k,UMX,3,1)/uL(i,j,k,URHO,3,1)
				v = uL(i,j,k,UMY,3,1)/uL(i,j,k,URHO,3,1)
				w = uL(i,j,k,UMZ,3,1)/uL(i,j,k,URHO,3,1)
				uL(i,j,k,UEINT,3,1) = uL(i,j,k,UEDEN,3,1) - 0.5d0*uL(i,j,k,URHO,3,1)*(u**2 + v**2 + w**2)

				u = uL(i,j,k,UMX,3,2)/uL(i,j,k,URHO,3,2)
				v = uL(i,j,k,UMY,3,2)/uL(i,j,k,URHO,3,2)
				w = uL(i,j,k,UMZ,3,2)/uL(i,j,k,URHO,3,2)
				uL(i,j,k,UEINT,3,2) = uL(i,j,k,UEDEN,3,2) - 0.5d0*uL(i,j,k,URHO,3,2)*(u**2 + v**2 + w**2)



	!Right Corrected States
				uR(i,j,k,URHO:UEDEN,1,1) = up(i+1,j,k,URHO:UEDEN,1) - dt/(3.d0*dx)*(flx(i,j+1,k,URHO:UEDEN,2) - flx(i,j,k,URHO:UEDEN,2))
				uR(i,j,k,URHO:UEDEN,1,2) = up(i+1,j,k,URHO:UEDEN,1) - dt/(3.d0*dx)*(flx(i,j,k+1,URHO:UEDEN,3) - flx(i,j,k,URHO:UEDEN,3))
				uR(i,j,k,URHO:UEDEN,2,1) = up(i,j+1,k,URHO:UEDEN,2) - dt/(3.d0*dy)*(flx(i+1,j,k,URHO:UEDEN,1) - flx(i,j,k,URHO:UEDEN,1))
				uR(i,j,k,URHO:UEDEN,2,2) = up(i,j+1,k,URHO:UEDEN,2) - dt/(3.d0*dy)*(flx(i,j,k+1,URHO:UEDEN,3) - flx(i,j,k,URHO:UEDEN,3))
				uR(i,j,k,URHO:UEDEN,3,1) = up(i,j,k+1,URHO:UEDEN,3) - dt/(3.d0*dz)*(flx(i+1,j,k,URHO:UEDEN,1) - flx(i,j,k,URHO:UEDEN,1))
				uR(i,j,k,URHO:UEDEN,3,2) = up(i,j,k+1,URHO:UEDEN,3) - dt/(3.d0*dz)*(flx(i,j+1,k,URHO:UEDEN,2) - flx(i,j,k,URHO:UEDEN,2))
				
				u = uR(i,j,k,UMX,1,1)/uR(i,j,k,URHO,1,1)
				v = uR(i,j,k,UMY,1,1)/uR(i,j,k,URHO,1,1)
				w = uR(i,j,k,UMZ,1,1)/uR(i,j,k,URHO,1,1)
				uR(i,j,k,UEINT,1,1) = uR(i,j,k,UEDEN,1,1) - 0.5d0*uR(i,j,k,URHO,1,1)*(u**2 + v**2 + w**2)

				u = uR(i,j,k,UMX,1,2)/uR(i,j,k,URHO,1,2)
				v = uR(i,j,k,UMY,1,2)/uR(i,j,k,URHO,1,2)
				w = uR(i,j,k,UMZ,1,2)/uR(i,j,k,URHO,1,2)
				uR(i,j,k,UEINT,1,2) = uR(i,j,k,UEDEN,1,2) - 0.5d0*uR(i,j,k,URHO,1,2)*(u**2 + v**2 + w**2)
			
				u = uR(i,j,k,UMX,2,1)/uR(i,j,k,URHO,2,1)
				v = uR(i,j,k,UMY,2,1)/uR(i,j,k,URHO,2,1)
				w = uR(i,j,k,UMZ,2,1)/uR(i,j,k,URHO,2,1)
				uR(i,j,k,UEINT,2,1) = uR(i,j,k,UEDEN,2,1) - 0.5d0*uR(i,j,k,URHO,2,1)*(u**2 + v**2 + w**2)

				u = uR(i,j,k,UMX,2,2)/uR(i,j,k,URHO,2,2)
				v = uR(i,j,k,UMY,2,2)/uR(i,j,k,URHO,2,2)
				w = uR(i,j,k,UMZ,2,2)/uR(i,j,k,URHO,2,2)
				uR(i,j,k,UEINT,2,2) = uR(i,j,k,UEDEN,2,2) - 0.5d0*uR(i,j,k,URHO,2,2)*(u**2 + v**2 + w**2)

				u = uR(i,j,k,UMX,3,1)/uR(i,j,k,URHO,3,1)
				v = uR(i,j,k,UMY,3,1)/uR(i,j,k,URHO,3,1)
				w = uR(i,j,k,UMZ,3,1)/uR(i,j,k,URHO,3,1)
				uR(i,j,k,UEINT,3,1) = uR(i,j,k,UEDEN,3,1) - 0.5d0*uR(i,j,k,URHO,3,1)*(u**2 + v**2 + w**2)

				u = uR(i,j,k,UMX,3,2)/uR(i,j,k,URHO,3,2)
				v = uR(i,j,k,UMY,3,2)/uR(i,j,k,URHO,3,2)
				w = uR(i,j,k,UMZ,3,2)/uR(i,j,k,URHO,3,2)
				uR(i,j,k,UEINT,3,2) = uR(i,j,k,UEDEN,3,2) - 0.5d0*uR(i,j,k,URHO,3,2)*(u**2 + v**2 + w**2)

			enddo
		enddo
	enddo
end subroutine corner_couple

!================================== Use 1D Electric Fields to Transverse correct the Temporary Magnetic Fields ===========================
subroutine corner_couple_mag(uL, uR, um, up,&
                            q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
       						E, dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : QVAR, QMAGX, QMAGY, QMAGZ

!Correction using Faraday's Law
implicit none

	integer, intent(in)			::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(inout)		::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
	real(rt), intent(inout)     ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3,2)
	real(rt), intent(in)    	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)    	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)	
	real(rt), intent(in)		::E(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,3)

	real(rt)					::dx, dy, dz, dt
	
	integer						:: i ,j ,k
   	do k = q_l3,q_h3-1
		do j = q_l2,q_h2-1
			do i = q_l1,q_h1-1	
		!Left State
				!X-direction 
				!-> Affected by Y flux
				uL(i,j,k,QMAGX,1,1) = um(i,j,k,QMAGX,1) + dt/(3.d0*dx)*&
										(E(i,j+1,k,3) - E(i,j,k,3))
				uL(i,j,k,QMAGZ,1,1) = um(i,j,k,QMAGZ,1) + dt/(6.d0*dx)*&
										((E(i,j+1,k+1,1) - E(i,j,k+1,1)) + &
										(E(i,j+1,k,1) - E(i,j,k,1)))
				!-> Affected by Z flux
				uL(i,j,k,QMAGX,1,2) = um(i,j,k,QMAGX,1) - dt/(3.d0*dx)*&
									 	(E(i,j,k+1,2) - E(i,j,k,2))
				uL(i,j,k,QMAGY,1,1) = um(i,j,k,QMAGY,1) - dt/(6.d0*dx)*&
									 	((E(i,j+1,k+1,1) - E(i,j+1,k,1)) + &
										(E(i,j,k+1,1) - E(i,j,k,1)))

				uL(i,j,k,QMAGY:QMAGZ,1,2) = um(i,j,k,QMAGY:QMAGZ,1)
				!Y-direction
				!-> Affected by X flux
				uL(i,j,k,QMAGY,2,1) = um(i,j,k,QMAGY,2) - dt/(3.d0*dy)*&
										(E(i+1,j,k,3) - E(i,j,k,3))
				uL(i,j,k,QMAGZ,2,1) = um(i,j,k,QMAGZ,2) - dt/(6.d0*dy)*&
									((E(i+1,j,k+1,2) - E(i,j,k+1,2)) + &
									(E(i+1,j,k,2) - E(i,j,k,2)))
				!-> Affected by Z flux
				uL(i,j,k,QMAGY,2,2) = um(i,j,k,QMAGY,2) + dt/(3.d0*dy)*&
										(E(i,j,k+1,1) - E(i,j,k,1))
				uL(i,j,k,QMAGX,2,1) = um(i,j,k,QMAGX,2) + dt/(6.d0*dy)*&
									((E(i+1,j,k+1,2) - E(i+1,j,k,2)) + &
									(E(i,j,k+1,2) - E(i,j,k,2)))

				uL(i,j,k,QMAGX,2,2) = um(i,j,k,QMAGX,2)
				uL(i,j,k,QMAGZ,2,2) = um(i,j,k,QMAGZ,2)

				!Z-Direction
				!-> Affected by X flux
				uL(i,j,k,QMAGZ,3,1) = um(i,j,k,QMAGZ,3) + dt/(3.d0*dz)*&
										(E(i+1,j,k,2) - E(i,j,k,2))
				uL(i,j,k,QMAGY,3,1) = um(i,j,k,QMAGZ,3) + dt/(6.d0*dz)*&
									((E(i+1,j+1,k,3) - E(i,j+1,k,3)) + &
									(E(i+1,j,k,3) - E(i,j,k,3)))
				!-> Affected by Y flux
				uL(i,j,k,QMAGZ,3,2) = um(i,j,k,QMAGZ,3) - dt/(3.d0*dz)*&
										(E(i,j+1,k,1) - E(i,j,k,1))
				uL(i,j,k,QMAGX,3,1) = um(i,j,k,QMAGX,3) - dt/(6.d0*dz)*&
									((E(i+1,j+1,k,3) - E(i+1,j,k,3)) + &
									(E(i,j+1,k,3) - E(i,j,k,3)))

				uL(i,j,k,QMAGX:QMAGY,3,2) = um(i,j,k,QMAGY:QMAGZ,3)
				uL(i,j,k,UEINT,1,1) = uL(i,j,k,UEDEN,1,1) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,1,1),uL(i,j,k,QMAGX:QMAGY,1,1))

				uL(i,j,k,UEINT,1,2) = uL(i,j,k,UEDEN,1,2) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,1,2),uL(i,j,k,QMAGX:QMAGY,1,2))

				uL(i,j,k,UEINT,2,1) = uL(i,j,k,UEDEN,2,1) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,2,1),uL(i,j,k,QMAGX:QMAGY,2,1))

				uL(i,j,k,UEINT,2,2) = uL(i,j,k,UEDEN,2,2) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,2,2),uL(i,j,k,QMAGX:QMAGY,2,2))

				uL(i,j,k,UEINT,3,1) = uL(i,j,k,UEDEN,3,1) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,3,1),uL(i,j,k,QMAGX:QMAGY,3,1))

				uL(i,j,k,UEINT,3,2) = uL(i,j,k,UEDEN,3,2) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,3,2),uL(i,j,k,QMAGX:QMAGY,3,2))
		!Right State
				!X-direction 
				!-> Affected by Y flux
				uR(i,j,k,QMAGX,1,1) = up(i,j,k,QMAGX,1) + dt/(3.d0*dx)*&
										(E(i+1,j+1,k,3) - E(i+1,j,k,3))
				uR(i,j,k,QMAGZ,1,1) = up(i,j,k,QMAGZ,1) + dt/(6.d0*dx)*&
										((E(i,j+1,k+1,1) - E(i,j,k+1,1)) + &
										(E(i,j+1,k,1) - E(i,j,k,1)))
				!-> Affected by Z flux
				uR(i,j,k,QMAGX,1,2) = up(i,j,k,QMAGX,1) - dt/(3.d0*dx)*&
									 	(E(i+1,j,k+1,2) - E(i+1,j,k,2))
				uR(i,j,k,QMAGY,1,1) = up(i,j,k,QMAGY,1) - dt/(6.d0*dx)*&
									 	((E(i,j+1,k+1,1) - E(i,j+1,k,1)) + &
										(E(i,j,k+1,1) - E(i,j,k,1)))

				uR(i,j,k,QMAGY:QMAGZ,1,2) = up(i,j,k,QMAGY:QMAGZ,1)
				!Y-direction
				!-> Affected by X flux
				uR(i,j,k,QMAGY,2,1) = up(i,j,k,QMAGY,2) - dt/(3.d0*dy)*&
										(E(i+1,j+1,k,3) - E(i,j+1,k,3))
				uR(i,j,k,QMAGZ,2,1) = up(i,j,k,QMAGZ,2) - dt/(6.d0*dy)*&
									((E(i+1,j,k+1,2) - E(i,j,k+1,2)) + &
									(E(i+1,j,k,2) - E(i,j,k,2)))
				!-> Affected by Z flux
				uR(i,j,k,QMAGY,2,2) = up(i,j,k,QMAGY,2) + dt/(3.d0*dy)*&
										(E(i,j+1,k+1,1) - E(i,j+1,k,1))
				uR(i,j,k,QMAGX,2,1) = up(i,j,k,QMAGX,2) + dt/(6.d0*dy)*&
									((E(i+1,j,k+1,2) - E(i+1,j,k,2)) + &
									(E(i,j,k+1,2) - E(i,j,k,2)))

				uR(i,j,k,QMAGX,2,2) = up(i,j,k,QMAGX,2)
				uR(i,j,k,QMAGZ,2,2) = up(i,j,k,QMAGZ,2)

				!Z-Direction
				!-> Affected by X flux
				uR(i,j,k,QMAGZ,3,1) = up(i,j,k,QMAGZ,3) + dt/(3.d0*dz)*&
										(E(i+1,j,k+1,2) - E(i,j,k+1,2))
				uR(i,j,k,QMAGY,3,1) = up(i,j,k,QMAGZ,3) + dt/(6.d0*dz)*&
									((E(i+1,j+1,k,3) - E(i,j+1,k,3)) + &
									(E(i+1,j,k,3) - E(i,j,k,3)))
				!-> Affected by Y flux
				uR(i,j,k,QMAGZ,3,2) = up(i,j,k,QMAGZ,3) - dt/(3.d0*dz)*&
										(E(i,j+1,k+1,1) - E(i,j,k+1,1))
				if(abs(uR(i,j,k,QMAGZ,3,2)).gt.1d-12.and.i.gt.0.and.j.gt.0.and.k.ge.0) then 
					print *, "Non trivial Bz = ",uR(i,j,k,QMAGZ,3,2) , "i j k = ", i, j, k
					print *, "Electric x",  E(i,j+1,k+1,1), E(i,j,k+1,1) 
					pause
				endif
					
				uR(i,j,k,QMAGX,3,1) = up(i,j,k,QMAGX,3) - dt/(6.d0*dz)*&
									((E(i+1,j+1,k,3) - E(i+1,j,k,3)) + &
									(E(i,j+1,k,3) - E(i,j,k,3)))

				uR(i,j,k,QMAGX:QMAGY,3,2) = up(i,j,k,QMAGY:QMAGZ,3)

				uR(i,j,k,UEINT,1,1) = uR(i,j,k,UEDEN,1,1) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,1,1),uR(i,j,k,QMAGX:QMAGY,1,1))

				uR(i,j,k,UEINT,1,2) = uR(i,j,k,UEDEN,1,2) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,1,2),uR(i,j,k,QMAGX:QMAGY,1,2))

				uR(i,j,k,UEINT,2,1) = uR(i,j,k,UEDEN,2,1) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,2,1),uR(i,j,k,QMAGX:QMAGY,2,1))

				uR(i,j,k,UEINT,2,2) = uR(i,j,k,UEDEN,2,2) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,2,2),uR(i,j,k,QMAGX:QMAGY,2,2))

				uR(i,j,k,UEINT,3,1) = uR(i,j,k,UEDEN,3,1) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,3,1),uR(i,j,k,QMAGX:QMAGY,3,1))

				uR(i,j,k,UEINT,3,2) = uR(i,j,k,UEDEN,3,2) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,3,2),uR(i,j,k,QMAGX:QMAGY,3,2))
			enddo
		enddo
	enddo

end subroutine corner_couple_mag

!====================================================== Final Conservative Corrections================================================================
subroutine half_step(uL, uR, um, up, q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,&
					   flx, flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : QVAR, URHO, UEDEN

implicit none
	
	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	
	real(rt), intent(in)	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3,2)

	real(rt), intent(out)	::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(out)	::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)

	real(rt)				:: dx, dy, dz, dt, u, v, w
	integer					:: i ,j ,k, n

	uL = um
	uR = up
   	do k = flx_l3,flx_h3-1
		do j = flx_l2,flx_h2-1
			do i = flx_l1,flx_h1-1
!left state				
				uL(i+1,j,k,URHO:UEDEN,1) = um(i+1,j,k,URHO:UEDEN,1)! - 0.5d0*dt/dx*(flx(i,j+1,k,URHO:UEDEN,2,2) - flx(i,j,k,URHO:UEDEN,2,2)) &
!											- 0.5d0*dt/dx*(flx(i,j,k+1,URHO:UEDEN,3,2) - flx(i,j,k,URHO:UEDEN,3,2))
				uL(i,j+1,k,URHO:UEDEN,2) = um(i,j+1,k,URHO:UEDEN,2)! - 0.5d0*dt/dy*(flx(i+1,j,k,URHO:UEDEN,1,2) - flx(i,j,k,URHO:UEDEN,1,2)) &
!											- 0.5d0*dt/dy*(flx(i,j,k+1,URHO:UEDEN,3,1) - flx(i,j,k,URHO:UEDEN,3,1))
				uL(i,j,k+1,URHO:UEDEN,3) = um(i,j,k+1,URHO:UEDEN,3)! - 0.5d0*dt/dz*(flx(i+1,j,k,URHO:UEDEN,1,1) - flx(i,j,k,URHO:UEDEN,1,1)) &
!											- 0.5d0*dt/dz*(flx(i,j+1,k,URHO:UEDEN,2,1) - flx(i,j,k,URHO:UEDEN,2,1))
				do n = 1,3
					u = uL(i,j,k,UMX,n)/uL(i,j,k,URHO,n)
					v = uL(i,j,k,UMY,n)/uL(i,j,k,URHO,n)
					w = uL(i,j,k,UMZ,n)/uL(i,j,k,URHO,n)
					uL(i,j,k,UEINT,n) = uL(i,j,k,UEDEN,n) - 0.5d0*uL(i,j,k,URHO,n)*(u**2 + v**2 + w**2)
				enddo
!right state				
				uR(i+1,j,k,URHO:UEDEN,1) = up(i+1,j,k,URHO:UEDEN,1)! - 0.5d0*dt/dx*(flx(i,j+1,k,URHO:UEDEN,2,2) - flx(i,j,k,URHO:UEDEN,2,2)) &
!											- 0.5d0*dt/dx*(flx(i,j,k+1,URHO:UEDEN,3,2) - flx(i,j,k,URHO:UEDEN,3,2))
				uR(i,j+1,k,URHO:UEDEN,2) = up(i,j+1,k,URHO:UEDEN,2)! - 0.5d0*dt/dy*(flx(i+1,j,k,URHO:UEDEN,1,2) - flx(i,j,k,URHO:UEDEN,1,2)) &
!											- 0.5d0*dt/dy*(flx(i,j,k+1,URHO:UEDEN,3,1) - flx(i,j,k,URHO:UEDEN,3,1))
				uR(i,j,k+1,URHO:UEDEN,3) = up(i,j,k+1,URHO:UEDEN,3)! - 0.5d0*dt/dz*(flx(i+1,j,k,URHO:UEDEN,1,1) - flx(i,j,k,URHO:UEDEN,1,1)) &
!											- 0.5d0*dt/dz*(flx(i,j+1,k,URHO:UEDEN,2,1) - flx(i,j,k,URHO:UEDEN,2,1))
				do n = 1,3
					u = uR(i,j,k,UMX,n)/uR(i,j,k,URHO,n)
					v = uR(i,j,k,UMY,n)/uR(i,j,k,URHO,n)
					w = uR(i,j,k,UMZ,n)/uR(i,j,k,URHO,n)
					uR(i,j,k,UEINT,n) = uR(i,j,k,UEDEN,n) - 0.5d0*uR(i,j,k,URHO,n)*(u**2 + v**2 + w**2)
				enddo
			enddo
		enddo
	enddo
end subroutine 

!================================================= Final Magnetic Corrections ========================================================================
subroutine half_step_mag(uL, uR, um, up, q_l1, q_l2, q_l3, q_h1, q_h2, q_h3, E, dx, dy, dz, dt)
use amrex_fort_module, only : rt => amrex_real
use meth_params_module, only : QVAR, QMAGX,QMAGY,QMAGZ

!Correction using Faraday's Law
implicit none

	integer, intent(in)			::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	real(rt), intent(inout)		::uL(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(inout)     ::uR(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)    	::um(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
	real(rt), intent(in)    	::up(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)	
	real(rt), intent(in)		::E(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,3)
	real(rt)					:: dx, dy, dz, dt
	integer						:: i ,j ,k, n

   	do k = q_l3,q_h3-1
		do j = q_l2,q_h2-1
			do i = q_l1,q_h1-1
		!---------------------------------------left state-----------------------------------------------------
				!X-Direction
				!Bx
				uL(i,j,k,QMAGX,1) = um(i,j,k,QMAGX,1) - 0.5d0*dt/dx*((E(i,j,k+1,2) - E(i,j,k,2)) &
								    - (E(i,j+1,k,3) - E(i,j,k,3)))
				!By
				uL(i,j,k,QMAGY,1) = um(i,j,k,QMAGY,1) - 0.5d0*dt/dx*((E(i,j+1,k+1,1) - E(i,j+1,k,1)) &
									+ (E(i,j,k+1,1) - E(i,j,k,1)) - (E(i+1,j+1,k,3) - E(i,j+1,k,3)) &
									- (E(i+1,j,k,3) - E(i,j,k,3)))
				!Bz
				uL(i,j,k,QMAGZ,1) = um(i,j,k,QMAGZ,1) + 0.5d0*dt/dx*((E(i,j+1,k+1,1) - E(i,j,k+1,1)) &
									+ (E(i,j+1,k,1) - E(i,j,k,1)) - (E(i+1,j,k+1,2) - E(i,j,k+1,2)) &
									- (E(i+1,j,k,2) - E(i,j,k,2)))
				if(i.eq.4.and.j.eq.14.and.k.eq.3) then 
					print *, "Un+1/2  x= ", uL(i,j,k,:,1)
					print *, "U in = ", um(i,j,k,:,1)
					print *, "E", E(i,j+1,k+1,1) , E(i,j,k+1,1), E(i,j+1,k,1) , E(i,j,k,1),  E(i+1,j,k+1,2) , E(i,j,k+1,2), E(i+1,j,k,2), E(i,j,k,2)
					pause
				endif
				!Y-Direction
				!Bx
				uL(i,j,k,QMAGX,2) = um(i,j,k,QMAGX,2) + 0.5d0*dt/dy*((E(i+1,j,k+1,2) - E(i+1,j,k,2)) &
								    + (E(i,j,k+1,2) - E(i,j,k,2)) - (E(i+1,j+1,k,3) - E(i+1,j,k,3)) &
									- (E(i,j+1,k,3) - E(i,j,k,3)))
				!By
				uL(i,j,k,QMAGY,2) = um(i,j,k,QMAGY,2) + 0.5d0*dt/dy*((E(i,j,k+1,1) - E(i,j,k,1)) &
								    - (E(i+1,j,k,3) - E(i,j,k,3)))
				!Bz
				uL(i,j,k,QMAGZ,2) = um(i,j,k,QMAGZ,2) - 0.5d0*dt/dy*((E(i+1,j,k+1,2) - E(i,j,k+1,2)) &
								    + (E(i+1,j,k,2) - E(i,j,k,2)) - (E(i,j+1,k+1,1) - E(i,j,k+1,1)) &
									- (E(i,j+1,k,1) - E(i,j,k,1)))		
				!Z-direction
				!Bx
				uL(i,j,k,QMAGX,3) = um(i,j,k,QMAGX,3) - 0.5d0*dt/dz*((E(i+1,j+1,k,3) - E(i+1,j,k,3)) &
								    + (E(i,j+1,k,3) - E(i,j,k,3)) - (E(i+1,j,k+1,2) - E(i+1,j,k,2)) &
									- (E(i,j+1,k,2) - E(i,j,k,2)))
				!By
				uL(i,j,k,QMAGY,3) = um(i,j,k,QMAGY,3) + 0.5d0*dt/dz*((E(i+1,j+1,k,3) - E(i,j+1,k,3)) &
								    +(E(i+1,j,k,3) - E(i,j,k,3)) - (E(i,j+1,k+1,1) - E(i,j+1,k,1)) &
									- (E(i,j,k+1,1) - E(i,j,k,1)))
				!Bz
				uL(i,j,k,QMAGZ,3) = um(i,j,k,QMAGZ,3) - 0.5d0*dt/dz*((E(i,j+1,k,1) - E(i,j,k,1)) &
									- (E(i+1,j,k,2) - E(i,j,k,2)))
				do n = 1,3
					uL(i,j,k,UEINT,n) = uL(i,j,k,UEDEN,n) - 0.5d0*dot_product(uL(i,j,k,QMAGX:QMAGY,n),uL(i,j,k,QMAGX:QMAGY,n))
				enddo

	!---------------------------------------right state-----------------------------------------------------
				!X-Direction
				!Bx
				uR(i,j,k,QMAGX,1) = up(i,j,k,QMAGX,1) - 0.5d0*dt/dx*((E(i+1,j,k+1,2) - E(i+1,j,k,2)) &
								    - (E(i+1,j+1,k,3) - E(i+1,j,k,3)))
				!By
				uR(i,j,k,QMAGY,1) = up(i,j,k,QMAGY,1) - 0.5d0*dt/dx*((E(i,j+1,k+1,1) - E(i,j+1,k,1)) &
									+ (E(i,j,k+1,1) - E(i,j,k,1)) - (E(i+1,j+1,k,3) - E(i,j+1,k,3)) &
									- (E(i+1,j,k,3) - E(i,j,k,3)))
				!Bz
				uR(i,j,k,QMAGZ,1) = up(i,j,k,QMAGZ,1) + 0.5d0*dt/dx*((E(i,j+1,k+1,1) - E(i,j,k+1,1)) &
									+ (E(i,j+1,k,1) - E(i,j,k,1)) - (E(i+1,j,k+1,2) - E(i,j,k+1,2)) &
									- (E(i+1,j,k,2) - E(i,j,k,2)))				
				!Y-Direction
				!Bx
				uR(i,j,k,QMAGX,2) = up(i,j,k,QMAGX,2) + 0.5d0*dt/dy*((E(i+1,j,k+1,2) - E(i+1,j,k,2)) &
								    + (E(i,j,k+1,2) - E(i,j,k,2)) - (E(i+1,j+1,k,3) - E(i+1,j,k,3)) &
									- (E(i,j+1,k,3) - E(i,j,k,3)))
				!By
				uR(i,j,k,QMAGY,2) = up(i,j,k,QMAGY,2) + 0.5d0*dt/dy*((E(i,j+1,k+1,1) - E(i,j+1,k,1)) &
								    - (E(i+1,j+1,k,3) - E(i,j+1,k,3)))
				!Bz
				uR(i,j,k,QMAGZ,2) = up(i,j,k,QMAGZ,2) - 0.5d0*dt/dy*((E(i+1,j,k+1,2) - E(i,j,k+1,2)) &
								    + (E(i+1,j,k,2) - E(i,j,k,2)) - (E(i,j+1,k+1,1) - E(i,j,k+1,1)) &
									- (E(i,j+1,k,1) - E(i,j,k,1)))		
				!Z-direction
				!Bx
				uR(i,j,k,QMAGX,3) = up(i,j,k,QMAGX,3) - 0.5d0*dt/dz*((E(i+1,j+1,k,3) - E(i+1,j,k,3)) &
								    +(E(i,j+1,k,3) - E(i,j,k,3)) - (E(i+1,j,k+1,2) - E(i+1,j,k,2)) &
									- (E(i,j,k+1,2) - E(i,j,k,2)))
				!By
				uR(i,j,k,QMAGY,3) = up(i,j,k,QMAGY,3) + 0.5d0*dt/dz*((E(i+1,j+1,k,3) - E(i,j+1,k,3)) &
								    +(E(i+1,j,k,3) - E(i,j,k,3)) - (E(i,j+1,k+1,1) - E(i,j+1,k,1)) &
									- (E(i,j,k+1,1) - E(i,j,k,1)))
				!Bz
				uR(i,j,k,QMAGZ,3) = up(i,j,k,QMAGZ,3) - 0.5d0*dt/dz*((E(i,j+1,k+1,1) - E(i,j,k+1,1)) &
									- (E(i+1,j,k+1,2) - E(i,j,k+1,2)))			
				do n = 1,3
					uR(i,j,k,UEINT,n) = uR(i,j,k,UEDEN,n) - 0.5d0*dot_product(uR(i,j,k,QMAGX:QMAGY,n),uR(i,j,k,QMAGX:QMAGY,n))
				enddo	
			enddo
		enddo
	enddo

end subroutine half_step_mag

!================================== Find the 2D corrected primitive variables =======================================

subroutine prim_half(q2D,q,q_l1,q_l2,q_l3,q_h1,q_h2,q_h3,flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3 &
					, dx, dy, dz, dt)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : QVAR

implicit none

	integer, intent(in)		::q_l1,q_l2,q_l3,q_h1,q_h2, q_h3
	integer, intent(in)		::flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	real(rt), intent(in)	::q(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)
	real(rt), intent(in) 	::flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR,3)
	real(rt), intent(out)	::q2D(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR)

	real(rt)				::flx_sum(QVAR)
	real(rt)				::qflx(QVAR)
	real(rt)				:: dx, dy, dz, dt	
	integer					::i, j, k
	q2D = q
	do k = flx_l3,flx_h3-1
		do j = flx_l2,flx_h2-1
			do i = flx_l1,flx_h1-1
				flx_sum = (flx(i+1,j,k,:,1) - flx(i,j,k,:,1)) + (flx(i,j+1,k,:,2) - flx(i,j,k,:,2)) + (flx(i,j,k+1,:,3) - flx(i,j,k,:,3)) 
				call qflux(qflx,flx_sum,q(i,j,k,:))
				q2D(i,j,k,:) = q(i,j,k,:) - 0.5d0*dt/dx*qflx
			enddo
		enddo
	enddo
end subroutine prim_half


!================================= Calculate the C to P Jacobian applied to the fluxes ===================================

subroutine qflux(qflx,flx,q)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only : QRHO, QU, QV, QW, QPRES, QMAGX, QMAGY, QMAGZ, QVAR, gamma_minus_1

implicit none

 real(rt), intent(in)		::flx(QVAR), q(QVAR)
 real(rt), intent(out)		::qflx(QVAR)
	qflx = 0.d0
	qflx(QRHO)  = flx(QRHO)
	qflx(QU)    = q(QU)*flx(QRHO) + q(QRHO)*flx(QU)
	qflx(QV)    = q(QV)*flx(QRHO) + q(QRHO)*flx(QV)
	qflx(QW)    = q(QW)*flx(QRHO) + q(QRHO)*flx(QW)
	qflx(QPRES) = flx(QRHO)*0.5d0*(q(QU)**2 + q(QV)**2 + q(QW)**2) + q(QRHO)*q(QU)*flx(QU) + q(QRHO)*q(QV)*flx(QV) &
				  + q(QRHO)*q(QW)*flx(QW) + 1.d0/gamma_minus_1*flx(QPRES) + q(QMAGX)*flx(QMAGX) + q(QMAGY)*flx(QMAGY) &
				  + q(QMAGZ)*flx(QMAGZ)
	qflx(QMAGX) = flx(QMAGX)
	qflx(QMAGY) = flx(QMAGY)
	qflx(QMAGZ) = flx(QMAGZ)

end subroutine qflux


!============================================ Debug code =====================================================
	subroutine checkisnanmult(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num)
	   use amrex_fort_module, only : rt => amrex_real

	implicit none
	integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3, num
	real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3,num)

	integer :: i,j,k,n


	do n = 1,num
		do k = uout_l3,uout_h3
			do j = uout_l2, uout_h2
				do i = uout_l1,uout_h1
					if(isnan(uout(i,j,k,n)).or.(abs(uout(i,j,k,n)).ge. 1d14)) then
						write(*,*) "Bad values ",  uout(i,j,k,:)
						write(*,*) "Failure to converge ", "i, j, k, n = ", i, j, k, n
						stop
					endif
				enddo
			enddo
		enddo
	enddo
	end subroutine checkisnanmult
!============ single =====================	

	subroutine checkisnans(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
	   use amrex_fort_module, only : rt => amrex_real

	implicit none
	integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
	real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

	integer :: i,j,k
!    write(*,*) "uout = ", uout 
!       pause

		do k = uout_l3,uout_h3
			do j = uout_l2, uout_h2
				do i = uout_l1,uout_h1
					if(isnan(uout(i,j,k)).or.(abs(uout(i,j,k)).ge. 1d14)) then
						write(*,*) "Bad values ",  uout(i,j,k)
						write(*,*) "Failure to converge ", "i, j, k = ", i, j, k
						stop
					endif
				enddo
			enddo
		enddo
	end subroutine checkisnans

!====================================== Density Check ========================================
	subroutine checknegdens(uout,uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3)
	   use amrex_fort_module, only : rt => amrex_real

	implicit none
	integer, intent(in)  :: uout_l1,uout_l2,uout_l3,uout_h1,uout_h2,uout_h3
	real(rt), intent(in) :: uout(uout_l1:uout_h1,uout_l2:uout_h2,uout_l3:uout_h3)

	integer :: i,j,k

		do k = uout_l3,uout_h3
			do j = uout_l2, uout_h2
				do i = uout_l1,uout_h1
					if(uout(i,j,k).le. 0.d0) then
						write(*,*) "Non-Positive Density ",  uout(i,j,k)
						write(*,*) "i, j, k = ", i, j, k
						stop
					endif
				enddo
			enddo
		enddo
	end subroutine checknegdens
end module ct_upwind
