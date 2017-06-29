module mhd_plm_module
!Module that gives a piecewise linear interpolation for the primitive variables 
!They are projected onto the characteristic variables for tracing. 
implicit none

  private vanleer, lvecx, lvecy, lvecz, rvecx, rvecy, rvecz, evals

  public plm 

contains 
  !
  ! characteristics based on u
  !
  !===========================================================================
  ! This is called from within threaded loops in advance_mhd_tile so *no* OMP here ...
  !===========================================================================

	subroutine plm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3, &
                   u,cspd,qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3, &
                   Ip,Im, ilo1,ilo2,ihi1,ihi2,dx,dy,dz,dt,k3d,kc,a_old,n)
    use amrex_fort_module, only : rt => amrex_real

	implicit none
	integer	, intent(in   ) ::   s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer	, intent(in   ) ::  qd_l1,qd_l2,qd_l3,qd_h1,qd_h2,qd_h3
    integer , intent(in   ) ::  ilo1,ilo2,ihi1,ihi2,n
 
    real(rt), intent(in   ) ::      s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3)
    real(rt), intent(in   ) ::      u(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3,3)
    real(rt), intent(in   ) ::   cspd(qd_l1:qd_h1,qd_l2:qd_h2,qd_l3:qd_h3)

    real(rt), intent(inout) :: Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)
    real(rt), intent(inout) :: Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,1:2,1:3,1:3)

    real(rt), intent(in   ) :: dx,dy,dz,dt,a_old
	real(rt) 				:: dW(QVAR)

    real(rt) :: dt_over_a
    integer          :: k3d,kc

    dt_over_a = dt / a_old

	end subroutine plm


!========================================= VanLeer TVD slope limiter =======================================
	subroutine vanleer(dW, WR, WL) 
	
	use amrex_fort_mudle, only : rt => amrex_real
	
	implicit none
	
	real(rt), intent(in )	::  WR, WL
	real(rt), intent(out)	::  dW
	dW = 0.0d0 	
	
	if( WR*WL .gt. 0.0d0 ) then 
	dW = 2.0d0*WR*WL/(WR + WL)

	end subroutine


!=========================================== Evals =========================================================

	subroutine evals(lam, Q)

	use amrex_fort_mudle, only : rt => amrex_real

	implicit none
	
	real(rt), intent(in)	:: Q(QVAR)
	real(rt), intent(out)	:: lam(3,7) !Three dimensions 7 waves

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Magnetic
	ca = (Q(QBX)**2 + Q(QBY)**2 + Q(QBZ)**2)/Q(QRHO)
	cax = (Q(QBX)**2)/Q(RHO)
	cay = (Q(QBY)**2)/Q(RHO)
	caz = (Q(QBZ)**2)/Q(RHO)
	!Sloooooooooow
	csx = 0.5d0*((a + ca) - sqrt((as + ca)**2 - 4.0d0*a*cax)
	csy = 0.5d0*((a + ca) - sqrt((as + ca)**2 - 4.0d0*a*cay)
	csz = 0.5d0*((a + ca) - sqrt((as + ca)**2 - 4.0d0*a*caz)
	!Fassssst
	cfx = 0.5d0*((a + ca) + sqrt((as + ca)**2 - 4.0d0*a*cax)
	cfy = 0.5d0*((a + ca) + sqrt((as + ca)**2 - 4.0d0*a*cay)
	cfz = 0.5d0*((a + ca) + sqrt((as + ca)**2 - 4.0d0*a*caz)
	
	!Ax eigenvalues
	lam(1,1) = Q(QU) - cfx
	lam(1,2) = Q(QU) - cax
	lam(1,3) = Q(QU) - csx
	lam(1,4) = Q(QU)
	lam(1,5) = Q(QU) + csx 
	lam(1,6) = Q(QU) + cax
	lam(1,7) = Q(QU) + cfx

	!Ay eigenvalues
	lam(2,1) = Q(QV) - cfy
	lam(2,2) = Q(QV) - cay
	lam(2,3) = Q(QV) - csy
	lam(2,4) = Q(QV)
	lam(2,5) = Q(QV) + csy 
	lam(2,6) = Q(QV) + cay
	lam(2,7) = Q(QV) + cfy

	!Az eigenvalues
	lam(3,1) = Q(QW) - cfz
	lam(3,2) = Q(QW) - caz
	lam(3,3) = Q(QW) - csz
	lam(3,4) = Q(QW)
	lam(3,5) = Q(QW) + csz 
	lam(3,6) = Q(QW) + caz
	lam(3,7) = Q(QW) + cfz
	end subroutine evals
	
!====================================== Left Eigenvectors ===============================================

!x direction
	subroutine lvecx(leig, Q) 
	use amrex_fort_mudle, only : rt => amrex_real
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cax, csx, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Magnetic
	ca = (Q(QBX)**2 + Q(QBY)**2 + Q(QBZ)**2)/Q(QRHO)
	cax = (Q(QBX)**2)/Q(RHO)
	!Sloooooooooow
	csx = 0.5d0*((a + ca) - sqrt((as + ca)**2 - 4.0d0*a*cax)
	!Fassssst
	cfx = 0.5d0*((a + ca) + sqrt((as + ca)**2 - 4.0d0*a*cax)
	!useful constants
	alf = (as - csx)/(cfx - csx)
	als = (cfx - as)/(cfx - csx)
	bety = Q(QBY)/(sqrt(Q(QBY)**2 + Q(QBZ)**2))
	betz = Q(QBZ)/(sqrt(Q(QBY)**2 + Q(QBZ)**2))
	cff = cfx*alf
	css = csx*als
	S = sign(1.0d0, Q(QBX))
	Qf = cfx*alf*S
	Qs = csx*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
	leig(1,:) = (/0.d0, -N*Cff, N*Qs*bety, N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/)
	leig(2,:) = (/0.d0,  0.d0, -0.5d0*betz, 0.5d0*bety, 0.d0, -0.5d0*betz/(sqrt(Q(QRHO))), 0.5d0*bety*S/(sqrt(Q(QRHO)))/)
	leig(3,:) = (/0.d0, -N*Css, -N*Qf*bety, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/)
	leig(4,:) = (/1.d0,  0.d0,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/)
	leig(5,:) = (/0.d0,  N*Css, N*Qf*bety, N*Qf*betz, N*als/Q(QRHO), -N*AAf*bety/Q(QRHO), -N*AAf*betz/Q(QRHO)/)
	leig(6,:) = (/0.d0,  0.d0, 0.5d0 betz, -0.5d0*bety, 0.d0, -0.5d0*betz*S/(sqrt(Q(QRHO))), 0.5d0*bety*S/(sqrt(Q(QRHO)))/)
	leig(7,:) = (/0.d0, N*Cff, -N*Qs*bety, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*bety/Q(QRHO), N*AAs*betz/Q(QRHO)/)
	
	
	end subroutine lvecx

!y direction
	subroutine lvecy(leig, Q) 
	use amrex_fort_mudle, only : rt => amrex_real
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfy, cay, csy, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Magnetic
	ca = (Q(QBX)**2 + Q(QBY)**2 + Q(QBZ)**2)/Q(QRHO)
	cay = (Q(QBY)**2)/Q(RHO)
	!Sloooooooooow
	csy = 0.5d0*((a + ca) - sqrt((as + ca)**2 - 4.0d0*a*cay)
	!Fassssst
	cfy = 0.5d0*((a + ca) + sqrt((as + ca)**2 - 4.0d0*a*cay)
	!useful constants
	alf = (as - csy)/(cfy - csy)
	als = (cfy - as)/(cfy - csy)
	betx = Q(QBX)/(sqrt(Q(QBX)**2 + Q(QBZ)**2))
	betz = Q(QBZ)/(sqrt(Q(QBX)**2 + Q(QBZ)**2))
	cff = cfy*alf
	css = csy*als
	S = sign(1.0d0, Q(QBY))
	Qf = cfy*alf*S
	Qs = csy*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
!Need to double check the rows
	leig(1,:) = (/0.d0, -N*Cff, N*Qs*betx, N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/)
	leig(2,:) = (/0.d0,  0.d0, -0.5d0*betz, 0.5d0*betx, 0.d0, -0.5d0*betz/(sqrt(Q(QRHO))), 0.5d0*betx*S/(sqrt(Q(QRHO)))/)
	leig(3,:) = (/0.d0, -N*Css, -N*Qf*betx, -N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/)
	leig(4,:) = (/1.d0,  0.d0,  0.d0, 0.d0, -1.d0/as, 0.d0, 0.d0/)
	leig(5,:) = (/0.d0,  N*Css, N*Qf*betx, N*Qf*betz, N*als/Q(QRHO), -N*AAf*betx/Q(QRHO), -N*AAf*betz/Q(QRHO)/)
	leig(6,:) = (/0.d0,  0.d0, 0.5d0 betz, -0.5d0*betx, 0.d0, -0.5d0*betz*S/(sqrt(Q(QRHO))), 0.5d0*betx*S/(sqrt(Q(QRHO)))/)
	leig(7,:) = (/0.d0, N*Cff, -N*Qs*betx, -N*Qs*betz, N*alf/Q(QRHO), N*AAs*betx/Q(QRHO), N*AAs*betz/Q(QRHO)/)
	
	
	end subroutine lvecy
end module mhd_plm_module
