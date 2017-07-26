module mhd_plm_module
!Module that gives a piecewise linear interpolation for the primitive variables 
!They are projected onto the characteristic variables for tracing. 

use meth_mhd_params_module
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

	subroutine plm(s,s_l1,s_l2,s_l3,s_h1,s_h2,s_h3,&	
				   bx, bxl1, bxl2, bxl3, bxh1, bxh2, bxh3, &
				   by, byl1, byl2, byl3, byh1, byh2, byh3, &
				   bz, bzl1, bzl2, bzl3, bzh1, bzh2, bzh3, &
                   Ip,Im, ilo1,ilo2,ilo3,ihi1,ihi2,ihi3,dx,dy,dz,dt,a_old)
    use amrex_fort_module, only : rt => amrex_real

	implicit none
	integer	, intent(in   ) ::  	s_l1, s_l2, s_l3, s_h1, s_h2, s_h3
    integer , intent(in   ) ::  	ilo1,ilo2,ilo3,ihi1,ihi2,ihi3
    integer , intent(in   ) ::  	bxl1, bxl2, bxl3, bxh1, bxh2, bxh3	
    integer , intent(in   ) ::  	byl1, byl2, byl3, byh1, byh2, byh3
    integer , intent(in   ) ::  	bzl1, bzl2, bzl3, bzh1, bzh2, bzh3
 
    real(rt), intent(in   ) ::      s( s_l1: s_h1, s_l2: s_h2, s_l3: s_h3, 1:QVAR) !Primitive Vars
    real(rt), intent(in   ) ::      bx(bxl1:bxh1, bxl2:bxh2, bxl3:bxh3)	!Face Centered Magnetic Fields
    real(rt), intent(in   ) ::      by(bxl1:bxh1, bxl2:bxh2, bxl3:bxh3)
    real(rt), intent(in   ) ::      bz(bxl1:bxh1, bxl2:bxh2, bxl3:bxh3)

    real(rt), intent(inout) :: 		Ip(ilo1-1:ihi1+1,ilo2-1:ihi2+1,ilo3-1:ihi3+1,1:QVAR,1:3)
    real(rt), intent(inout) :: 		Im(ilo1-1:ihi1+1,ilo2-1:ihi2+1,ilo3-1:ihi3+1,1:QVAR,1:3)

    real(rt), intent(in   ) :: 		dx,dy,dz,dt,a_old

	real(rt) 				:: 		dQL(7), dQR(7), dW(7), leig(7,7), reig(7,7), lam(7), summ(7)
	real(rt)				:: 		temp(s_l1:s_h1,s_l2:s_h2,s_l3:s_h3,8), smhd(7)
    real(rt)				:: 		dt_over_a
    integer          		:: 		ii,ibx,iby,ibz, i , j, k

    dt_over_a = dt / a_old
	Ip = 0.d0
	Im = 0.d0
	temp(:,:,:,1:5) = s(:,:,:,QRHO:QPRES) !Gas vars
	temp(:,:,:,6:8) = s(:,:,:,QMAGX:QMAGZ) !Mag vars
	ibx = 6
	iby = 7
	ibz = 8
	!PLM
	do k = ilo3, ihi3
		do j = ilo2, ihi2
			do i = ilo1, ihi1
	!============================================ X Direction ==============================================
				summ = 0.d0
				smhd = 0.d0
				dQL = 0.d0
				dQR = 0.d0
				dW = 0.d0
				!Skip Bx
				dQL(1:5) = 	temp(i,j,k,1:ibx-1) - temp(i-1,j,k,1:ibx-1) !gas
				dQL(6:7) = 	temp(i,j,k,ibx+1:8) - temp(i-1,j,k,ibx+1:8)	!mag			
				dQR(1:5) = 	temp(i+1,j,k,1:ibx-1) - temp(i,j,k,1:ibx-1)
				dQR(6:7) = 	temp(i+1,j,k,ibx+1:8) - temp(i,j,k,ibx+1:8)				
				do ii = 1,7
					call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
				enddo
				call evals(lam, dW, 1) !!X dir eigenvalues
				call lvecx(leig,dW)    !! left eigenvectors
				call rvecx(reig,dW)    !!right eigenvectors
				!!Using HLLD so sum over all eigenvalues
				do ii = 1,7
					summ(:) = summ(:) + lam(ii)*dot_product(leig(ii,:),dW)*reig(:,ii)
				enddo
	!MHD Source Terms 
				smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
				smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
				smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
				smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
				smhd(6) = temp(i,j,k,3)
				smhd(7) = temp(i,j,k,4)
				smhd 	= smhd*(bx(i,j,k) - bx(i-1,j,k))/dx !cross-talk of normal magnetic field direction
	!Interpolate
				Ip(i,j,k,QRHO:QPRES,1) 	 = temp(i,j,k,1:ibx-1) +0.5d0*(1.d0 - dt_over_a/dx*summ(1:5)) + 0.5d0*dt_over_a*smhd(1:5)
				Ip(i,j,k,QMAGX,1) 		 = temp(i+1,j,k,ibx) !! Bx stuff
				Ip(i,j,k,QMAGY:QMAGZ,1)  = temp(i,j,k,iby:ibz) +0.5d0*(1.d0 - dt_over_a/dx*summ(6:7)) + 0.5d0*dt_over_a*smhd(6:7)
				Im(i,j,k,QRHO:QPRES,1)	 = temp(i,j,k,1:ibx-1) -0.5d0*(1.d0 - dt_over_a/dx*summ(1:5)) - 0.5d0*dt_over_a*smhd(1:5)
				Im(i,j,k,QMAGX,1)		 = temp(i-1,j,k,ibx) !! Bx stuff
				Im(i,j,k,QMAGY:QMAGZ,1)  = temp(i,j,k,iby:ibz) -0.5d0*(1.d0 - dt_over_a/dx*summ(6:7)) - 0.5d0*dt_over_a*smhd(6:7)
				
	!========================================= Y Direction ================================================				
				summ = 0.d0
				smhd = 0.d0
				dQL = 0.d0
				dQR = 0.d0
				dW = 0.d0
				!Skip By
				dQL(1:6) = 	temp(i,j,k,1:ibx) - temp(i,j-1,k,1:ibx) !gas + bx
				dQL(7) = 	temp(i,j,k,8) - temp(i,j-1,k,iby+1)		!bz			
				dQR(1:6) = 	temp(i,j+1,k,1:ibx) - temp(i,j,k,1:ibx)
				dQR(7) = 	temp(i,j+1,k,ibz) - temp(i,j,k,ibz)				
				do ii = 1,7
					call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
				enddo
				call evals(lam, dW, 2) !!Y dir eigenvalues
				call lvecy(leig,dW)    !!left eigenvectors
				call rvecy(reig,dW)    !!right eigenvectors
				!!Using HLLD so sum over all eigenvalues
				do ii = 1,7
					summ(:) = summ(:) + lam(ii)*dot_product(leig(ii,:),dW)*reig(:,ii)
				enddo
	!MHD Source Terms 
				smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
				smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
				smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
				smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
				smhd(6) = temp(i,j,k,2)
				smhd(7) = temp(i,j,k,4)
				smhd 	= smhd*(by(i,j,k) - by(i-1,j,k))/dy !cross-talk of normal magnetic field direction
	!Interpolate
				Ip(i,j,k,QRHO:QPRES,2) 	= temp(i,j,k,1:ibx-1) +0.5d0*(1.d0 - dt_over_a/dy*summ(1:5)) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
				Ip(i,j,k,QMAGX,2) 		= temp(i,j,k,ibx) + 0.5d0*(1.d0 - dt_over_a/dy*summ(6)) + 0.5d0*dt_over_a*smhd(6)
				Ip(i,j,k,QMAGY,2) 		= temp(i,j+1,k,iby) !! By stuff
				Ip(i,j,k,QMAGZ,2)  		= temp(i,j,k,ibz) + 0.5d0*(1.d0 - dt_over_a/dy*summ(7)) + 0.5d0*dt_over_a*smhd(7)
				Im(i,j,k,QRHO:QPRES,2)	= temp(i,j,k,1:ibx-1) - 0.5d0*(1.d0 - dt_over_a/dy*summ(1:5)) - 0.5d0*dt_over_a*smhd(1:5) !!GAS
				Im(i,j,k,QMAGX,2) 		= temp(i,j,k,ibx) - 0.5d0*(1.d0 - dt_over_a/dy*summ(6)) + 0.5d0*dt_over_a*smhd(6)
				Im(i,j,k,QMAGY,2)		= temp(i,j-1,k,iby) !! By stuff
				Im(i,j,k,QMAGZ,2) 		= temp(i,j,k,ibz) -0.5d0*(1.d0 - dt_over_a/dy*summ(7)) - 0.5d0*dt_over_a*smhd(7)
				
	!========================================= Z Direction ================================================				
				summ = 0.d0
				smhd = 0.d0
				dQL = 0.d0
				dQR = 0.d0
				dW = 0.d0
				!Skip Bz
				dQL(1:7) = 	temp(i,j,k,1:iby) - temp(i,j,k-1,1:iby) 
				dQR(1:7) = 	temp(i,j,k+1,1:iby) - temp(i,j,k,1:iby)
				do ii = 1,7
					call vanleer(dW(ii),dQL(ii),dQR(ii)) !!slope limiting
				enddo
				call evals(lam, dW, 3) !!Z dir eigenvalues
				call lvecz(leig,dW)    !!left eigenvectors
				call rvecz(reig,dW)    !!right eigenvectors
				!!Using HLLD so sum over all eigenvalues
				do ii = 1,7
					summ(:) = summ(:) + lam(ii)*dot_product(leig(ii,:),dW(:))*reig(:,ii)
				enddo
	!MHD Source Terms 
				smhd(2) = temp(i,j,k,ibx)/temp(i,j,k,1)
				smhd(3) = temp(i,j,k,iby)/temp(i,j,k,1)
				smhd(4) = temp(i,j,k,ibz)/temp(i,j,k,1)
				smhd(5) = temp(i,j,k,ibx)*temp(i,j,k,2) + temp(i,j,k,iby)*temp(i,j,k,3) + temp(i,j,k,ibz)*temp(i,j,k,4)
				smhd(6) = temp(i,j,k,2)
				smhd(7) = temp(i,j,k,3)
				smhd 	= smhd*(bz(i,j,k) - bz(i-1,j,k))/dz !cross-talk of normal magnetic field direction
	!Interpolate
				Ip(i,j,k,QRHO:QPRES,3) 	= temp(i,j,k,1:ibx-1) + 0.5d0*(1.d0 - dt_over_a/dz*summ(1:5)) + 0.5d0*dt_over_a*smhd(1:5) !!GAS
				Ip(i,j,k,QMAGX:QMAGY,3)	= temp(i,j,k,ibx:iby) + 0.5d0*(1.d0 - dt_over_a/dz*summ(6:7)) + 0.5d0*dt_over_a*smhd(6:7)
				Ip(i,j,k,QMAGZ,3) 		= temp(i,j,k+1,ibz) !! Bz stuff
				Im(i,j,k,QRHO:QPRES,3)	= temp(i,j,k,1:ibx-1) - 0.5d0*(1.d0 - dt_over_a/dz*summ(1:5)) - 0.5d0*dt_over_a*smhd(1:5) !!GAS
				Im(i,j,k,QMAGX:QMAGY,3) = temp(i,j,k,ibx:iby) - 0.5d0*(1.d0 - dt_over_a/dz*summ(6:7)) + 0.5d0*dt_over_a*smhd(6:7)
				Im(i,j,k,QMAGZ,3)		= temp(i,j,k-1,ibz) !! Bz stuff
			enddo
		enddo
	enddo

!Need to add source terms, heating cooling, gravity, etc.

	end subroutine plm


!========================================= VanLeer TVD slope limiter =======================================
	subroutine vanleer(dW, WR, WL) 
	
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	real(rt), intent(in )	::  WR, WL
	real(rt), intent(out)	::  dW
	dW = 0.0d0 	
	
	if( WR*WL .gt. 0.0d0 ) then 
	dW = 2.0d0*WR*WL/(WR + WL)
	endif 

	end subroutine


!=========================================== Evals =========================================================

	subroutine evals(lam, Q, dir)

	use amrex_fort_module, only : rt => amrex_real

	implicit none
	
	real(rt), intent(in)	:: Q(QVAR)
	real(rt), intent(out)	:: lam(7) !7 waves
	integer, intent(in)		:: dir !Choose direction, 1 for x, 2 for y, 3 for z

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cfy, cfz, cax, cay, caz, csx, csy, csz, ca, as

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca  = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))

	if(dir.eq.1) then	
		!Ax eigenvalues
		lam(1) = Q(QU) - cfx
		lam(2) = Q(QU) - cax
		lam(3) = Q(QU) - csx
		lam(4) = Q(QU)
		lam(5) = Q(QU) + csx 
		lam(6) = Q(QU) + cax
		lam(7) = Q(QU) + cfx
	elseif(dir.eq.2) then
		!Ay eigenvalues
		lam(1) = Q(QV) - cfy
		lam(2) = Q(QV) - cay
		lam(3) = Q(QV) - csy
		lam(4) = Q(QV)
		lam(5) = Q(QV) + csy 
		lam(6) = Q(QV) + cay
		lam(7) = Q(QV) + cfy
	else
		!Az eigenvalues
		lam(1) = Q(QW) - cfz
		lam(2) = Q(QW) - caz
		lam(3) = Q(QW) - csz
		lam(4) = Q(QW)
		lam(5) = Q(QW) + csz 
		lam(6) = Q(QW) + caz
		lam(7) = Q(QW) + cfz
	endif
	end subroutine evals
	
!====================================== Left Eigenvectors ===============================================

!x direction
	subroutine lvecx(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ax
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cax, csx, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	!useful constants
	alf = (as - csx)/(cfx - csx)
	als = (cfx - as)/(cfx - csx)
	bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	cff = cfx*alf
	css = csx*als
	S = sign(1.0d0, Q(QMAGX))
	Qf = cfx*alf*S
	Qs = csx*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
	leig(1,:) = (/0.d0, -N*Cff	, N*Qs*bety		, N*Qs*betz		, N*alf/Q(QRHO)	, N*AAs*bety/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) !u - cf
	leig(2,:) = (/0.d0,  0.d0	, -0.5d0*betz	, 0.5d0*bety	, 0.d0			, -0.5d0*betz/(sqrt(Q(QRHO)))	, 0.5d0*bety*S/(sqrt(Q(QRHO)))	/) !u - cAx
	leig(3,:) = (/0.d0, -N*Css	, -N*Qf*bety	, -N*Qf*betz	, N*als/Q(QRHO)	, -N*AAf*bety/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) !u - cs
	leig(4,:) = (/1.d0,  0.d0	,  0.d0			, 0.d0			, -1.d0/as		, 0.d0							, 0.d0							/) !u 
	leig(5,:) = (/0.d0,  N*Css	, N*Qf*bety		, N*Qf*betz		, N*als/Q(QRHO)	, -N*AAf*bety/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) !u + cs
	leig(6,:) = (/0.d0,  0.d0	, 0.5d0*betz	, -0.5d0*bety	, 0.d0			, -0.5d0*betz*S/(sqrt(Q(QRHO)))	, 0.5d0*bety*S/(sqrt(Q(QRHO)))	/) !u + cAx
	leig(7,:) = (/0.d0, N*Cff	, -N*Qs*bety	, -N*Qs*betz	, N*alf/Q(QRHO)	, N*AAs*bety/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) !u + cf
	
	
	end subroutine lvecx

!y direction
	subroutine lvecy(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Ay
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfy, cay, csy, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	!Sloooooooooow
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	!Fassssst
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	!useful constants
	alf = (as - csy)/(cfy - csy)
	als = (cfy - as)/(cfy - csy)
	betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	cff = cfy*alf
	css = csy*als
	S = sign(1.0d0, Q(QMAGY))
	Qf = cfy*alf*S
	Qs = csy*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
!Need to double check the rows
	leig(1,:) = (/0.d0, -N*Cff , N*Qs*betx		, N*Qs*betz		, N*alf/Q(QRHO)	, N*AAs*betx/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) ! v - cf
	leig(2,:) = (/0.d0,  0.d0  , -0.5d0*betz	, 0.5d0*betx	, 0.d0			, -0.5d0*betz/(sqrt(Q(QRHO)))	, 0.5d0*betx*S/(sqrt(Q(QRHO)))	/) ! v - cAy
	leig(3,:) = (/0.d0, -N*Css , -N*Qf*betx		, -N*Qf*betz	, N*als/Q(QRHO)	, -N*AAf*betx/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) ! v - cs
	leig(4,:) = (/1.d0,  0.d0  ,  0.d0			, 0.d0			, -1.d0/as		, 0.d0							, 0.d0							/) ! v 
	leig(5,:) = (/0.d0,  N*Css , N*Qf*betx		, N*Qf*betz		, N*als/Q(QRHO)	, -N*AAf*betx/Q(QRHO)			, -N*AAf*betz/Q(QRHO)			/) ! v + cs
	leig(6,:) = (/0.d0,  0.d0  , 0.5d0*betz		, -0.5d0*betx	, 0.d0			, -0.5d0*betz*S/(sqrt(Q(QRHO)))	, 0.5d0*betx*S/(sqrt(Q(QRHO)))	/) ! v + cAy
	leig(7,:) = (/0.d0, N*Cff  , -N*Qs*betx		, -N*Qs*betz	, N*alf/Q(QRHO)	, N*AAs*betx/Q(QRHO)			, N*AAs*betz/Q(QRHO)			/) ! v + cf
	
	
	end subroutine lvecy

!z direction
	subroutine lvecz(leig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes Leig, where the rows are the left eigenvectors of the characteristic matrix Az
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::leig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfz, caz, csz, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
	!useful constants
	alf = (as - csz)/(cfz - csz)
	als = (cfz - as)/(cfz - csz)
	betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	cff = cfz*alf
	css = csz*als
	S = sign(1.0d0, Q(QMAGZ))
	Qf = cfz*alf*S
	Qs = csz*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	N = 0.5d0/as
	
!Need to double check the order
	leig(1,:) = (/0.d0, -N*Cff, N*Qs*betx	, N*Qs*bety		, N*alf/Q(QRHO) , N*AAs*betx/Q(QRHO)		    , N*AAs*bety/Q(QRHO)			/) !w - cf
	leig(2,:) = (/0.d0,  0.d0 , -0.5d0*bety , 0.5d0*betx	, 0.d0			, -0.5d0*bety/(sqrt(Q(QRHO)))   , 0.5d0*betx*S/(sqrt(Q(QRHO)))	/) !w - cAz
	leig(3,:) = (/0.d0, -N*Css, -N*Qf*betx  , -N*Qf*bety	, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO)		    , -N*AAf*bety/Q(QRHO)			/) !w - cs
	leig(4,:) = (/1.d0,  0.d0 ,  0.d0	    , 0.d0			, -1.d0/as		, 0.d0						    , 0.d0							/) !w
 	leig(5,:) = (/0.d0,  N*Css, N*Qf*betx   , N*Qf*bety		, N*als/Q(QRHO) , -N*AAf*betx/Q(QRHO)		    , -N*AAf*bety/Q(QRHO)			/) !w + cs
	leig(6,:) = (/0.d0,  0.d0 , 0.5d0*bety  , -0.5d0*betx	, 0.d0			, -0.5d0*betz*S/(sqrt(Q(QRHO))) , 0.5d0*betx*S/(sqrt(Q(QRHO)))  /) !w + cAz
	leig(7,:) = (/0.d0, N*Cff, -N*Qs*betx   , -N*Qs*bety	, N*alf/Q(QRHO) , N*AAs*bety/Q(QRHO)			, N*AAs*bety/Q(QRHO)			/) !w + cf
	end subroutine lvecz

!====================================== Right Eigenvectors ===============================================
!x direction
	subroutine rvecx(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ax
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfx, cax, csx, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, bety, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cax = (Q(QMAGX)**2)/Q(QRHO)
	!Sloooooooooow
	csx = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cax))
	!Fassssst
	cfx = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cax))
	!useful constants
	alf = (as - csx)/(cfx - csx)
	als = (cfx - as)/(cfx - csx)
	bety = Q(QMAGY)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	betz = Q(QMAGZ)/(sqrt(Q(QMAGY)**2 + Q(QMAGZ)**2))
	cff = cfx*alf
	css = csx*als
	S = sign(1.0d0, Q(QMAGX))
	Qf = cfx*alf*S
	Qs = csx*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	
				!   u - cf 				u - Cax 				u - cs			u 		u + cs			u + Cax					u + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*bety			, -betz					, -Qf*bety		, 0.d0  , Qf*bety		, betz 					, -Qs*bety		/)
	reig(4,:) = (/	Qs*betz			, bety 					, -Qf*betz		, 0.d0  , Qf*betz		, -bety 				, -Qs*betz		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*bety		, -betz*S*sqrt(Q(QRHO))	, -AAf*bety		, 0.d0  , -AAf*bety		, -betz*S*sqrt(Q(QRHO)) , AAs*bety		/)
	reig(7,:) = (/	AAs*betz		, bety*S*sqrt(Q(QRHO))	, -AAf*betz		, 0.d0	, -AAf*betz		, bety*S*sqrt(Q(QRHO))	, AAs*betz		/)
	
	
	end subroutine rvecx

!y direction
	subroutine rvecy(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Ay
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfy, cay, csy, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, betz

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	cay = (Q(QMAGY)**2)/Q(QRHO)
	!Sloooooooooow
	csy = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*cay))
	!Fassssst
	cfy = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*cay))
	!useful constants
	alf = (as - csy)/(cfy - csy)
	als = (cfy - as)/(cfy - csy)
	betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	betz = Q(QMAGZ)/(sqrt(Q(QMAGX)**2 + Q(QMAGZ)**2))
	cff = cfy*alf
	css = csy*als
	S = sign(1.0d0, Q(QMAGY))
	Qf = cfy*alf*S
	Qs = csy*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
	
				!   v - cf 				v - Cay 				v - cs			v 		v + cs			v + Cay					v + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*betx			, -betz					, -Qf*betx		, 0.d0  , Qf*betx		, betz 					, -Qs*betx		/)
	reig(4,:) = (/	Qs*betz			, betx 					, -Qf*betz		, 0.d0  , Qf*betz		, -betx 				, -Qs*betz		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*betx		, -betz*S*sqrt(Q(QRHO))	, -AAf*betx		, 0.d0  , -AAf*betx		, -betz*S*sqrt(Q(QRHO)) , AAs*betx		/)
	reig(7,:) = (/	AAs*betz		, betx*S*sqrt(Q(QRHO))	, -AAf*betz		, 0.d0	, -AAf*betz		, betx*S*sqrt(Q(QRHO))	, AAs*betz		/)
	
	
	end subroutine rvecy

!z direction
	subroutine rvecz(reig, Q) 
	use amrex_fort_module, only : rt => amrex_real
	
	implicit none
	
	!returnes reig, where the cols are the right eigenvectors of the characteristic matrix Az
	real(rt), intent(in	)	::Q(QVAR)
	real(rt), intent(out)	::reig(7,7)

	!The characteristic speeds of the system 
	real(rt)				:: cfz, caz, csz, ca, as, S, N
	real(rt)				:: cff, css, Qf, Qs, AAf, AAs, alf, als, betx, bety

	!Speeeeeeeedssssss
	as = gamma_const * Q(QPRES)/Q(QRHO)
	!Alfven
	ca = (Q(QMAGX)**2 + Q(QMAGY)**2 + Q(QMAGZ)**2)/Q(QRHO)
	caz = (Q(QMAGZ)**2)/Q(QRHO)
	!Sloooooooooow
	csz = 0.5d0*((as + ca) - sqrt((as + ca)**2 - 4.0d0*as*caz))
	!Fassssst
	cfz = 0.5d0*((as + ca) + sqrt((as + ca)**2 - 4.0d0*as*caz))
	!useful constants
	alf = (as - csz)/(cfz - csz)
	als = (cfz - as)/(cfz - csz)
	betx = Q(QMAGX)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	bety = Q(QMAGY)/(sqrt(Q(QMAGX)**2 + Q(QMAGY)**2))
	cff = cfz*alf
	css = csz*als
	S = sign(1.0d0, Q(QMAGZ))
	Qf = cfz*alf*S
	Qs = csz*als*S
	AAf = sqrt(as)*alf*sqrt(Q(QRHO))
	AAs = sqrt(as)*als*sqrt(Q(QRHO))
				!   w - cf 				w - Caz 				w - cs			w 		w + cs			w + Caz					w + cf
	reig(1,:) = (/	Q(QRHO)*alf		, 0.d0					, Q(QRHO)*als	, 1.d0	, Q(QRHO)*als	, 0.d0					, Q(QRHO)*alf	/)
	reig(2,:) = (/	-cff    	    , 0.d0  				, -css			, 0.d0  , css			, 0.d0  				, 	cff			/)
	reig(3,:) = (/	Qs*betx			, -bety					, -Qf*betx		, 0.d0  , Qf*betx		, bety 					, -Qs*betx		/)
	reig(4,:) = (/	Qs*bety			, betx 					, -Qf*bety		, 0.d0  , Qf*bety		, -betx 				, -Qs*bety		/)
	reig(5,:) = (/	Q(QRHO)*as*alf	, 0.d0	 				, Q(QRHO)*as*als, 0.d0  , Q(QRHO)*as*als, 0.d0  				, Q(QRHO)*as*alf/)
	reig(6,:) = (/	AAs*betx		, -bety*S*sqrt(Q(QRHO))	, -AAf*betx		, 0.d0  , -AAf*betx		, -bety*S*sqrt(Q(QRHO)) , AAs*betx		/)
	reig(7,:) = (/	AAs*bety		, betx*S*sqrt(Q(QRHO))	, -AAf*bety		, 0.d0	, -AAf*bety		, betx*S*sqrt(Q(QRHO))	, AAs*bety		/)
	
	
	end subroutine rvecz
end module mhd_plm_module
