module hlld_solver

implicit none

private hlldx, hlldy, hlldz, primtofluxx, primtofluxy, primtofluxz
public hlld

contains

subroutine hlld(qm,qp,qRd_l1,qRd_l2,qRd_l3,qRd_h1,qRd_h2,qRd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
				dir)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module!, only: QVAR
implicit none 

	integer, intent(in)   :: qRd_l1,qRd_l2,qRd_l3,qRd_h1,qRd_h2,qRd_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	integer, intent(in)   :: dir

	real(rt), intent(in)  :: qm(qRd_l1:qRd_h1,qRd_l2:qRd_h2,qRd_l3:qRd_h3,QVAR)
	real(rt), intent(in)  :: qp(qRd_l1:qRd_h1,qRd_l2:qRd_h2,qRd_l3:qRd_h3,QVAR)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR)
	
	integer				  :: i, j, k

 !flx = 0.d0 
 if(dir.eq.1) then
	do k = flx_l3, flx_h3
		do j = flx_l2, flx_h2
			do i = flx_l1-1, flx_h1-1
				call hlldx(i,j,k,qp(i,j,k,:),qm(i+1,j,k,:),flx(i+1,j,k,:))
			enddo
		enddo
	enddo

 elseif(dir.eq.2) then

	do k = flx_l3, flx_h3
		do j = flx_l2-1, flx_h2-1
			do i = flx_l1, flx_h1
				call hlldy(i,j,k,qp(i,j,k,:),qm(i,j+1,k,:),flx(i,j+1,k,:))
			enddo
		enddo
	enddo
 else 
	do k = flx_l3-1, flx_h3-1
		do j = flx_l2, flx_h2
			do i = flx_l1, flx_h1
				call hlldz(i,j,k,qp(i,j,k,:),qm(i,j,k+1,:),flx(i,j,k+1,:))
			enddo
		enddo
	enddo
 endif
end subroutine hlld

!================================================= X Direction =======================================================
subroutine hlldx(i,j,k,qL,qR,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bx respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	integer , intent(in)  :: i,j,k
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, caxL
	real(rt) 			  :: caR, caxR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)              :: uL(QVAR), uR(QVAR)	
	real(rt)			  :: UsL(QVAR), FsL(QVAR)
	real(rt)			  :: UsR(QVAR), FsR(QVAR)
	real(rt)			  :: UssL(QVAR), FssL(QVAR)
	real(rt)			  :: UssR(QVAR), FssR(QVAR)
	character(len=10)	  :: choice
!Riemann solve
	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	UsL = 0.d0
	UsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	UssL = 0.d0
	UssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0
	
    call PToC(qL,uL)
    call PToC(qR,uR)	
	call primtofluxx(qL, FL)
	call primtofluxx(qR, FR)	
	
	eL   = (qL(QPRES) -0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)
	eR   = (qR(QPRES) -0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)
	asL  = gamma_const * (qL(QPRES) - 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/qL(QRHO)
	asR  = gamma_const * (qR(QPRES) - 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/qR(QRHO)
	caL  = (qL(QMAGX)**2 + qL(QMAGY)**2 + qL(QMAGZ)**2)/qL(QRHO) !Magnetic Speeds
	caR  = (qR(QMAGX)**2 + qR(QMAGY)**2 + qR(QMAGZ)**2)/qR(QRHO)
	caxL = (qL(QMAGX)**2)/qL(QRHO)
	caxR = (qR(QMAGX)**2)/qR(QRHO)
	!Catch the fastest waves, brah
	cfL  = sqrt(0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*caxL)))
	cfR  = sqrt(0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*caxR)))

	!Riemann Speeds
	sL   = min(qL(QU) - cfL, qR(QU) - cfR)
	sR 	 = max(qL(QU) + cfL,qR(QU) + cfR)
	sM   = ((sR - qR(QU))*qR(QRHO)*qR(QU) - (sL - qL(QU))*qL(QRHO)*qL(QU) - qR(QPRES) + qL(QPRES))/((sR - qR(QU))*qR(QRHO) - (sL - qL(QU))*qL(QRHO))

	!Pressures in the Riemann Fan
	ptL  = qL(QPRES)
	ptR  = qR(QPRES)
	pst  = (sR - qR(QU))*qR(QRHO)*ptL - (sL - qL(QU))*qL(QRHO)*ptR + qL(QRHO)*qR(QRHO)*(sR - qR(QU))*(sL - qL(QU))*(qR(QU) - qL(QU))
	pst  = pst/((sR - qR(QU))*qR(QRHO) - (sL - qL(QU))*qL(QRHO))

	!------------------------------------------- * states-------------------------------------------------------------------------
	!density
	UsL(QRHO) = qL(QRHO)*((sL - qL(QU))/(sL - sM))
	UsR(QRHO) = qR(QRHO)*((sR - qR(QU))/(sR - sM))
	!velocities
	!X dir
	UsL(QU)    = sM
	UsR(QU)    = sM
	!Y dir
	UsL(QV)    = qL(QV) - qL(QMAGX)*qL(QMAGY)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2))
	UsR(QV)    = qR(QV) - qR(QMAGX)*qR(QMAGY)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qL(QMAGX)**2))
	!Z dir
	UsL(QW)    = qL(QW) - qL(QMAGX)*qL(QMAGZ)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2))
	UsR(QW)    = qR(QW) - qR(QMAGX)*qR(QMAGZ)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qL(QMAGX)**2))
		
	UsL(QU:QW) = UsL(QU:QW)*UsL(QRHO)
	UsR(QU:QW) = UsR(QU:QW)*UsR(QRHO)

	!Magnetic Fields
	!X dir
	UsL(QMAGX) = qL(QMAGX)
	UsR(QMAGX) = qL(QMAGX) 

	!Y dir
	UsL(QMAGY) = qL(QMAGY)*(qL(QRHO)*(sL - qL(QU))**2 - qL(QMAGX)**2)/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2)
	UsR(QMAGY) = qR(QMAGY)*(qR(QRHO)*(sR - qR(QU))**2 - qR(QMAGX)**2)/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qR(QMAGX)**2)

	!Z dir
	UsL(QMAGZ) = qL(QMAGZ)*(qL(QRHO)*(sL - qL(QU))**2 - qL(QMAGX)**2)/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2)
	UsR(QMAGZ) = qR(QMAGZ)*(qR(QRHO)*(sR - qR(QU))**2 - qR(QMAGX)**2)/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qR(QMAGX)**2)
	
	!Energy *Stored in Pressure slot
	UsL(QPRES) = (sL - qL(QU))*eL - ptL*qL(QU) + pst*sM + qL(QMAGX)*(dot_product(qL(QU:QW),qL(QMAGX:QMAGZ)) &
				  - dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)))
	UsL(QPRES) = UsL(QPRES)/(sL - sM)
	UsR(QPRES) = (sR - qR(QU))*eR - ptR*qR(QU) + pst*sM + qR(QMAGX)*(dot_product(qR(QU:QW),qR(QMAGX:QMAGZ)) &
				  - dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)))
	UsR(QPRES) = UsR(QPRES)/(sR - sM)
	!speeds
	ssL = sM - abs(qL(QMAGX))/sqrt(UsL(QRHO))
	ssR = sM + abs(qR(QMAGX))/sqrt(UsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	UssL(QRHO)  = UsL(QRHO)
	UssR(QRHO)  = UsR(QRHO)
	!u
	UssL(QU)    = sM
	UssR(QU)    = sM

	!v
	UssL(QV)    = (sqrt(UsL(QRHO))*UsL(QV)/UsL(QRHO) + &
                       sqrt(UsR(QRHO))*UsR(QV)/UsR(QRHO) + &
		    + (UsR(QMAGY) - UsL(QMAGY)) * sign(1.d0,qL(QMAGX)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QV)    = UssL(QV)

	!w
	UssL(QW)    = (sqrt(UsL(QRHO))*UsL(QW)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QW)/UsR(QRHO) &
		    + (UsR(QMAGZ) - UsL(QMAGZ)) * sign(1.d0,qL(QMAGX)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QW)    = UssL(QW)
	UssL(QU:QW) = UssL(QU:QW)*UssL(QRHO)
	UssR(QU:QW) = UssR(QU:QW)*UssR(QRHO)

	!Bx
	UssL(QMAGX) = UsL(QMAGX)
	UssR(QMAGX) = UsR(QMAGX)

	!By
	UssL(QMAGY) = (sqrt(UsL(QRHO))*UsR(QMAGY) + &
                       sqrt(UsR(QRHO))*UsL(QMAGY) + &
                       sqrt(UsL(QRHO)*UsR(QRHO)) * (UsR(QV)/UsR(QRHO) - UsL(QV)/UsL(QRHO)) * sign(1.d0,UsR(QMAGX)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGY) = UssL(QMAGY)

	!Bz
	UssL(QMAGZ) = (sqrt(UsL(QRHO))*UsR(QMAGZ) + &
                       sqrt(UsR(QRHO))*UsL(QMAGZ) + &
                       sqrt(UsL(QRHO)*UsR(QRHO)) * (UsR(QW)/UsR(QRHO) - UsL(QW)/UsL(QRHO)) * sign(1.d0,UsR(QMAGX)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGZ) = UssL(QMAGZ)

	!Energy *Stored in Pressure slot

	UssL(QPRES) = UsL(QPRES) - sqrt(UsL(QRHO))*(dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)) &
				- dot_product(UssL(QU:QW)/UssL(QRHO),UssL(QMAGX:QMAGZ))) * sign(1.d0, UsL(QMAGX))

	UssR(QPRES) = UsR(QPRES) + sqrt(UsR(QRHO))*(dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)) &
				- dot_product(UssR(QU:QW)/UssR(QRHO),UssR(QMAGX:QMAGZ))) * sign(1.d0, UsR(QMAGX))

	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(UsL - uL)
	FssL = FsL + ssL*(UssL - UsL)
	FsR  = FR + sR*(UsR - uR)
	FssR = FsR + ssR*(UssR - UsR)
	!Solve the RP
	if(sL .gt. 0.d0) then
	flx = FL
	choice = "FL"
	elseif(sL .le. 0.d0 .and. ssL .gt. 0.d0) then
	flx = FsL
	choice = "FsL"
	elseif(ssl .le. 0.d0 .and. sM .gt. 0.d0) then
	flx = FssL
	choice = "FssL"
	elseif(sM .le. 0.d0 .and. ssR .gt. 0.d0) then
	flx = FssR
	choice = "FssR"
	elseif(ssR .le. 0.d0 .and. sR .gt. 0.d0) then
	flx = FsR
	choice = "FsR"
	else 
	flx = FR
	choice = "FR"
	endif

        if (abs(flx(QU)) .gt. 0.) print *,'BAD FLX(QU) IN HLLDX ',i,j,k,flx(QU) 
        if (abs(flx(QV)) .gt. 0.) print *,'BAD FLX(QV) IN HLLDX ',i,j,k,flx(QV) 
        if (abs(flx(QW)) .gt. 0.) print *,'BAD FLX(QW) IN HLLDX ',i,j,k,flx(QW) 

!		if (i.ge.15 .and. i.le.18 .and. j.eq.-3 .and. k.eq.-2) then
!			write(*,*) "Flux at i = ", i
!			write(*,*) "Flux = ", choice
!			write(*,*) "FL = ", FL(QMAGY)
!			write(*,*) " sL = ",  sL
!			write(*,*) " sR = ",  sR
!			write(*,*) "QLmag = ", qL(QMAGX)
!			write(*,*) " QR(QU) = ",  qR(QU)
!			write(*,*) " QL(QRHO) = ",  qL(QRHO)
!			write(*,*) " QR(QRHO) = ",  qR(QRHO)
!			write(*,*) " QL(MAGY) = ",  qL(QMAGY)
!			write(*,*) " QR(MAGY) = ",  qR(QMAGY)
!			write(*,*) " QL(PRES) = ",  qL(QPRES)
!			write(*,*) " QR(PRES) = ",  qR(QPRES)
!			write(*,*) " ssL = ", ssL
!			write(*,*) " uL = ",  uL(QMAGY) 
!			write(*,*) " sR = ",  sR
!			write(*,*) " sM = ",  sM
!			write(*,*) "UssL = ", UssL(QRHO)
!			write(*,*) "UssR = ", UssR(QRHO)
!		endif

end subroutine hlldx

!============================================================= Y Direction =================================================================

subroutine hlldy(i,j,k,qL,qR,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/By respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	integer , intent(in)  :: i,j,k
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)	  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, cayL
	real(rt) 	  :: caR, cayR, asL, asR, ptL, ptR, eL, eR
	real(rt)	  :: FL(QVAR), FR(QVAR)
        real(rt)          :: uL(QVAR), uR(QVAR)
	real(rt)	  :: UsL(QVAR), FsL(QVAR)
	real(rt)	  :: UsR(QVAR), FsR(QVAR)
	real(rt)	  :: UssL(QVAR), FssL(QVAR)
	real(rt)	  :: UssR(QVAR), FssR(QVAR)

	character(len=10)	  :: choice
!Riemann solve

	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	UsL = 0.d0
	UsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	UssL = 0.d0
	UssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0

       call PToC(qL,uL)
       call PToC(qR,uR)

       call primtofluxy(i,j,k,qL, FL)
       call primtofluxy(i,j,k,qR, FR)

       if (i.eq.-3.and.j.ge.28.and.k.eq.-3) then
          print *,'CALLING PRIMTO FLUX WITH J = ',j
          print *,' QL ', qL
          print *,' QR ', qR
          print *,' FL ', FL
          print *,' FR ', FR
          if (j.ge.32) stop
       end if

	eL   = (qL(QPRES) -0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)
	eR   = (qR(QPRES) -0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)
	asL  = gamma_const * (qL(QPRES) - 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/qL(QRHO)
	asR  = gamma_const * (qR(QPRES) - 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/qR(QRHO)
	caL  = (qL(QMAGX)**2 + qL(QMAGY)**2 + qL(QMAGZ)**2)/qL(QRHO) !Magnetic Speeds
	caR  = (qR(QMAGX)**2 + qR(QMAGY)**2 + qR(QMAGZ)**2)/qR(QRHO)
	cayL = (qL(QMAGY)**2)/qL(QRHO)
	cayR = (qR(QMAGY)**2)/qR(QRHO)
	!Catch the fastest waves, brah
	cfL  = sqrt(0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*cayL)))
	cfR  = sqrt(0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*cayR)))
	!Riemann Speeds
	sL   = min(qL(QV) - cfL,qR(QV) - cfR)
	sR 	 = max(qL(QV) + cfL,qR(QV) + cfR)
	ptL  = qL(QPRES)
	ptR  = qR(QPRES)
	sM   = ((sR - qR(QV))*qR(QRHO)*qR(QV) - (sL - qL(QV))*qL(QRHO)*qL(QV) - ptR + ptL)/((sR - qR(QV))*qR(QRHO) - (sL - qL(QV))*qL(QRHO))
	!Pressures in the Riemann Fan
	pst  = (sR - qR(QV))*qR(QRHO)*ptL - (sL - qL(QV))*qL(QRHO)*ptR + qL(QRHO)*qR(QRHO)*(sR - qR(QV))*(sL - qL(QV))*(qR(QV) - qL(QV))
	pst  = pst/((sR - qR(QV))*qR(QRHO) - (sL - qL(QV))*qL(QRHO))

	!------------------------------------------- * states-------------------------------------------------------------------------
	!density
	UsL(QRHO) = qL(QRHO)*((sL - qL(QV))/(sL - sM))
	UsR(QRHO) = qR(QRHO)*((sR - qR(QV))/(sR - sM))
	!velocities
	!X dir
	if(abs(qL(QMAGY)*qL(QMAGX)*(sM - qL(QV))).le. 1d-14) then
		UsL(QU) = qL(QU)
	else
		UsL(QU)    = qL(QU) - qL(QMAGY)*qL(QMAGX)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2))
	endif
	if(abs(qR(QMAGY)*qR(QMAGX)*(sM - qR(QV))).le. 1d-14) then
		UsR(QU) = qR(QU)
	else
		UsR(QU)    = qR(QU) - qR(QMAGY)*qR(QMAGX)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qR(QMAGY)**2))
	endif
	!Y dir
	UsL(QV)    = sM
	UsR(QV)    = sM
	!Z dir
	if(abs(qL(QMAGY)*qL(QMAGZ)*(sM - qL(QV))).lt. 1d-14) then
		UsL(QW)    = qL(QW)
	else
		UsL(QW)    = qL(QW) - qL(QMAGY)*qL(QMAGZ)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2))
	endif
	if(abs(qR(QMAGY)*qR(QMAGZ)*(sM - qR(QV))).lt. 1d-14) then
		UsR(QW) = qR(QW)
	else
		UsR(QW)    = qR(QW) - qR(QMAGY)*qR(QMAGZ)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qL(QMAGY)**2))
	endif
	UsL(QU:QW) = UsL(QU:QW)*UsL(QRHO)
	UsR(QU:QW) = UsR(QU:QW)*UsR(QRHO)
	
	!Magnetic Fields
	!X dir
	if(abs(qL(QMAGX)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)).lt. 1d-14) then
		UsL(QMAGX) = qL(QMAGX)
	else
		UsL(QMAGX) = qL(QMAGX)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2)
	endif
	if(abs(qR(QMAGX)*(qR(QRHO)*(sR - qR(QV))**2 - qR(QMAGY)**2)).lt.1d-14) then 
		UsR(QMAGX) = qR(QMAGX)
	else
		UsR(QMAGX) = qR(QMAGX)*(qR(QRHO)*(sR - qR(QV))**2 - qR(QMAGY)**2)/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qR(QMAGY)**2)
	endif
	!Y dir
	UsL(QMAGY) = qL(QMAGY)
	UsR(QMAGY) = qL(QMAGY) 
	!Z dir
	if(abs(qL(QMAGZ)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)).lt.1d-14) then
		UsL(QMAGZ) = qL(QMAGZ)
	else
		UsL(QMAGZ) = qL(QMAGZ)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2)
	endif
	if(abs(qR(QMAGZ)*(qR(QRHO)*(sR - qR(QV))**2 - qL(QMAGY)**2)).lt.1d-14) then 
		UsR(QMAGZ) = qR(QMAGZ)
	else
		UsR(QMAGZ) = qR(QMAGZ)*(qR(QRHO)*(sR - qR(QV))**2 - qL(QMAGY)**2)/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qL(QMAGY)**2)
	endif
	
	!Energy *Stored in Pressure slot
	UsL(QPRES) = (sL - qL(QV))*eL - ptL*qL(QV) + pst*sM + qL(QMAGY)*(dot_product(qL(QU:QW),qL(QMAGX:QMAGZ)) &
				  - dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)))
	UsL(QPRES) = UsL(QPRES)/(sL - sM)
	UsR(QPRES) = (sR - qR(QV))*eR - ptR*qR(QV) + pst*sM + qR(QMAGY)*(dot_product(qR(QU:QW),qR(QMAGX:QMAGZ)) &
				  - dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)))
	UsR(QPRES) = UsR(QPRES)/(sR - sM)

	!speeds
	ssL = sM - abs(qL(QMAGY))/sqrt(UsL(QRHO))
	ssR = sM + abs(qR(QMAGY))/sqrt(UsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	UssL(QRHO)  = UsL(QRHO)
	UssR(QRHO)  = UsR(QRHO)

	!u
	UssL(QU)    = (sqrt(UsL(QRHO))*UsL(QU)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QU)/UsR(QRHO)&
		    + (UsR(QMAGX) - UsL(QMAGX))*sign(1.d0,qL(QMAGY)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QU)    = UssL(QU)

	!v
	UssL(QV)    = sM
	UssR(QV)    = sM

	!w
	UssL(QW)    = (sqrt(UsL(QRHO))*UsL(QW)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QW)/UsR(QRHO)&
		    + (UsR(QMAGZ) - UsL(QMAGZ))*sign(1.d0,qL(QMAGY)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QW)    = UssL(QW)

	UssL(QU:QW) = UssL(QU:QW)*UssL(QRHO)
	UssR(QU:QW) = UssR(QU:QW)*UssR(QRHO)

	!Bx
	UssL(QMAGX) = (sqrt(UsL(QRHO))*UsR(QMAGX) + sqrt(UsR(QRHO))*UsL(QMAGX) + sqrt(UsL(QRHO)*UsR(QRHO)) * &
                      ( UsR(QU)/UsR(QRHO) - UsL(QU)/UsR(QRHO)) * sign(1.d0,UsR(QMAGY)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGX) = UssL(QMAGX)
	!By
	UssL(QMAGY) = UsL(QMAGY)
	UssR(QMAGY) = UsR(QMAGY)
	!Bz
	UssL(QMAGZ) = (sqrt(UsL(QRHO))*UsR(QMAGZ) + sqrt(UsR(QRHO))*UsL(QMAGZ) + &
                       sqrt(UsL(QRHO)*UsR(QRHO))*( UsR(QW)/UsR(QRHO) - UsL(QW)/UsL(QRHO) ) * sign(1.d0,UsR(QMAGY)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGZ) = UssL(QMAGZ)

	!Energy *Stored in Pressure slot

	UssL(QPRES) = UsL(QPRES) - sqrt(UsL(QRHO))*(dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)) &
				 - dot_product(UssL(QU:QW)/UssL(QRHO),UssL(QMAGX:QMAGZ)))*sign(1.d0, UsL(QMAGY))

	UssR(QPRES) = UsR(QPRES) + sqrt(UsR(QRHO))*(dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)) &
				 - dot_product(UssR(QU:QW)/UssR(QRHO),UssR(QMAGX:QMAGZ)))*sign(1.d0, UsR(QMAGY))

	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(UsL - uL)
	FssL = FsL + ssL*(UssL - UsL)
	FsR  = FR + sR*(UsR - uR)
	FssR = FsR + ssR*(UssR - UsR)
	!Solve the RP
	if(sL .gt. 0.d0) then
	flx = FL
	choice = "FL"
	elseif(sL .le. 0.d0 .and. ssL .gt. 0.d0) then
	flx = FsL
	choice = "FsL"
	elseif(ssL .le. 0.d0 .and. sM .gt. 0.d0) then
	flx = FssL
	choice = "FssL"
	elseif(sM .le. 0.d0 .and. ssR .gt. 0.d0) then
	flx = FssR
	choice = "FssR"
	elseif(ssR .le. 0.d0 .and. sR .gt. 0.d0) then
	flx = FsR
	choice = "FsR"
	else 
	flx = FR
	choice = "FR"
	endif

        if (i.eq.-3.and.j.ge.28.and.k.eq.-2) then
 			write(*,*) "Flux at j = ",j
 			write(*,*) "Flux = ", choice
 			write(*,*) "QL = ", QL(QMAGX)
 			write(*,*) "QR = ", QR(QMAGX)
 			write(*,*) "FL = ", FL(QMAGX)
 			write(*,*) "FR = ", FR(QMAGX)
 			write(*,*) "QL = ", qL(QMAGX)
 			write(*,*) "QR = ", qR(QMAGX)
 			write(*,*) "UsL = ", UsL(QMAGX)
 			write(*,*) "UsR = ", UsR(QMAGX)
 			write(*,*) "UssL = ", UssL(QMAGX)
 			write(*,*) "UssR = ", UssR(QMAGX)
 			write(*,*) " "
       endif

        if (abs(flx(QU)) .gt. 0.) print *,'BAD FLX(QU) IN HLLDY ',i,j,k,flx(QU) 
        if (abs(flx(QV)) .gt. 0.) print *,'BAD FLX(QV) IN HLLDY ',i,j,k,flx(QV) 
        if (abs(flx(QW)) .gt. 0.) print *,'BAD FLX(QW) IN HLLDY ',i,j,k,flx(QW) 
 
end subroutine hlldy

!============================================================= Z Direction =================================================================

subroutine hlldz(i,j,k,qL,qR,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bz respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	integer , intent(in)  :: i,j,k
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, cazL
	real(rt) 			  :: caR, cazR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)              :: uL(QVAR), uR(QVAR)
	real(rt)			  :: UsL(QVAR), FsL(QVAR)
	real(rt)			  :: UsR(QVAR), FsR(QVAR)
	real(rt)			  :: UssL(QVAR), FssL(QVAR)
	real(rt)			  :: UssR(QVAR), FssR(QVAR)
	character(len=10)	  :: choice

!Riemann solve
	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	UsL = 0.d0
	UsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	UssL = 0.d0
	UssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0
    
    call PToC(qL,uL)
    call PToC(qR,uR)
	call primtofluxz(qL, FL)
	call primtofluxz(qR, FR)
	
	eL   = (qL(QPRES) -0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)
	eR   = (qR(QPRES) -0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/(gamma_minus_1) + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
			+ 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)
	asL  = gamma_const * (qL(QPRES) - 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)))/qL(QRHO)
	asR  = gamma_const * (qR(QPRES) - 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)))/qR(QRHO)
	caL  = (qL(QMAGX)**2 + qL(QMAGY)**2 + qL(QMAGZ)**2)/qL(QRHO) !Magnetic Speeds
	caR  = (qR(QMAGX)**2 + qR(QMAGY)**2 + qR(QMAGZ)**2)/qR(QRHO)
	cazL = (qL(QMAGZ)**2)/qL(QRHO)
	cazR = (qR(QMAGZ)**2)/qR(QRHO)
	!Catch the fastest waves, brah
	cfL  = sqrt(0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*cazL)))
	cfR  = sqrt(0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*cazR)))
	!Riemann Speeds
	sL   = min(qL(QW) - cfL,qR(QW) - cfR)
	sR 	 = max(qL(QW) + cfL,qR(QW) + cfR)
	sM   = ((sR - qR(QW))*qR(QRHO)*qR(QW) - (sL - qL(QW))*qL(QRHO)*qL(QW) - qR(QPRES) + qL(QPRES))/((sR - qR(QW))*qR(QRHO) - (sL - qL(QW))*qL(QRHO))
	!Pressures in the Riemann Fan
	ptL  = qL(QPRES)
	ptR  = qR(QPRES)
	pst  = (sR - qR(QW))*qR(QRHO)*ptL - (sL - qL(QW))*qL(QRHO)*ptR + qL(QRHO)*qR(QRHO)*(sR - qR(QW))*(sL - qL(QW))*(qR(QW) - qL(QW))
	pst  = pst/((sR - qR(QW))*qR(QRHO) - (sL - qL(QW))*qL(QRHO))

	!------------------------------------------- * states-------------------------------------------------------------------------
	!density
	UsL(QRHO) = qL(QRHO)*((sL - qL(QW))/(sL - sM))
	UsR(QRHO) = qR(QRHO)*((sR - qR(QW))/(sR - sM))
	!velocities
	!X dir
	UsL(QU)    = qL(QU) - qL(QMAGZ)*qL(QMAGX)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2))
	UsR(QU)    = qR(QU) - qR(QMAGZ)*qR(QMAGX)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2))
	!Y dir
	UsL(QV)    = qL(QV) - qL(QMAGZ)*qL(QMAGY)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2))
	UsR(QV)    = qR(QV) - qR(QMAGZ)*qR(QMAGY)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2))
	!Z dir
	UsL(QW)    = sM
	UsR(QW)    = sM
	
	UsL(QU:QW) = UsL(QU:QW)*UsL(QRHO)
	UsR(QU:QW) = UsR(QU:QW)*UsR(QRHO)

	!Magnetic Fields
	!X dir
	UsL(QMAGX) = qL(QMAGX)*(qL(QRHO)*(sL - qL(QW))**2 - qL(QMAGZ)**2)/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2)
	UsR(QMAGX) = qR(QMAGX)*(qR(QRHO)*(sR - qR(QW))**2 - qL(QMAGZ)**2)/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2)
	!Y dir
	UsL(QMAGY) = qL(QMAGY)*(qL(QRHO)*(sL - qL(QW))**2 - qL(QMAGZ)**2)/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2)
	UsR(QMAGY) = qR(QMAGY)*(qR(QRHO)*(sR - qR(QW))**2 - qL(QMAGZ)**2)/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2)
	!Z dir
	UsL(QMAGZ) = qL(QMAGZ)
	UsR(QMAGZ) = qL(QMAGZ) 
	
	!Energy *Stored in Pressure slot
	UsL(QPRES) = (sL - qL(QW))*eL - ptL*qL(QW) + pst*sM + qL(QMAGZ)*(dot_product(qL(QU:QW),qL(QMAGX:QMAGZ))&
			   - dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)))
	UsL(QPRES) = UsL(QPRES)/(sL - sM)
	UsR(QPRES) = (sR - qR(QW))*eR - ptR*qR(QW) + pst*sM + qR(QMAGZ)*(dot_product(qR(QU:QW),qR(QMAGX:QMAGZ))&
			   - dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)))
	UsR(QPRES) = UsR(QPRES)/(sR - sM)

	!speeds
	ssL = sM - abs(qL(QMAGZ))/sqrt(UsL(QRHO))
	ssR = sM + abs(qR(QMAGZ))/sqrt(UsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	UssL(QRHO)  = UsL(QRHO)
	UssR(QRHO)  = UsR(QRHO)
	!u
	UssL(QU)    = (sqrt(UsL(QRHO))*UsL(QU)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QU)/UsR(QRHO)&
		     + (UsR(QMAGX) - UsL(QMAGX))*sign(1.d0,qL(QMAGZ)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QU)    = UssL(QU)
	!v
	UssL(QV)    = (sqrt(UsL(QRHO))*UsL(QV)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QV)/UsR(QRHO)&
		    + (UsR(QMAGY) - UsL(QMAGY))*sign(1.d0,qL(QMAGZ)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QV)    = UssL(QW)
	!w
	UssL(QW)    = sM
	UssR(QW)    = sM

	UssL(QU:QW) = UssL(QU:QW)*UssL(QRHO)
	UssR(QU:QW) = UssR(QU:QW)*UssR(QRHO)

	!Bx
	UssL(QMAGX) = (sqrt(UsL(QRHO))*UsR(QMAGX) + sqrt(UsR(QRHO))*UsL(QMAGX) + &
                       sqrt(UsL(QRHO)*UsR(QRHO))*(UsR(QU) - UsL(QU))*sign(1.d0,UsR(QMAGZ)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGX) = UssL(QMAGX)
	!By
	UssL(QMAGY) = (sqrt(UsL(QRHO))*UsR(QMAGY) + sqrt(UsR(QRHO))*UsL(QMAGY) + &
                       sqrt(UsL(QRHO)*UsR(QRHO))*(UsR(QV) - UsL(QV))*sign(1.d0,UsR(QMAGZ)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
	UssR(QMAGY) = UssL(QMAGY)
	!Bz
	UssL(QMAGZ) = UsL(QMAGZ)
	UssR(QMAGZ) = UsR(QMAGZ)

	!Energy *Stored in Pressure slot

	UssL(QPRES) = UsL(QPRES) - sqrt(UsL(QRHO))*(dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ))&
				- dot_product(UssL(QU:QW)/UssL(QRHO),UssL(QMAGX:QMAGZ)))*sign(1.d0, UsL(QMAGZ))
	UssR(QPRES) = UsR(QPRES) + sqrt(UsR(QRHO))*(dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ))&
				- dot_product(UssR(QU:QW)/UssR(QRHO),UssR(QMAGX:QMAGZ)))*sign(1.d0, UsR(QMAGZ))

	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(UsL - uL)
	FssL = FsL + ssL*(UssL - UsL)
	FsR  = FR + sR*(UsR - uR)
	FssR = FsR + ssR*(UssR - UsR)
	!Solve the RP
	if(sL .gt. 0.d0) then
	flx = FL
	choice = "FL"
	elseif(sL .le. 0.d0 .and. ssL .gt. 0.d0) then
	flx = FsL
	choice = "FsL"
	elseif(ssl .le. 0.d0 .and. sM .gt. 0.d0) then
	flx = FssL
	choice = "FssL"
	elseif(sM .le. 0.d0 .and. ssR .gt. 0.d0) then
	flx = FssR
	choice = "FssR"
	elseif(ssR .le. 0.d0 .and. sR .gt. 0.d0) then
	flx = FsR
	choice = "FsR"
	else 
	flx = FR
	choice = "FR"
	endif
!	do i = 1, QVAR
!		if(abs(flx(UMX)).gt. 1.d-14) then
!			write(*,*) "Flux is nan in", i, "component"
!			write(*,*) "Flux = ", choice
!			write(*,*) "FL = ", FL
!			write(*,*) "FR = ", FR
!			write(*,*) "QL = ", qL
!			write(*,*) "QR = ", qR
!			write(*,*) "UsL = ", UsL
!			write(*,*) "UsR = ", UsR
!			write(*,*) "UssL = ", UssL
!			write(*,*) "UssR = ", UssR
!			pause
!			return
!		endif
!	enddo
!	flx = 0.5*(FL + FR) - 0.5*sR*(UR - UL)

        if (abs(flx(QU)) .gt. 0.) print *,'BAD FLX(QU) IN HLLDZ ',i,j,k,flx(QU) 
        if (abs(flx(QV)) .gt. 0.) print *,'BAD FLX(QV) IN HLLDZ ',i,j,k,flx(QV) 
        if (abs(flx(QW)) .gt. 0.) print *,'BAD FLX(QW) IN HLLDZ ',i,j,k,flx(QW) 
end subroutine hlldz

!====================================================== Fluxes ================================================================================

!----------------------------------------- X Direction ---------------------------------------------------------

subroutine primtofluxx(Q, F)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
implicit none

	real(rt), intent(in)  :: Q(QVAR)
	real(rt), intent(out) :: F(QVAR)
	real(rt)			  :: e

	F = 0.d0
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(Q(QRHO)*gamma_minus_1)& 
			   + 0.5d0*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0/Q(QRHO)*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QU)
	F(UMX)	 = Q(QRHO)*Q(QU)**2 + Q(QPRES) - Q(QMAGX)**2
	F(UMY)   = Q(QRHO)*Q(QU)*Q(QV) - Q(QMAGX)*Q(QMAGY)
	F(UMZ)   = Q(QRHO)*Q(QU)*Q(QW) - Q(QMAGX)*Q(QMAGZ)
	F(UEDEN) = Q(QU)*(Q(QRHO)*e + Q(QPRES)) -Q(QMAGX)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = 0.d0
	F(QMAGY) = Q(QU)*Q(QMAGY) - Q(QMAGX)*Q(QV)
	F(QMAGZ) = Q(QU)*Q(QMAGZ) - Q(QMAGX)*Q(QW)

end subroutine primtofluxx

!-------------------------------------- Y Direction ------------------------------------------------------------

subroutine primtofluxy(i,j,k,Q, F)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
implicit none

	real(rt), intent(in)  :: Q(QVAR)
	real(rt), intent(out) :: F(QVAR)
	integer :: i,j,k
	real(rt)			  :: e

	F = 0.d0
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(Q(QRHO)*gamma_minus_1)& 
			   + 0.5d0*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0/Q(QRHO)*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QV)
	F(UMX)	 = Q(QRHO)*Q(QU)*Q(QV) - Q(QMAGX)*Q(QMAGY)
	F(UMY)   = Q(QRHO)*Q(QV)**2 + Q(QPRES) - Q(QMAGY)**2
	F(UMZ)   = Q(QRHO)*Q(QV)*Q(QW) - Q(QMAGY)*Q(QMAGZ)
	F(UEDEN) = Q(QV)*(Q(QRHO)*e + Q(QPRES)) -Q(QMAGY)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = Q(QV)*Q(QMAGX) - Q(QMAGY)*Q(QU)
	F(QMAGY) = 0.d0
	F(QMAGZ) = Q(QV)*Q(QMAGZ) - Q(QMAGY)*Q(QW)

        if (abs(F(QMAGX)) .gt. 0.) then
           print *,'BAD BX AT ', i,j,k
           print *,'Q(QU)   = ', Q(QU)
           print *,'Q(QV)   = ', Q(QV)
           print *,'Q(MAGX) = ', Q(QMAGX)
           print *,'Q(MAGY) = ', Q(QMAGY)
        end if

end subroutine primtofluxy

!-------------------------------------- Z Direction ------------------------------------------------------------

subroutine primtofluxz(Q, F)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
implicit none

	real(rt), intent(in)  :: Q(QVAR)
	real(rt), intent(out) :: F(QVAR)
	real(rt)			  :: e

	F = 0.d0
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(Q(QRHO)*gamma_minus_1)& 
			   + 0.5d0*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0/Q(QRHO)*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QW)
	F(UMX)	 = Q(QRHO)*Q(QW)*Q(QU) - Q(QMAGX)*Q(QMAGZ)
	F(UMY)   = Q(QRHO)*Q(QW)*Q(QV) - Q(QMAGY)*Q(QMAGZ)
	F(UMZ)   = Q(QRHO)*Q(QW)**2 + Q(QPRES) - Q(QMAGZ)**2
	F(UEDEN) = Q(QW)*(Q(QRHO)*e + Q(QPRES)) -Q(QMAGZ)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = Q(QW)*Q(QMAGX) - Q(QMAGZ)*Q(QU)
	F(QMAGY) = Q(QW)*Q(QMAGY) - Q(QMAGZ)*Q(QV)
	F(QMAGZ) = 0.d0

end subroutine primtofluxz

!================================================= Calculate the Conservative Variables ===============================================

subroutine PToC(q, u)

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none

	real(rt), intent(in)	::q(QVAR)
	real(rt), intent(out)	::u(QVAR)
			u = 0.d0

 			u(URHO)   = q(QRHO)
 			u(UMX)    = q(QRHO)*q(QU)
 			u(UMY)    = q(QRHO)*q(QV)
 			u(UMZ)    = q(QRHO)*q(QW)
			u(UEINT)  = (q(QPRES) - 0.5d0*dot_product(q(QMAGX:QMAGZ),q(QMAGX:QMAGZ)))/(gamma_minus_1)
			u(UEDEN)  = u(UEINT)  + 0.5d0*q(QRHO)*dot_product(q(QU:QW),q(QU:QW)) &
							      + 0.5d0*(dot_product(q(QMAGX:QMAGZ),q(QMAGX:QMAGZ)))
			u(QMAGX:QMAGZ) = q(QMAGX:QMAGZ)
end subroutine PToC

end module hlld_solver
