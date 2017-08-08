module hlld_solver

implicit none

private hlldx, hlldy, hlldz, primtofluxx, primtofluxy, primtofluxz
public hlld

contains

subroutine hlld(qm,qp,qRd_l1,qRd_l2,qRd_l3,qRd_h1,qRd_h2,qRd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
				dir)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module, only: QVAR
implicit none 

	integer, intent(in)   :: qRd_l1,qRd_l2,qRd_l3,qRd_h1,qRd_h2,qRd_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	integer, intent(in)   :: dir

	real(rt), intent(in)  :: qm(qRd_l1:qRd_h1,qRd_l2:qRd_h2,qRd_l3:qRd_h3,QVAR)
	real(rt), intent(in)  :: qp(qRd_l1:qRd_h1,qRd_l2:qRd_h2,qRd_l3:qRd_h3,QVAR)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR)
	
	integer				  :: i, j, k

 flx = 0.d0 

 if(dir.eq.1) then
	do k = flx_l3, flx_h3
		do j = flx_l2, flx_h2
			do i = flx_l1, flx_h1
				call hlldx(qp(i,j,k,:),qm(i+1,j,k,:),flx(i,j,k,:))
				if(isnan(flx(i,j,k,2)).or. abs(flx(i,j,k,2)).ge. 1d14) then
					write(*,*) "Nan flux at ", i, j, k
					write(*,*) "q - =", qm(i+1,j,k,:)
					write(*,*) "q + =", qp(i,j,k,:)
					pause
				endif
			enddo
		enddo
	enddo
 elseif(dir.eq.2) then
	do k = flx_l3, flx_h3
		do j = flx_l2, flx_h2
			do i = flx_l1, flx_h1
				call hlldy(qp(i,j,k,:),qm(i,j+1,k,:),flx(i,j,k,:))
				if(isnan(flx(i,j,k,2))) then
					write(*,*) "Nan flux_y at ", i, j, k
					write(*,*) "q - =", qm(i,j+1,k,:)
					write(*,*) "q + =", qp(i,j,k,:)
					pause
				endif
			enddo
		enddo
	enddo
 else 
	do k = flx_l3, flx_h3
		do j = flx_l2, flx_h2
			do i = flx_l1, flx_h1
				call hlldz(qp(i,j,k,:),qm(i,j,k+1,:),flx(i,j,k,:))
			enddo
		enddo
	enddo
 endif
end subroutine hlld

!================================================= X Direction =======================================================
subroutine hlldx(qL,qR,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bx respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, caxL
	real(rt) 			  :: caR, caxR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)			  :: QsL(QVAR), FsL(QVAR)
	real(rt)			  :: QsR(QVAR), FsR(QVAR)
	real(rt)			  :: QssL(QVAR), FssL(QVAR)
	real(rt)			  :: QssR(QVAR), FssR(QVAR)
	integer				  :: i
	character(len=10)	  :: choice
!Riemann solve
	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	QsL = 0.d0
	QsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	QssL = 0.d0
	QssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0
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
	sL   = min(qL(QU) - cfL,qR(QU) - cfR)
	sR 	 = max(qL(QU) + cfL,qR(QU) + cfR)
	sM   = ((sR - qR(QU))*qR(QRHO)*qR(QU) - (sL - qL(QU))*qL(QRHO)*qL(QU) - qR(QPRES) + qL(QPRES))/((sR - qR(QU))*qR(QRHO) - (sL - qL(QU))*qL(QRHO))
	!Pressures in the Riemann Fan
	ptL  = qL(QPRES)
	ptR  = qR(QPRES)
	pst  = (sR - qR(QU))*qR(QRHO)*ptL - (sL - qL(QU))*qL(QRHO)*ptR + qL(QRHO)*qR(QRHO)*(sR - qR(QU))*(sL - qL(QU))*(qR(QU) - qL(QU))
	pst  = pst/((sR - qR(QU))*qR(QRHO) - (sL - qL(QU))*qL(QRHO))

	!------------------------------------------- * states-------------------------------------------------------------------------
	!density
	QsL(QRHO) = qL(QRHO)*((sL - qL(QU))/(sL - sM))
	QsR(QRHO) = qR(QRHO)*((sR - qR(QU))/(sR - sM))
	!velocities
	!X dir
	QsL(QU)    = sM
	QsR(QU)    = sM
	!Y dir
	QsL(QV)    = qL(QV) - qL(QMAGX)*qL(QMAGY)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2))
	QsR(QV)    = qR(QV) - qR(QMAGX)*qR(QMAGY)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qL(QMAGX)**2))
	!Z dir
	QsL(QW)    = qL(QW) - qL(QMAGX)*qL(QMAGZ)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2))
	QsR(QW)    = qR(QW) - qR(QMAGX)*qR(QMAGZ)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qL(QMAGX)**2))
	
	!Magnetic Fields
	!X dir
	QsL(QMAGX) = qL(QMAGX)
	QsR(QMAGX) = qL(QMAGX) 
	!Y dir
	QsL(QMAGY) = qL(QMAGY)*(qL(QRHO)*(sL - qL(QU))**2 - qL(QMAGX)**2)/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2)
	QsR(QMAGY) = qR(QMAGY)*(qR(QRHO)*(sR - qR(QU))**2 - qR(QMAGX)**2)/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qR(QMAGX)**2)
	!Z dir
	QsL(QMAGZ) = qL(QMAGZ)*(qL(QRHO)*(sL - qL(QU))**2 - qL(QMAGX)**2)/(qL(QRHO)*(sL - qL(QU))*(sL - sM) - qL(QMAGX)**2)
	QsR(QMAGZ) = qR(QMAGZ)*(qR(QRHO)*(sR - qR(QU))**2 - qR(QMAGX)**2)/(qR(QRHO)*(sR - qR(QU))*(sR - sM) - qR(QMAGX)**2)
	
	!Energy *Stored in Pressure slot
	QsL(QPRES) = (sL - qL(QU))*eL - ptL*qL(QU) + pst*sM + qL(QMAGX)*(qL(QU)*qL(QMAGX) + qL(QV)*qL(QMAGY) + qL(QW)*qL(QMAGZ) &
				  - (QsL(QU)*QsL(QMAGX) + QsL(QV)*QsL(QMAGY) + QsL(QW)*QsL(QMAGZ)))
	QsL(QPRES) = QsL(QPRES)/(sL - sM)
	QsR(QPRES) = (sR - qR(QU))*eR - ptR*qR(QU) + pst*sM + qR(QMAGX)*(qR(QU)*qR(QMAGX) + qR(QV)*qR(QMAGY) + qR(QW)*qR(QMAGZ) &
				  - (QsR(QU)*QsR(QMAGX) + QsR(QV)*QsR(QMAGY) + QsR(QW)*QsR(QMAGZ)))
	QsR(QPRES) = QsR(QPRES)/(sR - sM)
	!Hack
	do i = 1, QVAR
		if (isnan(QsL(i))) then
			QsL(i) = qL(i)
		endif
		if(isnan(QsR(i))) then
			QsR(i) = qR(i)
		endif
	enddo
	!speeds
	ssL = sM - abs(qL(QMAGX))/sqrt(QsL(QRHO))
	ssR = sM + abs(qR(QMAGX))/sqrt(QsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	QssL(QRHO)  = QsL(QRHO)
	QssR(QRHO)  = QsR(QRHO)
	!u
	QssL(QU)    = sM
	QssR(QU)    = sM
	!v
	QssL(QV)    = (sqrt(QsL(QRHO))*QsL(QV) + sqrt(QsR(QRHO))*qsR(QV) + (QsR(QMAGY) - QsL(QMAGY))*sign(1.d0,qL(QMAGX)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QV)    = QssL(QV)
	!w
	QssL(QW)    = (sqrt(QsL(QRHO))*QsL(QW) + sqrt(QsR(QRHO))*qsR(QW) + (QsR(QMAGZ) - QsL(QMAGZ))*sign(1.d0,qL(QMAGX)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QW)    = QssL(QW)
	!Bx
	QssL(QMAGX) = QsL(QMAGX)
	QssR(QMAGX) = QsR(QMAGX)
	!By
	QssL(QMAGY) = (sqrt(QsL(QRHO))*QsR(QMAGY) + sqrt(QsR(QRHO))*QsL(QMAGY) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QV) - QsL(QV))*sign(1.d0,QsR(QMAGX)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGY) = QssL(QMAGY)
	!Bz
	QssL(QMAGZ) = (sqrt(QsL(QRHO))*QsR(QMAGZ) + sqrt(QsR(QRHO))*QsL(QMAGZ) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QW) - QsL(QW))*sign(1.d0,QsR(QMAGX)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGZ) = QssL(QMAGZ)
	!Energy *Stored in Pressure slot
	QssL(QPRES) = QsL(QPRES) - sqrt(QsL(QRHO))*(dot_product(QsL(QU:QW),QsL(QMAGX:QMAGZ)) - dot_product(QssL(QU:QW),QssL(QMAGX:QMAGZ)))*sign(1.d0, QsL(QMAGX))
	QssR(QPRES) = QsR(QPRES) + sqrt(QsR(QRHO))*(dot_product(QsR(QU:QW),QsR(QMAGX:QMAGZ)) - dot_product(QssR(QU:QW),QssR(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGX))
	!Hack
	do i = 1, QVAR
		if (isnan(QssL(i))) then
			QssL(i) = QsL(i)
		endif
		if(isnan(QssR(i))) then
			QssR(i) = QsR(i)
		endif
	enddo
	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(QsL - qL)
	FssL = FL + ssL*QssL - (ssL - sL)*QsL - sL*qL
	FsR  = FR + sR*(QsR - qR)
	FssR = FR + ssR*QssR - (ssR - sR)*QsR - sR*qR
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
!		if(isnan(flx(i))) then
!			write(*,*) "Flux is nan in", i, "component"
!			write(*,*) "Flux = ", choice
!			write(*,*) "FL = ", FL
!			write(*,*) "FR = ", FR
!			write(*,*) "QL = ", qL
!			write(*,*) "QR = ", qR
!			write(*,*) "QsL = ", QsL
!			write(*,*) "QsR = ", QsR
!			write(*,*) "QssL = ", QssL
!			write(*,*) "QssR = ", QssR
!			pause
!			return
!		endif
!	enddo

	!Rusanof flux
	flx = 0.5d0*(FL + FR) + 0.5d0*sM*(QL - QR)
end subroutine hlldx

!============================================================= Y Direction =================================================================

subroutine hlldy(qR,qL,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/By respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, cayL
	real(rt) 			  :: caR, cayR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)			  :: QsL(QVAR), FsL(QVAR)
	real(rt)			  :: QsR(QVAR), FsR(QVAR)
	real(rt)			  :: QssL(QVAR), FssL(QVAR)
	real(rt)			  :: QssR(QVAR), FssR(QVAR)
	integer				  :: i
	character(len=10)	  :: choice
!Riemann solve

	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	QsL = 0.d0
	QsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	QssL = 0.d0
	QssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0
	call primtofluxy(qL, FL)
	call primtofluxy(qR, FR)
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
	QsL(QRHO) = qL(QRHO)*((sL - qL(QV))/(sL - sM))
	QsR(QRHO) = qR(QRHO)*((sR - qR(QV))/(sR - sM))
	!velocities
	!X dir
	QsL(QU)    = qL(QU) - qL(QMAGY)*qL(QMAGX)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2))
	QsR(QU)    = qR(QU) - qR(QMAGY)*qR(QMAGX)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qR(QMAGY)**2))
	!Y dir
	QsL(QV)    = sM
	QsR(QV)    = sM
	!Z dir
	QsL(QW)    = qL(QW) - qL(QMAGY)*qL(QMAGZ)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2))
	QsR(QW)    = qR(QW) - qR(QMAGY)*qR(QMAGZ)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qL(QMAGY)**2))
	
	!Magnetic Fields
	!X dir
	QsL(QMAGX) = qL(QMAGX)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2)
	QsR(QMAGX) = qR(QMAGX)*(qR(QRHO)*(sR - qR(QV))**2 - qL(QMAGY)**2)/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qL(QMAGY)**2)
	!Y dir
	QsL(QMAGY) = qL(QMAGY)
	QsR(QMAGY) = qL(QMAGY) 
	!Z dir
	QsL(QMAGZ) = qL(QMAGZ)*(qL(QRHO)*(sL - qL(QV))**2 - qL(QMAGY)**2)/(qL(QRHO)*(sL - qL(QV))*(sL - sM) - qL(QMAGY)**2)
	QsR(QMAGZ) = qR(QMAGZ)*(qR(QRHO)*(sR - qR(QV))**2 - qL(QMAGY)**2)/(qR(QRHO)*(sR - qR(QV))*(sR - sM) - qL(QMAGY)**2)
	
	!Energy *Stored in Pressure slot
	QsL(QPRES) = (sL - qL(QV))*eL - ptL*qL(QV) + pst*sM + qL(QMAGY)*(qL(QU)*qL(QMAGX) + qL(QV)*qL(QMAGY) + qL(QW)*qL(QMAGZ) &
				  - (QsL(QU)*QsL(QMAGX) + QsL(QV)*QsL(QMAGY) + QsL(QW)*QsL(QMAGZ)))
	QsL(QPRES) = QsL(QPRES)/(sL - sM)
	QsR(QPRES) = (sR - qR(QV))*eR - ptR*qR(QV) + pst*sM + qR(QMAGY)*(qR(QU)*qR(QMAGX) + qR(QV)*qR(QMAGY) + qR(QW)*qR(QMAGZ) &
				  - (QsR(QU)*QsR(QMAGX) + QsR(QV)*QsR(QMAGY) + QsR(QW)*QsR(QMAGZ)))
	QsR(QPRES) = QsR(QPRES)/(sR - sM)
	!Hack
	do i = 1, QVAR
		if (isnan(QsL(i))) then
			QsL(i) = qL(i)
		endif
		if(isnan(QsR(i))) then
			QsR(i) = qR(i)
		endif
	enddo
	!speeds
	ssL = sM - abs(qL(QMAGY))/sqrt(QsL(QRHO))
	ssR = sM + abs(qR(QMAGY))/sqrt(QsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	QssL(QRHO)  = QsL(QRHO)
	QssR(QRHO)  = QsR(QRHO)
	!u
	QssL(QU)    = (sqrt(QsL(QRHO))*QsL(QU) + sqrt(QsR(QRHO))*qsR(QU) + (QsR(QMAGX) - QsL(QMAGX))*sign(1.d0,qL(QMAGY)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QU)    = QssL(QU)
	!v
	QssL(QV)    = sM
	QssR(QV)    = sM
	!w
	QssL(QW)    = (sqrt(QsL(QRHO))*QsL(QW) + sqrt(QsR(QRHO))*qsR(QW) + (QsR(QMAGZ) - QsL(QMAGZ))*sign(1.d0,qL(QMAGY)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QW)    = QssL(QW)
	!Bx
	QssL(QMAGX) = (sqrt(QsL(QRHO))*QsR(QMAGX) + sqrt(QsR(QRHO))*QsL(QMAGX) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QU) - QsL(QU))*sign(1.d0,QsR(QMAGY)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGX) = QssL(QMAGX)
	!By
	QssL(QMAGY) = QsL(QMAGY)
	QssR(QMAGY) = QsR(QMAGY)
	!Bz
	QssL(QMAGZ) = (sqrt(QsL(QRHO))*QsR(QMAGZ) + sqrt(QsR(QRHO))*QsL(QMAGZ) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QW) - QsL(QW))*sign(1.d0,QsR(QMAGY)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGZ) = QssL(QMAGZ)
	!Energy *Stored in Pressure slot
	QssL(QPRES) = QsL(QPRES) - sqrt(QsL(QRHO))*(dot_product(QsL(QU:QW),QsL(QMAGX:QMAGZ)) - dot_product(QssL(QU:QW),QssL(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGY))
	QssR(QPRES) = QsR(QPRES) + sqrt(QsR(QRHO))*(dot_product(QsR(QU:QW),QsR(QMAGX:QMAGZ)) - dot_product(QssR(QU:QW),QssR(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGY))
	
	!Hack
	do i = 1, QVAR
		if (isnan(QssL(i))) then
			QssL(i) = QsL(i)
		endif
		if(isnan(QssR(i))) then
			QssR(i) = QsR(i)
		endif
	enddo
	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(QsL - qL)
	FssL = FL + ssL*QssL - (ssL - sL)*QsL - sL*qL
	FsR  = FR + sR*(QsR - qR)
	FssR = FR + ssR*QssR - (ssR - sR)*QsR - sR*qR
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
	!Rusanof flux
	flx = 0.5d0*(FL + FR) + 0.5d0*sM*(QL - QR)
end subroutine hlldy

!============================================================= Z Direction =================================================================

subroutine hlldz(qR,qL,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bz respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module

implicit none
	real(rt), intent(in)  :: qL(QVAR)
	real(rt), intent(in)  :: qR(QVAR)
	real(rt), intent(inout) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, cazL
	real(rt) 			  :: caR, cazR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)			  :: QsL(QVAR), FsL(QVAR)
	real(rt)			  :: QsR(QVAR), FsR(QVAR)
	real(rt)			  :: QssL(QVAR), FssL(QVAR)
	real(rt)			  :: QssR(QVAR), FssR(QVAR)
	integer				  :: i
	character(len=10)	  :: choice

!Riemann solve
	flx = 0.d0
	FL  = 0.d0
	FR  = 0.d0
	QsL = 0.d0
	QsR = 0.d0
	FsL = 0.d0
	FsR = 0.d0
	QssL = 0.d0
	QssR = 0.d0
	FssL = 0.d0
	FssR = 0.d0

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
	QsL(QRHO) = qL(QRHO)*((sL - qL(QW))/(sL - sM))
	QsR(QRHO) = qR(QRHO)*((sR - qR(QW))/(sR - sM))
	!velocities
	!X dir
	QsL(QU)    = qL(QU) - qL(QMAGZ)*qL(QMAGX)*((sM - qL(QU))/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2))
	QsR(QU)    = qR(QU) - qR(QMAGZ)*qR(QMAGX)*((sM - qR(QU))/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2))
	!Y dir
	QsL(QV)    = qL(QV) - qL(QMAGZ)*qL(QMAGY)*((sM - qL(QV))/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2))
	QsR(QV)    = qR(QV) - qR(QMAGZ)*qR(QMAGY)*((sM - qR(QV))/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2))
	!Z dir
	QsL(QW)    = sM
	QsR(QW)    = sM
	
	!Magnetic Fields
	!X dir
	QsL(QMAGX) = qL(QMAGX)*(qL(QRHO)*(sL - qL(QW))**2 - qL(QMAGZ)**2)/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2)
	QsR(QMAGX) = qR(QMAGX)*(qR(QRHO)*(sR - qR(QW))**2 - qL(QMAGZ)**2)/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2)
	!Y dir
	QsL(QMAGY) = qL(QMAGY)*(qL(QRHO)*(sL - qL(QW))**2 - qL(QMAGZ)**2)/(qL(QRHO)*(sL - qL(QW))*(sL - sM) - qL(QMAGZ)**2)
	QsR(QMAGY) = qR(QMAGY)*(qR(QRHO)*(sR - qR(QW))**2 - qL(QMAGZ)**2)/(qR(QRHO)*(sR - qR(QW))*(sR - sM) - qL(QMAGZ)**2)
	!Z dir
	QsL(QMAGZ) = qL(QMAGZ)
	QsR(QMAGZ) = qL(QMAGZ) 
	
	!Energy *Stored in Pressure slot
	QsL(QPRES) = (sL - qL(QW))*eL - ptL*qL(QW) + pst*sM + qL(QMAGZ)*(qL(QU)*qL(QMAGX) + qL(QV)*qL(QMAGY) + qL(QW)*qL(QMAGZ) &
				  - (QsL(QU)*QsL(QMAGX) + QsL(QV)*QsL(QMAGY) + QsL(QW)*QsL(QMAGZ)))
	QsL(QPRES) = QsL(QPRES)/(sL - sM)
	QsR(QPRES) = (sR - qR(QW))*eR - ptR*qR(QW) + pst*sM + qR(QMAGZ)*(qR(QU)*qR(QMAGX) + qR(QV)*qR(QMAGY) + qR(QW)*qR(QMAGZ) &
				  - (QsR(QU)*QsR(QMAGX) + QsR(QV)*QsR(QMAGY) + QsR(QW)*QsR(QMAGZ)))
	QsR(QPRES) = QsR(QPRES)/(sR - sM)
	!Hack
	do i = 1, QVAR
		if (isnan(QsL(i))) then
			QsL(i) = qL(i)
		endif
		if(isnan(QsR(i))) then
			QsR(i) = qR(i)
		endif
	enddo
	!speeds
	ssL = sM - abs(qL(QMAGZ))/sqrt(QsL(QRHO))
	ssR = sM + abs(qR(QMAGZ))/sqrt(QsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	QssL(QRHO)  = QsL(QRHO)
	QssR(QRHO)  = QsR(QRHO)
	!u
	QssL(QU)    = (sqrt(QsL(QRHO))*QsL(QU) + sqrt(QsR(QRHO))*qsR(QU) + (QsR(QMAGX) - QsL(QMAGX))*sign(1.d0,qL(QMAGZ)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QU)    = QssL(QU)
	!v
	QssL(QV)    = (sqrt(QsL(QRHO))*QsL(QV) + sqrt(QsR(QRHO))*qsR(QV) + (QsR(QMAGY) - QsL(QMAGY))*sign(1.d0,qL(QMAGZ)))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QV)    = QssL(QW)
	!w
	QssL(QW)    = sM
	QssR(QW)    = sM
	!Bx
	QssL(QMAGX) = (sqrt(QsL(QRHO))*QsR(QMAGX) + sqrt(QsR(QRHO))*QsL(QMAGX) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QU) - QsL(QU))*sign(1.d0,QsR(QMAGZ)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGX) = QssL(QMAGX)
	!By
	QssL(QMAGY) = (sqrt(QsL(QRHO))*QsR(QMAGY) + sqrt(QsR(QRHO))*QsL(QMAGY) + sqrt(QsL(QRHO)*QsR(QRHO))*(QsR(QV) - QsL(QV))*sign(1.d0,QsR(QMAGZ)))&
				   /(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QMAGY) = QssL(QMAGY)
	!Bz
	QssL(QMAGZ) = QsL(QMAGZ)
	QssR(QMAGZ) = QsR(QMAGZ)
	!Energy *Stored in Pressure slot
	QssL(QPRES) = QsL(QPRES) - sqrt(QsL(QRHO))*(dot_product(QsL(QU:QW),QsL(QMAGX:QMAGZ)) - dot_product(QssL(QU:QW),QssL(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGZ))
	QssR(QPRES) = QsR(QPRES) + sqrt(QsR(QRHO))*(dot_product(QsR(QU:QW),QsR(QMAGX:QMAGZ)) - dot_product(QssR(QU:QW),QssR(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGZ))
	!Hack
	do i = 1, QVAR
		if (isnan(QssL(i))) then
			QssL(i) = QsL(i)
		endif
		if(isnan(QssR(i))) then
			QssR(i) = QsR(i)
		endif
	enddo
	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(QsL - qL)
	FssL = FL + ssL*QssL - (ssL - sL)*QsL - sL*qL
	FsR  = FR + sR*(QsR - qR)
	FssR = FR + ssR*QssR - (ssR - sR)*QsR - sR*qR
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
!		if(isnan(flx(i))) then
!			write(*,*) "Flux is nan in", i, "component"
!			write(*,*) "Flux = ", choice
!			write(*,*) "FL = ", FL
!			write(*,*) "FR = ", FR
!			write(*,*) "QL = ", qL
!			write(*,*) "QR = ", qR
!			write(*,*) "QsL = ", QsL
!			write(*,*) "QsR = ", QsR
!			write(*,*) "QssL = ", QssL
!			write(*,*) "QssR = ", QssR
!			pause
!			return
!		endif
!	enddo
	!Rusanof flux
	flx = 0.5d0*(FL + FR) + 0.5d0*sM*(QL - QR)
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
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(gamma_minus_1)& 
			   + 0.5d0*Q(QRHO)*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QU)
	F(UMX)	 = Q(QRHO)*Q(QU)**2 + Q(QPRES) - Q(QMAGX)**2
	F(UMY)   = Q(QRHO)*Q(QU)*Q(QV) - Q(QMAGX)*Q(QMAGY)
	F(UMZ)   = Q(QRHO)*Q(QU)*Q(QW) - Q(QMAGX)*Q(QMAGZ)
	F(UEDEN) = Q(QU)*(e + Q(QPRES)) -Q(QMAGX)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = 0.d0
	F(QMAGY) = Q(QU)*Q(QMAGY) - Q(QMAGX)*Q(QV)
	F(QMAGZ) = Q(QU)*Q(QMAGZ) - Q(QMAGX)*Q(QW)

end subroutine primtofluxx

!-------------------------------------- Y Direction ------------------------------------------------------------

subroutine primtofluxy(Q, F)
 use amrex_fort_module, only : rt => amrex_real
 use meth_params_module
implicit none

	real(rt), intent(in)  :: Q(QVAR)
	real(rt), intent(out) :: F(QVAR)
	real(rt)			  :: e

	F = 0.d0
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(gamma_minus_1)& 
			   + 0.5d0*Q(QRHO)*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QV)
	F(UMX)	 = Q(QRHO)*Q(QU)*Q(QV) - Q(QMAGX)*Q(QMAGY)
	F(UMY)   = Q(QRHO)*Q(QV)**2 + Q(QPRES) - Q(QMAGY)**2
	F(UMZ)   = Q(QRHO)*Q(QV)*Q(QW) - Q(QMAGY)*Q(QMAGZ)
	F(UEDEN) = Q(QV)*(e + Q(QPRES)) -Q(QMAGY)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = Q(QV)*Q(QMAGX) - Q(QMAGY)*Q(QU)
	F(QMAGY) = 0.d0
	F(QMAGZ) = Q(QV)*Q(QMAGZ) - Q(QMAGY)*Q(QW)

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
	e 		 = (Q(QPRES)-0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ)))/(gamma_minus_1)& 
			   + 0.5d0*Q(QRHO)*dot_product(Q(QU:QW),Q(QU:QW)) + 0.5d0*dot_product(Q(QMAGX:QMAGZ),Q(QMAGX:QMAGZ))
	F(URHO)  = Q(QRHO)*Q(QW)
	F(UMX)	 = Q(QRHO)*Q(QW)*Q(QU) - Q(QMAGX)*Q(QMAGZ)
	F(UMY)   = Q(QRHO)*Q(QW)*Q(QV) - Q(QMAGY)*Q(QMAGZ)
	F(UMZ)   = Q(QRHO)*Q(QW)**2 + Q(QPRES) - Q(QMAGZ)**2
	F(UEDEN) = Q(QW)*(e + Q(QPRES)) -Q(QMAGZ)*dot_product(Q(QMAGX:QMAGZ),Q(QU:QW))
	F(QMAGX) = Q(QW)*Q(QMAGX) - Q(QMAGZ)*Q(QU)
	F(QMAGY) = Q(QW)*Q(QMAGY) - Q(QMAGZ)*Q(QV)
	F(QMAGZ) = 0.d0

end subroutine primtofluxz

end module hlld_solver
