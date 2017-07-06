module hlld_solver

implicit none

private hlldx,hlldy,hlldz, primtoflux
public hlld


subroutine hlld(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
				dir)
 use amrex_fort_module, only : rt => amrex_real
implicit none 

	integer, intent(in)   :: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	integer, intent(in)   :: dir

	real(rt), intent(in)  :: qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
	real(rt), intent(in)  :: qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)

 if(dir.eq.1) then
	do k = qpd_l3, qpd_h3
		do j = qpd_l2, qpd_h2
			do i = qpq_l1, qpd_h1
				call hlldx(qm(i,j,k,:),qp(i,j,k,:),flx(i,j,k,:))
			enddo
		enddo
	enddo
 elseif(dir.eq.2) then
	do k = qpd_l3, qpd_h3
		do j = qpd_l2, qpd_h2
			do i = qpq_l1, qpd_h1
				call hlldy(qm(i,j,k,:),qp(i,j,k,:),flx(i,j,k,:))
			enddo
		enddo
	enddo
 else 
	do k = qpd_l3, qpd_h3
		do j = qpd_l2, qpd_h2
			do i = qpq_l1, qpd_h1
				call hlldz(qm(i,j,k,:),qp(i,j,k,:),flx(i,j,k,:))
			enddo
		enddo
	enddo
 endif
end subroutine hlld

!================================================= X Direction =======================================================
subroutine hlldx(qm,qp,flx)

!Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bx respectively. 
!Total Pressure is constant throughout the Riemann fan, pst!

 use amrex_fort_module, only : rt => amrex_real
 use meth_mhd_params_module

implicit none
	real(rt), intent(in)  :: qm(QVAR)
	real(rt), intent(in)  :: qp(QVAR)
	real(rt), intent(out) :: flx(QVAR)

	real(rt)			  :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, caxL
	real(rt) 			  :: caR, caxR, asL, asR, ptL, ptR, eL, eR
	real(rt)			  :: FL(QVAR), FR(QVAR)
	real(rt)			  :: QsL(QVAR), FsL(QVAR)
	real(rt)			  :: QsR(QVAR), FsR(QVAR)
	real(rt)			  :: QssL(QVAR), FssL(QVAR)
	real(rt)			  :: QssR(QVAR), FssR(QVAR)

!Riemann solve
	call primtoflux(qm, FL)
	call primtoflux(qp, FR)
	
	eL   = qm(QPRES)/(qm(QRHO)*gamma_minus_1)
	eR   = qp(QPRES)/(qp(QRHO)*gamma_minus_1)
	asL  = gamma_const * qm(QPRES)/qm(QRHO)
	asR  = gamma_const * qp(QPRES)/qp(QRHO)
	caL  = (qm(QMAGX)**2 + qm(QMAGY)**2 + qm(QMAGZ)**2)/qm(QRHO) !Magnetic Speeds
	caR  = (qp(QMAGX)**2 + qp(QMAGY)**2 + qp(QMAGZ)**2)/qp(QRHO)
	caxL = (qm(QMAGX)**2)/qm(QRHO)
	caxR = (qp(QMAGX)**2)/qp(QRHO)
	!Catch the fastest waves, brah
	cfL  = 0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*caL)
	cfR  = 0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*caR)
	!Riemann Speeds
	sL   = min(qm(QU),qp(QU)) - max(cfL,cfR)
	sR 	 = max(qm(QU),qp(QU)) + max(cfL,cfR)
	sM   = ((sR - qp(QU))*qp(QRHO)*qp(QU) - (sL - qm(QU))*qm(QRHO)*qm(QU) - qp(QPRES) + qm(QPRES))/((sR - qp(QU))*qp(QRHO) - (sL - qm(QU))*qm(QRHO))
	!Pressures in the Riemann Fan
	ptL  = qm(PRES) + 0.5d0*(qm(QMAGX)**2 + qm(QMAGY)**2 + qm(QMAGZ)**2)
	ptR  = qp(PRES) + 0.5d0*(qp(QMAGX)**2 + qp(QMAGY)**2 + qp(QMAGZ)**2)
	pst  = (sR - qp(QU))*qp(QRHO)*ptL - (sL - qm(QU))*qm(QRHO)*ptR + qm(QRHO)*qp(QRHO)*(sR - qp(QU))*(sL - qm(QU))*(qp(QU) - qm(QU))
	pst  = pst/((sR - qp(QU))*qp(QRHO) - (sL - qm(QU))*qm(QRHO))

	!------------------------------------------- * states-------------------------------------------------------------------------
	!density
	QsL(QRHO) = qm(QRHO)*((sL - qm(QU))/(sL - sM))
	QsR(QRHO) = qp(QRHO)*((sR - qp(QU))/(sR - sM))
	!velocities
	!X dir
	QsL(QU)    = sM
	QsR(QU)    = sM
	!Y dir
	QsL(QV)    = qm(QV) - qm(QMAGX)*qm(QMAGY)*((sM - qm(QU))/(qm(QRHO)*(sL - sM) - qm(QMAGX)**2)
	QsR(QV)    = qp(QV) - qp(QMAGX)*qp(QMAGY)*((sM - qp(QU))/(qp(QRHO)*(sR - sM) - qm(QMAGX)**2)
	!Z dir
	QsL(QW)    = qm(QW) - qm(QMAGX)*qm(QMAGZ)*((sM - qm(QU))/(qm(QRHO)*(sL - sM) - qm(QMAGX)**2)
	QsR(QW)    = qp(QW) - qp(QMAGX)*qp(QMAGZ)*((sM - qp(QU))/(qp(QRHO)*(sR - sM) - qm(QMAGX)**2)
	
	!Magnetic Fields
	!X dir
	QsL(QMAGX) = qm(QMAGX)
	QsR(QMAGX) = qm(QMAGX) 
	!Y dir
	QsL(QMAGY) = qm(QMAGY)*(qm(QRHO)*(sL - qm(QU))**2 - qm(QMAGX)**2)/(qm(QRHO)*(sL - qm(QU))*(sL - sM) - qm(QMAGX)**2)
	QsR(QMAGY) = qp(QMAGY)*(qp(QRHO)*(sR - qp(QU))**2 - qm(QMAGX)**2)/(qp(QRHO)*(sR - qp(QU))*(sR - sM) - qm(QMAGX)**2)
	!Z dir
	QsL(QMAGZ) = qm(QMAGZ)*(qm(QRHO)*(sL - qm(QU))**2 - qm(QMAGX)**2)/(qm(QRHO)*(sL - qm(QU))*(sL - sM) - qm(QMAGX)**2)
	QsR(QMAGZ) = qp(QMAGZ)*(qp(QRHO)*(sR - qp(QU))**2 - qm(QMAGX)**2)/(qp(QRHO)*(sR - qp(QU))*(sR - sM) - qm(QMAGX)**2)
	
	!Energy *Stored in Pressure slot
	QsL(QPRES) = (sL - qm(QU))*eL - ptL*qm(QU) + pst*sM + qm(QMAGX)*(qm(QU)*qm(QMAGX) + qm(QV)*qm(QMAGY) + qm(QW)*qm(QMAGZ) &
				  - (QsL(QU)*QsL(QMAGX) + QsL(QV)*QsL(QMAGY) + QsL(QW)*QsL(QMAGZ))
	QsL(QPRES) = QsL(QPRES)/(sL - sM)
	QsR(QPRES) = (sR - qp(QU))*eR - ptR*qp(QU) + pst*sM + qp(QMAGX)*(qp(QU)*qp(QMAGX) + qp(QV)*qp(QMAGY) + qp(QW)*qp(QMAGZ) &
				  - (QsR(QU)*QsR(QMAGX) + QsR(QV)*QsR(QMAGY) + QsR(QW)*QsR(QMAGZ))
	QsR(QPRES) = QsR(QPRES)/(sR - sM)

	!speeds
	ssL = sM - abs(qm(QMAGX))/sqrt(QsL(QRHO))
	ssR = sM - abs(qp(QMAGX))/sqrt(QsR(QRHO))

	!----------------------------------------- ** states ------------------------------------------------------------------------------
	!Dens
	QssL(QRHO)  = QsL(QRHO)
	QssR(QRHO)  = QsR(QRHO)
	!u
	QssL(QU)    = sM
	QssR(QU)    = sM
	!v
	QssL(QV)    = (sqrt(QsL(QRHO))*QsL(QV) + sqrt(QsR(QRHO))*qsR(QV) + (QsR(QMAGY) - QsL(QMAGY))*sign(1.d0,qm(QMAGX))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
	QssR(QV)    = QssL(QV)
	!w
	QssL(QW)    = (sqrt(QsL(QRHO))*QsL(QW) + sqrt(QsR(QRHO))*qsR(QW) + (QsR(QMAGZ) - QsL(QMAGZ))*sign(1.d0,qm(QMAGX))/(sqrt(QsL(QRHO)) + sqrt(QsR(QRHO)))
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
	QssL(QPRES) = QsL(QPRES) - sqrt(QsL(QRHO))*(dot_product(QsL(QU:QW),QsL(QMAGX:QMAGZ)) - dot_product(QssL(QU:QW),QssL(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGX))
	QssR(QPRES) = QsR(QPRES) - sqrt(QsR(QRHO))*(dot_product(QsR(QU:QW),QsR(QMAGX:QMAGZ)) - dot_product(QssR(QU:QW),QssR(QMAGX:QMAGZ)))*sign(1.d0, QsR(QMAGX))

	!--------------------------------------------------------- Fluxes ----------------------------------------------------------------------
	FsL  = FL + sL*(QsL - qm)
	FssL = FL + ssL*UssL - (ssL - sL)*UsL - sL*qm
	FsR  = FR + sR*(QsR - qp)
	FssR = FR + ssR*UssR - (ssR - sR)*UsR - sR*qp
	!Solve the RP
	if(sL .gt. 0.d0) then
	flux = FL
	elseif(sL .le. 0.d0 .and. ssL .ge. 0.d0) then
	flux = FsL
	elseif(ssl .le. 0.d0 .and. sM .ge. 0.d0) then
	flux = FssL
	elseif(sM .le. 0.d0 .and. ssR .ge. 0.d0) then
	flux = FssR
	elseif(ssR .le. 0.d0 .and. sR .ge. 0.d0) then
	flux = FsR
	else 
	flux = FR
	endif
end subroutine hlldx
end module hlld_solver
