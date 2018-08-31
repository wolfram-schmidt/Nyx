module hlld_solver

   implicit none
   public  hlld
   private PToC
contains

subroutine hlld(work_lo, work_hi, qm ,qp ,q_l1 ,q_l2 ,q_l3 ,q_h1 ,q_h2 ,q_h3 , &
                flx ,flx_l1 ,flx_l2 ,flx_l3 ,flx_h1 ,flx_h2 ,flx_h3 ,dir)

  !Riemann solve:
  !Main assumption, the normal velocity/Mag field is constant in the Riemann fan, and is sM/Bn respectively. 
  !Total Pressure is constant throughout the Riemann fan, pst!

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module

   integer, intent(in)   :: q_l1,q_l2,q_l3,q_h1,q_h2,q_h3
   integer, intent(in)   :: work_lo(3), work_hi(3)
   integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
   integer, intent(in)   :: dir

   real(rt), intent(in)  :: qm(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
   real(rt), intent(in)  :: qp(q_l1:q_h1,q_l2:q_h2,q_l3:q_h3,QVAR,3)
   real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,QVAR)

   real(rt)   :: cfL, cfR, sL, sR, sM, ssL, ssR, pst, caL, canL
   real(rt)   :: caR, canR, asL, asR, ptL, ptR, eL, eR
   real(rt)   :: QL(QVAR), QR(QVAR)
   real(rt)   :: FL(QVAR), FR(QVAR)
   real(rt)   :: uL(QVAR), uR(QVAR)
   real(rt)   :: UsL(QVAR), FsL(QVAR)
   real(rt)   :: UsR(QVAR), FsR(QVAR)
   real(rt)   :: UssL(QVAR), FssL(QVAR)
   real(rt)   :: UssR(QVAR), FssR(QVAR)

   integer           :: QVELN, QVELP1, QVELP2
   integer           :: QMAGN, QMAGP1, QMAGP2
   integer           :: UMN  , UMP1  , UMP2
   integer           :: i,j,k
   
   if (dir .eq. 1) then
      QMAGN  = QMAGX
      QMAGP1 = QMAGY
      QMAGP2 = QMAGZ
      QVELN  = QU
      QVELP1 = QV
      QVELP2 = QW
      UMN    = UMX
      UMP1   = UMY
      UMP2   = UMZ
   else if (dir .eq. 2) then
      QMAGN  = QMAGY
      QMAGP1 = QMAGZ
      QMAGP2 = QMAGX
      QVELN  = QV
      QVELP1 = QW
      QVELP2 = QU
      UMN    = UMY
      UMP1   = UMZ
      UMP2   = UMX
   else if (dir .eq. 3) then
      QMAGN  = QMAGZ
      QMAGP1 = QMAGX
      QMAGP2 = QMAGY
      QVELN  = QW
      QVELP1 = QU
      QVELP2 = QV
      UMN    = UMZ
      UMP1   = UMX
      UMP2   = UMY
   end if

   do k = work_lo(3), work_hi(3)
    do j = work_lo(2), work_hi(2)
     do i = work_lo(1), work_hi(1)

      
      if (dir .eq. 1) then
         qL(:) = qp(i-1,j,k,:,dir)
      else if (dir .eq. 2) then
         qL(:) = qp(i,j-1,k,:,dir)
      else if (dir .eq. 3) then
         qL(:) = qp(i,j,k-1,:,dir)
      end if

      qR(:) = qm(i,j,k,:,dir)

      flx(i,j,k,:) = 0.d0  
      FL  = 0.d0; FR = 0.d0; UsL = 0.d0; UsR = 0.d0; FsL = 0.d0; FsR = 0.d0; UssL = 0.d0; UssR = 0.d0; FssL = 0.d0; FssR = 0.d0

      call PToC(qL,uL)
      call PToC(qR,uR)

      ! Total Pressures
      ptL  = qL(QPRES) + 0.5d0*(qL(QMAGX)**2 + qL(QMAGY)**2 + qL(QMAGZ)**2)
      ptR  = qR(QPRES) + 0.5d0*(qR(QMAGX)**2 + qR(QMAGY)**2 + qR(QMAGZ)**2)

      ! Note this is actually (rho E)
!      eL   = qL(QPRES)/(gamma_minus_1) &
!             + 0.5d0*dot_product(qL(QMAGX:QMAGZ),qL(QMAGX:QMAGZ)) &
!             + 0.5d0*dot_product(qL(QU:QW),qL(QU:QW))*qL(QRHO)

      eL         = uL(UEDEN)
      FL(URHO)   = qL(QRHO)*qL(QVELN)
      FL(UMN)    = qL(QRHO)*qL(QVELN)**2 + ptL - qL(QMAGN)**2
      FL(UMP1)   = qL(QRHO)*qL(QVELN)*qL(QVELP1) - qL(QMAGN)*qL(QMAGP1)
      FL(UMP2)   = qL(QRHO)*qL(QVELN)*qL(QVELP2) - qL(QMAGN)*qL(QMAGP2)
      FL(UEDEN)  = qL(QVELN)*(eL + ptL) - qL(QMAGN)*dot_product(qL(QMAGX:QMAGZ),qL(QU:QW))
      FL(QMAGN)  = 0.d0
      FL(QMAGP1) = qL(QVELN)*qL(QMAGP1) - qL(QVELP1)*qL(QMAGN)
      FL(QMAGP2) = qL(QVELN)*qL(QMAGP2) - qL(QVELP2)*qL(QMAGN)  

      ! Note this is actually (rho E)
!      eR   = qR(QPRES)/(gamma_minus_1) &
!             + 0.5d0*dot_product(qR(QMAGX:QMAGZ),qR(QMAGX:QMAGZ)) &
!             + 0.5d0*dot_product(qR(QU:QW),qR(QU:QW))*qR(QRHO)

      eR        = uR(UEDEN)
      FR(URHO)  = qR(QRHO)*qR(QVELN)
      FR(UMN)   = qR(QRHO)*qR(QVELN)**2 + ptR - qR(QMAGN)**2
      FR(UMP1)  = qR(QRHO)*qR(QVELN)*qR(QVELP1) - qR(QMAGN)*qR(QMAGP1)
      FR(UMP2)  = qR(QRHO)*qR(QVELN)*qR(QVELP2) - qR(QMAGN)*qR(QMAGP2)
      FR(UEDEN) = qR(QVELN)*(eR + ptR) - qR(QMAGN)*dot_product(qR(QMAGX:QMAGZ),qR(QU:QW))
      FR(QMAGN) = 0.d0
      FR(QMAGP1) = qR(QVELN)*qR(QMAGP1) - qR(QVELP1)*qR(QMAGN)
      FR(QMAGP2) = qR(QVELN)*qR(QMAGP2) - qR(QVELP2)*qR(QMAGN)  

      asL  = gamma_const*qL(QPRES)/qL(QRHO)
      asR  = gamma_const*qR(QPRES)/qR(QRHO)

      caL  = (qL(QMAGN)**2 + qL(QMAGP1)**2 + qL(QMAGP2)**2)/qL(QRHO)!Magnetic Speeds squared
      caR  = (qR(QMAGN)**2 + qR(QMAGP1)**2 + qR(QMAGP2)**2)/qR(QRHO)

      canL = (qL(QMAGN)**2)/qL(QRHO)
      canR = (qR(QMAGN)**2)/qR(QRHO)

!Catch the fastest waves, brah
      cfL  = sqrt(0.5d0*((asL + caL) + sqrt((asL + caL)**2 - 4.0d0*asL*canL)))
      cfR  = sqrt(0.5d0*((asR + caR) + sqrt((asR + caR)**2 - 4.0d0*asR*canR)))

!Riemann Speeds
      sL   = min(qL(QVELN) - cfL,qR(QVELN) - cfR)
      sR   = max(qL(QVELN) + cfL,qR(QVELN) + cfR)
!      sL   = min(qL(QVELN),qR(QVELN)) - max(cfL, cfR)
!      sR   = max(qL(QVELN),qR(QVELN)) + max(cfL, cfR)
      sM   = (sR - qR(QVELN))*qR(QRHO)*qR(QVELN) - (sL - qL(QVELN))*qL(QRHO)*qL(QVELN) - ptR + ptL
      sM   = sM/((sR - qR(QVELN))*qR(QRHO) - (sL - qL(QVELN))*qL(QRHO))

!Pressures in the Riemann Fan
      pst  = (sR - qR(QVELN))*qR(QRHO)*ptL - (sL - qL(QVELN))*qL(QRHO)*ptR & 
           + qL(QRHO)*qR(QRHO)*(sR - qR(QVELN))*(sL - qL(QVELN))*(qR(QVELN) - qL(QVELN))
      pst  = pst/((sR - qR(QVELN))*qR(QRHO) - (sL - qL(QVELN))*qL(QRHO))

!------------------------------------------- Density * states-------------------------------------------------------------------------

! Density
      UsL(QRHO) = qL(QRHO)*((sL - qL(QVELN))/(sL - sM))
      UsR(QRHO) = qR(QRHO)*((sR - qR(QVELN))/(sR - sM))

!--------------------------------------------------------- Vel * states----------------------------------------------------------------

!Normal dir
      UsL(QVELN)    = sM
      UsR(QVELN)    = sM

!Perpendicular dir
      if(abs(qL(QMAGN)*qL(QMAGP1)*(sM - qL(QVELN))).lt. 1d-14) then
        UsL(QVELP1)    = qL(QVELP1)
      else
        UsL(QVELP1)    = qL(QVELP1) - qL(QMAGN)*qL(QMAGP1)*((sM - qL(QVELN))/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2))
      endif
      if(abs(qR(QMAGN)*qR(QMAGP1)*(sM - qR(QVELN))).lt. 1d-14) then
        UsR(QVELP1) = qR(QVELP1)
      else
        UsR(QVELP1)    = qR(QVELP1) - qR(QMAGN)*qR(QMAGP1)*((sM - qR(QVELN))/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2))
      endif
    
! Second Perpendicular dir
      if(abs(qL(QMAGN)*qL(QMAGP2)*(sM - qL(QVELN))).le. 1d-14) then
        UsL(QVELP2) = qL(QVELP2)
      else
        UsL(QVELP2)    = qL(QVELP2) - qL(QMAGN)*qL(QMAGP2)*((sM - qL(QVELN))/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2))
      endif
      if(abs(qR(QMAGN)*qR(QMAGP2)*(sM - qR(QVELN))).le. 1d-14) then
        UsR(QVELP2) = qR(QVELP2)
      else
        UsR(QVELP2)    = qR(QVELP2) - qR(QMAGN)*qR(QMAGP2)*((sM - qR(QVELN))/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2))
      endif

      UsL(QU:QW) = UsL(QU:QW)*UsL(QRHO)
      UsR(QU:QW) = UsR(QU:QW)*UsR(QRHO)

!--------------------------------------------------------- B * states ----------------------------------------------------------------------

!Magnetic Fields

!Normal dir
      UsL(QMAGN) = qL(QMAGN)
      UsR(QMAGN) = qR(QMAGN) 

!Perpendicular dir
      if(abs(qL(QMAGP1)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)).lt.1d-14) then
        UsL(QMAGP1) = qL(QMAGP1)
      else
        UsL(QMAGP1) = qL(QMAGP1)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2)
      endif
      if(abs(qR(QMAGP1)*(qR(QRHO)*(sR - qR(QVELN))**2 - qL(QMAGN)**2)).lt.1d-14) then 
        UsR(QMAGP1) = qR(QMAGP1)
      else
        UsR(QMAGP1) = qR(QMAGP1)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2)
      endif

! Second Perpendicular dir
      if(abs(qL(QMAGP2)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)).lt. 1d-14) then
        UsL(QMAGP2) = qL(QMAGP2)
      else
        UsL(QMAGP2) = qL(QMAGP2)*(qL(QRHO)*(sL - qL(QVELN))**2 - qL(QMAGN)**2)/(qL(QRHO)*(sL - qL(QVELN))*(sL - sM) - qL(QMAGN)**2)
      endif
      if(abs(qR(QMAGP2)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)).lt.1d-14) then 
        UsR(QMAGP2) = qR(QMAGP2)
      else
        UsR(QMAGP2) = qR(QMAGP2)*(qR(QRHO)*(sR - qR(QVELN))**2 - qR(QMAGN)**2)/(qR(QRHO)*(sR - qR(QVELN))*(sR - sM) - qR(QMAGN)**2)
      endif

!Energy *Stored in Pressure slot
      UsL(QPRES) = (sL - qL(QVELN))*eL - ptL*qL(QVELN) + pst*sM + &
                    qL(QMAGN)*(dot_product(qL(QU:QW),qL(QMAGX:QMAGZ)) &
                 -  dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)))
      UsL(QPRES) = UsL(QPRES)/(sL - sM)
      UsR(QPRES) = (sR - qR(QVELN))*eR - ptR*qR(QVELN) + pst*sM + &
                    qR(QMAGN)*( dot_product(qR(QU:QW),qR(QMAGX:QMAGZ)) &
                 -  dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)))
      UsR(QPRES) = UsR(QPRES)/(sR - sM)

!speeds
      ssL = sM - abs(qL(QMAGN))/sqrt(UsL(QRHO))
      ssR = sM + abs(qR(QMAGN))/sqrt(UsR(QRHO))

!----------------------------------------- ** states ------------------------------------------------------------------------------
!Dens
      UssL(QRHO)  = UsL(QRHO)
      UssR(QRHO)  = UsR(QRHO)

!Vel in normal direction
      UssL(QVELN)    = sM
      UssR(QVELN)    = sM

!Vel in perpendicular direction
      UssL(QVELP1)    = (sqrt(UsL(QRHO))*UsL(QVELP1)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QVELP1)/UsR(QRHO)&
                      + (UsR(QMAGP1) - UsL(QMAGP1))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
      UssR(QVELP1)    = UssL(QVELP1)

!Vel in second perpendicular direction
      UssL(QVELP2)    = (sqrt(UsL(QRHO))*UsL(QVELP2)/UsL(QRHO) + sqrt(UsR(QRHO))*UsR(QVELP2)/UsR(QRHO)&
                      + (UsR(QMAGP2) - UsL(QMAGP2))*sign(1.d0,qL(QMAGN)))/(sqrt(UsL(QRHO)) + sqrt(UsR(QRHO)))
      UssR(QVELP2)    = UssL(QVELP2)

      UssL(QU:QW) = UssL(QU:QW)*UssL(QRHO)
      UssR(QU:QW) = UssR(QU:QW)*UssR(QRHO)

!B in normal direction
      UssL(QMAGN) = UsL(QMAGN)
      UssR(QMAGN) = UsR(QMAGN)

!B in perpendicular direction
      UssL(QMAGP1) = (sqrt(UsL(QRHO))*UsR(QMAGP1) + sqrt(UsR(QRHO))*UsL(QMAGP1) + sqrt(UsL(QRHO)*UsR(QRHO)) &
                   * (UsR(QVELP1)/UsR(QRHO) - UsL(QVELP1)/UsL(QRHO))*sign(1.d0,UsL(QMAGN)))/(sqrt(UsL(QRHO))&
                   + sqrt(UsR(QRHO)))
      UssR(QMAGP1) = UssL(QMAGP1)

!B in second perpendicular direction
      UssL(QMAGP2) = (sqrt(UsL(QRHO))*UsR(QMAGP2) + sqrt(UsR(QRHO))*UsL(QMAGP2) + sqrt(UsL(QRHO)*UsR(QRHO)) &
                   * (UsR(QVELP2)/UsR(QRHO) - UsL(QVELP2)/UsR(QRHO))*sign(1.d0,UsR(QMAGN)))/(sqrt(UsL(QRHO))&
                   + sqrt(UsR(QRHO)))
      UssR(QMAGP2) = UssL(QMAGP2)

!Energy **Stored in Pressure slot
      UssL(QPRES) = UsL(QPRES) - sqrt(UsL(QRHO))*(dot_product(UsL(QU:QW)/UsL(QRHO),UsL(QMAGX:QMAGZ)) &
                               - dot_product(UssL(QU:QW)/UssL(QRHO),UssL(QMAGX:QMAGZ)))*sign(1.d0, UsL(QMAGN))
      UssR(QPRES) = UsR(QPRES) + sqrt(UsR(QRHO))*(dot_product(UsR(QU:QW)/UsR(QRHO),UsR(QMAGX:QMAGZ)) &
                               - dot_product(UssR(QU:QW)/UssR(QRHO),UssR(QMAGX:QMAGZ)))*sign(1.d0, UsR(QMAGN))

!--------------------------------------------------------- Fluxes ----------------------------------------------

      FsL  = FL + sL*UsL - sL*uL
      FssL = FsL + ssL*UssL - ssL*UsL
      FsR  = FR + sR*UsR - sR*uR
      FssR = FsR + ssR*UssR - ssR*UsR

!Solve the RP
      if(sL .gt. 0.d0) then
        flx(i,j,k,:) = FL
      elseif(sL .le. 0.d0 .and. ssL .gt. 0.d0) then
        flx(i,j,k,:) = FsL
      elseif(ssL .le. 0.d0 .and. sM .gt. 0.d0) then
        flx(i,j,k,:) = FssL
      elseif(sM .le. 0.d0 .and. ssR .gt. 0.d0) then
        flx(i,j,k,:) = FssR
      elseif(ssR .le. 0.d0 .and. sR .gt. 0.d0) then
        flx(i,j,k,:) = FsR
      else 
        flx(i,j,k,:) = FR
      endif

     end do
    end do
   end do

end subroutine hlld

!================================================= Calculate the Conservative Variables ===============================================

subroutine PToC(q, u)

   use amrex_fort_module, only : rt => amrex_real
   use meth_params_module

   implicit none

   real(rt), intent(in ) ::q(QVAR)
   real(rt), intent(out) ::u(QVAR)
   u = 0.d0

   u(URHO)        = q(QRHO)
   u(UMX)         = q(QRHO)*q(QU)
   u(UMY)         = q(QRHO)*q(QV)
   u(UMZ)         = q(QRHO)*q(QW)
   u(UEINT)       = q(QPRES)/(gamma_minus_1)
   u(UEDEN)       = u(UEINT)  + 0.5d0*q(QRHO)*dot_product(q(QU:QW),q(QU:QW)) &
                  + 0.5d0*(dot_product(q(QMAGX:QMAGZ),q(QMAGX:QMAGZ)))
   u(QMAGX:QMAGZ) = q(QMAGX:QMAGZ)

end subroutine PToC

end module hlld_solver

