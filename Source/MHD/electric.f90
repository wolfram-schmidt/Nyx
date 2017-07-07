subroutine electric(Q, E) !Use ideal Ohm's Law
use amrex_fort_module, only : rt => amrex_real
use meth_mhd_params_module, only : QVAR, QU,QV, QW, QMAGX, QMAGY, QMAGZ

implicit none

 real(rt), intent(in)	::Q(QVAR)
 real(rt), intent(out) 	::E(3)

		!E = -v X B
	E(1)	= -Q(QV)*Q(QMAGZ) + Q(QW)*Q(QMAGY)
	E(2)	= -Q(QW)*Q(QMAGX) + Q(QU)*Q(QMAGZ)
	E(3)	= -Q(QU)*Q(QMAGY) + Q(QV)*Q(QMAGX)
 
end subroutine electric
