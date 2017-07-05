module hlld_solver

implicit none

private hlldx,hlldy,hlldz
public hlld


subroutine hlld(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3, &
				dir)

implicit none 

	integer, intent(in)   :: qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3
	integer, intent(in)   :: flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3
	integer, intent(in)   :: dir

	real(rt), intent(in)  :: qm(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
	real(rt), intent(in)  :: qp(qpd_l1:qpd_h1,qpd_l2:qpd_h2,qpd_l3:qpd_h3,QVAR)
	real(rt), intent(out) :: flx(flx_l1:flx_h1,flx_l2:flx_h2,flx_l3:flx_h3,NVAR)

if(dir.eq.1) then
	call hlldx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)
elseif(dir.eq.2) then
	call hlldy(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)
else 
	call hlldz(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)
endif
end subroutine hlld

subroutine hlldx(qm,qp,qpd_l1,qpd_l2,qpd_l3,qpd_h1,qpd_h2,qpd_h3, &
                 flx,flx_l1,flx_l2,flx_l3,flx_h1,flx_h2,flx_h3)

end subroutine hlldx
end module hlld_solver
