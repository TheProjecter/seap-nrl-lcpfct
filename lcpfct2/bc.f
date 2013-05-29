	subroutine bcfcn(icout,icin,jctop,nr,nx,ny,bc_wall,bc_outf,relax,gammam,rhoamb,preamb)
	include 'prm.h'
	integer icout,icin,nr,nx,bc_wall,bc_outf
	real relax,gammam,rhoamb,preamb
	real a(3,2),b(3,2),c(2,3),d(2,3)
 	real bc(4,0:npt),bcx(4,0:npt,0:npt),bcv(2,0:npt,0:npt,2),bce(2,0:npt,0:npt)
	integer ibc(0:npt,0:npt)
	integer ibox_l,ibox_r,icyl_lls,icyl_lrs,icyl_rls,icyl_rrs
	common /fct_bc/ bc,bcx,bcv,bce,ibc
!	density
	a(1,bc_wall)=1.0
	b(1,bc_wall)=0.0
	a(1,bc_outf)=1.0-relax
	b(1,bc_outf)=relax*rhoamb
!	velocity
	a(2,bc_wall)=1.0
	b(2,bc_wall)=0.0
	a(2,bc_outf)=1.0-relax
	b(2,bc_outf)=0.0
	c(bc_wall,:)=1.0
	c(bc_outf,:)=1.0-relax
	d(bc_wall,:)=0.0
	d(bc_outf,:)=0.0
!	energy
	a(3,bc_wall)=1.0
	b(3,bc_wall)=0.0
	a(3,bc_outf)=1.0-relax
	b(3,bc_outf)=relax*preamb/gammam


	ibox_l=   0
	ibox_r=nx+1
	icyl_lls=nr-icout
	icyl_lrs=nr-icin 
	icyl_rls=nr+icin 
	icyl_rrs=nr+icout

	bcx(1,0:nx+1,0:ny+1)=1.0
	bcx(2,0:nx+1,0:ny+1)=0.0
	bcx(3,0:nx+1,0:ny+1)=0.0
	bcx(4,0:nx+1,0:ny+1)=0.0

	bcx(1,           0 ,1:   ny+1)=a(1,bc_outf)
	bcx(2,           0 ,1:   ny+1)=b(1,bc_outf)
	bcx(1,icyl_lls,1:jctop)=a(1,bc_wall)
	bcx(2,icyl_lls,1:jctop)=b(1,bc_wall)
	bcx(1,icyl_lrs,1:jctop)=a(1,bc_wall)
	bcx(2,icyl_lrs,1:jctop)=b(1,bc_wall)
	bcx(1,icyl_rls,1:jctop)=a(1,bc_wall)
	bcx(2,icyl_rls,1:jctop)=b(1,bc_wall)
	bcx(1,icyl_rrs,1:jctop)=a(1,bc_wall)
	bcx(2,icyl_rrs,1:jctop)=b(1,bc_wall)
	bcx(1,           nx+1,1:   ny+1)=a(1,bc_outf)
	bcx(2,           nx+1,1:   ny+1)=b(1,bc_outf)

	bcx(3,              1,1:   ny)=1.0
	bcx(4,              1,1:   ny)=1.0
	bcx(3,icyl_lls,1:jctop)=1.0
	bcx(4,icyl_lls,1:jctop)=1.0
	bcx(3,icyl_lrs,1:jctop)=1.0
	bcx(4,icyl_lrs,1:jctop)=1.0
	bcx(3,icyl_rls,1:jctop)=1.0
	bcx(4,icyl_rls,1:jctop)=1.0
	bcx(3,icyl_rrs,1:jctop)=1.0
	bcx(4,icyl_rrs,1:jctop)=1.0
	bcx(3,           nx+1,1:   ny)=1.0
	bcx(4,           nx+1,1:   ny)=1.0

C	bcx(3:4,iwall_lls+1:iwall_lrs-1,1:jctop)=1.0
C	bcx(3:4,iwall_rls+1:iwall_rrs-1,1:jctop)=1.0

	bcx(1:2,icyl_lls+1:icyl_lrs-1,1:jctop)=1.0
	bcx(1:2,icyl_rls+1:icyl_rrs-1,1:jctop)=0.0
	print*,'%bcx, ',-(icout+1)+nr,-(icin -1)+nr,(icin -1)+nr,(icout+1)+nr
	print*,'%bcy, ',ny,jctop


	bcv(1,0:nx+1,0:ny+1,:)=1.0
	bcv(2,0:nx+1,0:ny+1,:)=0.0

	bcv(1,           0 ,1:   ny+1,:)=a(2,bc_outf)
	bcv(2,           0 ,1:   ny+1,:)=b(2,bc_outf)
	bcv(1,icyl_lls,1:jctop,:)=a(2,bc_wall)
	bcv(2,icyl_lls,1:jctop,:)=b(2,bc_wall)
	bcv(1,icyl_lrs,1:jctop,:)=a(2,bc_wall)
	bcv(2,icyl_lrs,1:jctop,:)=b(2,bc_wall)
	bcv(1,icyl_rls,1:jctop,:)=a(2,bc_wall)
	bcv(2,icyl_rls,1:jctop,:)=b(2,bc_wall)
	bcv(1,icyl_rrs,1:jctop,:)=a(2,bc_wall)
	bcv(2,icyl_rrs,1:jctop,:)=b(2,bc_wall)
	bcv(1,           nx+1,1:   ny+1,:)=a(2,bc_outf)
	bcv(2,           nx+1,1:   ny+1,:)=b(2,bc_outf)

	bcv(1:2,icyl_lls+1:icyl_lrs-1,1:jctop,:)=1.0
	bcv(1:2,icyl_rls+1:icyl_rrs-1,1:jctop,:)=0.0

	ibc=0

	ibc(ibox_l  ,1:ny+1 )=bc_outf
	ibc(icyl_lls,1:jctop)=bc_wall
	ibc(icyl_lrs,1:jctop)=bc_wall
	ibc(icyl_rls,1:jctop)=bc_wall
	ibc(icyl_rrs,1:jctop)=bc_wall
	ibc(ibox_r  ,1:ny+1 )=bc_outf

	bce(1,0:nx+1,0:ny+1)=1.0
	bce(2,0:nx+1,0:ny+1)=0.0

	bce(1,           0   ,1:   ny+1)=a(3,bc_outf)
	bce(2,           0   ,1:   ny+1)=b(3,bc_outf)
	bce(1,icyl_lls,1:jctop  )=a(3,bc_wall)
	bce(2,icyl_lls,1:jctop  )=b(3,bc_wall)
	bce(1,icyl_lrs,1:jctop  )=a(3,bc_wall)
	bce(2,icyl_lrs,1:jctop  )=b(3,bc_wall)
	bce(1,icyl_rls,1:jctop  )=a(3,bc_wall)
	bce(2,icyl_rls,1:jctop  )=b(3,bc_wall)
	bce(1,icyl_rrs,1:jctop  )=a(3,bc_wall)
	bce(2,icyl_rrs,1:jctop  )=b(3,bc_wall)
	bce(1,           nx+1,1:   ny+1)=a(3,bc_outf)
	bce(2,           nx+1,1:   ny+1)=b(3,bc_outf)

	bce(1:2,icyl_lls+1:icyl_lrs-1,1:jctop)=1.0
	bce(1:2,icyl_rls+1:icyl_rrs-1,1:jctop)=0.0

	return
	end
