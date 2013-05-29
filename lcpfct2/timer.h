	real, dimension(100) :: tprofiler,eprofiler,cprofiler
	character(40), dimension(100) :: aprofiler
	character(10) atime
	character( 8) adate
	character(40) amsg
	integer, dimension(100) :: ncall
	integer, dimension(100) :: iprofiler,jprofiler,level,dprofiler
	integer 		:: ielapsed,npr,iunit,icount,jcount,iprof
	common /tau_cmn/ tprofiler,eprofiler,cprofiler,			   &
	&		 ncall,iprofiler,jprofiler,level,ilevel,	   &
	&		 dprofiler,ic0,jc0,nprof,irate,ielapsed,npr,iunit, &
	&		 icount,jcount,iprof,				   &
	&		 aprofiler, adate,atime
!$omp threadprivate (/tau_cmn/)
