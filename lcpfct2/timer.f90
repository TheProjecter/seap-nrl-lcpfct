	subroutine tau_profile_init()
	include 'timer.h'
	nprof=0
	tprofiler=0.0
	eprofiler=0.0
	cprofiler=0.0
	aprofiler=' '
	ncall=0
	ilevel=0
	jelapsed=0
	dprofiler=0
	call date_and_time(adate,atime)
	call system_clock(ic0,irate,imax)
	return
	end

	subroutine tau_profile_timer(profiler,a)
	integer profiler(2)
	include 'timer.h'
	character*(*) a
	nprof=nprof+1
	 profiler(    2)=nprof
	aprofiler(nprof)=a
	return
	end

	subroutine tau_profile_start(profiler)
	integer profiler(2),ierr
	include 'timer.h'
	call system_clock(icount,irate,imax)
	profiler(1)=icount
	iprof=profiler(2)
	ncall(iprof)=ncall(iprof)+1
	iprofiler(iprof)=icount
	level(iprof)=ilevel+1
	if (dprofiler(iprof).eq.0) ilevel=ilevel+1
	return
	end

	subroutine tau_profile_stop (profiler)
	integer profiler(2),ierr
	include 'timer.h'
	icount=profiler(1)
	iprof=profiler(2)
	call system_clock(jcount,irate,imax)
	jprofiler(iprof)=jcount
	itime=jcount-icount
	tprofiler(iprof)=tprofiler(iprof)+itime
	if (dprofiler(iprof).eq.0) ilevel=ilevel-1
	return
	end

	subroutine tau_profile_exit(bmsg)
	integer profiler(2)
	character( 8) bdate
	character(10) btime
	character(*)  bmsg

	include 'timer.h'
	character(4) anpr
	lmsg=len(bmsg)
	lmsg=min(lmsg,40)
	amsg=' '
	amsg(1:lmsg)=bmsg(1:lmsg)
	call date_and_time(bdate,btime)
	write(anpr,'(i4.4)') npr
	iunit=50+npr
	open (iunit,file='profile.dat.'//anpr,status='unknown',form='unformatted')
	write(iunit) adate,atime,bdate,btime,amsg
	write(iunit) nprof
	write(iunit) level(1:nprof)
	write(iunit) ncall(1:nprof)
	write(iunit) iprofiler(1:nprof)
	write(iunit) jprofiler(1:nprof)
	write(iunit) tprofiler(1:nprof)
	write(iunit) aprofiler(1:nprof)
	write(iunit) dprofiler(1:nprof)
	close(iunit)

	return
	end	

	subroutine tau_profile_set_node(inum)
	integer inum
	include 'timer.h'
	npr=inum
	return
	end

	subroutine tau_profile_timer_dynamic(profiler,a)
	integer profiler(2)
	include 'timer.h'
	character*(*) a
	nprof=nprof+1
	 profiler(    2)=nprof
	aprofiler(nprof)=a
	dprofiler(nprof)=1
	return
	end

	subroutine tau_profile_fexit(bmsg)
	integer profiler(2)
	character( 8) bdate
	character(10) btime
	character(*)  bmsg

	include 'timer.h'
	character(4) anpr
	lmsg=len(bmsg)
	lmsg=min(lmsg,40)
	amsg=' '
	amsg(1:lmsg)=bmsg(1:lmsg)
	call date_and_time(bdate,btime)
	write(anpr,'(i4.4)') npr
	iunit=50+npr
	open (iunit,file='profile.txt.'//anpr,status='unknown',form='formatted')
	write(iunit,*) adate,atime,bdate,btime,amsg
	write(iunit,*) nprof
	write(iunit,*) level(1:nprof)
	write(iunit,*) ncall(1:nprof)
	write(iunit,*) iprofiler(1:nprof)
	write(iunit,*) jprofiler(1:nprof)
	write(iunit,*) tprofiler(1:nprof)
	write(iunit,*) aprofiler(1:nprof)
	write(iunit,*) dprofiler(1:nprof)
	close(iunit)

	return
	end	
