all :
	make clean
	make f
	make c

c:
	make fast2dc
	make gmonc.out
	make gprofc

f :
	make fast2df
	make gmonf.out
	make gproff

fast2df : fast2d.o gasdyn.o gasdynp.o fctp.o bc.o blk.o fct.o makegrid.o residiff.o sources.o timer.o velocity.o zerodiff.o zeroflux.o
	pgf90 -pg -Mextend -r8 -o fast2df fast2d.o gasdyn.o gasdynp.o fctp.o bc.o blk.o fct.o makegrid.o residiff.o sources.o timer.o velocity.o zerodiff.o zeroflux.o
	
fast2dc : fast2d.o gasdyn.o gasdynp.o lcpfctc.o bc.o blk.o fct.o makegrid.o residiff.o sources.o timer.o velocity.o zerodiff.o zeroflux.o
	pgf90 -pg -Mextend -r8 -o fast2dc fast2d.o gasdyn.o gasdynp.o lcpfctc.o bc.o blk.o fct.o makegrid.o residiff.o sources.o timer.o velocity.o zerodiff.o zeroflux.o

fast2d.o : fast2d.f prm.h
	pgf90 -pg -Mextend -r8 -c fast2d.f

gasdyn.o : gasdyn.f prm.h
	pgf90 -pg -Mextend -r8 -c gasdyn.f

gasdynp.o : gasdyn.f prm.h
	pgf90 -pg -Mextend -r8 -c gasdynp.f

lcpfctc.o : lcpfctc.c
	gcc -pg -c lcpfctc.c

fctp.o : fctp.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c fctp.f

bc.o : bc.f prm.h
	pgf90 -pg -Mextend -r8 -c bc.f

blk.o : blk.f
	pgf90 -pg -Mextend -r8 -c blk.f

fct.o : fct.f fct.h
	pgf90 -pg -Mextend -r8 -c fct.f

makegrid.o : makegrid.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c makegrid.f

residiff.o : residiff.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c residiff.f

sources.o : sources.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c sources.f

timer.o : timer.f90 timer.h
	pgf90 -pg -Mextend -r8 -c timer.f90

velocity.o : velocity.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c velocity.f

zerodiff.o : zerodiff.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c zerodiff.f

zeroflux.o : zeroflux.f prm.h fct.h
	pgf90 -pg -Mextend -r8 -c zeroflux.f

clean :
	-rm fast2df fast2dc *.o gmon*.out gprof*.out fast2d*.dat profile* fast2d*.out

gmonc.out : fast2dc
	fast2dc >> fast2dc.out
	mv gmon.out gmonc.out
	mv fast2d.dat fast2dc.dat
	mv fast2d.0000.dat fast2dc.0000.dat
	mv profile.dat.0000 profilec.dat.0000

gmonf.out : fast2df
	fast2df >> fast2df.out
	mv gmon.out gmonf.out
	mv fast2d.dat fast2df.dat
	mv fast2d.0000.dat fast2df.0000.dat
	mv profile.dat.0000 profilef.dat.0000

gprofc : fast2dc gmonc.out
	gprof fast2dc gmonc.out > gprofc.out

gproff : fast2df gmonf.out
	gprof fast2df gmonf.out > gproff.out

.PNONY : gmon.out
