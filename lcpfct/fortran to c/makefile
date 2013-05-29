all :
	make clean
	make fast2d
	make gmon.out
	make gprof

fast2d : fast2d.o gasdyn.o lcpfctc.o lcpfctf.o lcpfct_others.o fctblk_init.o
	g77 -pg -o fast2dc fast2d.o gasdyn.o lcpfctc.o lcpfct_others.o fctblk_init.o
	g77 -pg -o fast2df fast2d.o gasdyn.o lcpfctf.o lcpfct_others.o fctblk_init.o

fast2d.o : fast2d.f
	g77 -pg -c fast2d.f

gasdyn.o : gasdyn.f
	g77 -pg -c gasdyn.f

lcpfctc.o : lcpfctc.c
	gcc -pg -c lcpfctc.c

lcpfctf.o : lcpfctf.f
	g77 -pg -c lcpfctf.f

lcpfct_others.o : lcpfct_others.f
	g77 -pg -c lcpfct_others.f

fctblk_init.o : fctblk_init.f
	g77 -pg -c fctblk_init.f

clean :
	-rm fast2dc fast2df fast2d.o gasdyn.o lcpfctc.o lcpfctf.o lcpfct_others.o fctblk_init.o gmonc.out gmonf.out

gmon.out : fast2dc fast2df
	fast2dc
	mv gmon.out gmonc.out
	fast2df
	mv gmon.out gmonf.out

gprof : fast2dc fast2df gmonc.out gmonf.out
	gprof fast2dc gmonc.out > gprofc.out
	gprof fast2df gmonf.out > gproff.out

.PNONY : gmon.out