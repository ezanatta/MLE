      program ffinder

      parameter(m=10000)

      double precision xb(m),yb(m),sb(m),fb(m)
      double precision xd(m),yd(m),sd(m),fd(m),ftot(m),f(m)
      double precision sub(m),areb(m),sut(m),aret(m),f2(m)
      double precision xf(m),yf(m),sf(m),ff(m),suf(m),aref(m),idf(m)
      integer idb(m),idd(m)


c      open(1,file='N2768B_star.txt',status='old')
c      open(2,file='N2768T_star.txt',status='old')
      open(4,file='N2768F_GC.txt',status='old')
      open(3,file='N2768_f.GC.dat',status='unknown')

      do i=1,m
         n=n+1

         

c         read(1,*,end=11)xb(i),yb(i),sb(i),fb(i),sub(i),areb(i),idb(i)
c         read(2,*,end=22)xd(i),yd(i),sd(i),ftot(i),sut(i),aret(i),idd(i)
         read(4,*,end=44)xf(i),yf(i),sf(i),ff(i),suf(i),aref(i),idf(i)

c         ftot(i)=fb(i)+fd(i)

         f(i)=fb(i)/ftot(i)
         f2(i)=ff(i)/aref(i)

         write(3,*)f(i),f2(i),idb(i)

         enddo

 11      close(1)
 22      close(2)
         close(3)
 44      close(4)

         stop
         end



