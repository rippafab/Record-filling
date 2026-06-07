program legge_cat
  
  parameter(NN=1E7)
  real*8 time(NN),q(NN)
  integer n(-100:1000)
  
  !open(10,file='cat2.dat',status='old')
  open(10,file='catq2.dat',status='old')
  !open(10,file='fort.199',status='old')
  !open(10,file='cat_filled.dat',status='old')
  !open(10,file='c4.dat',status='old')

  kk=1000
  bin=0.1
  nnp=0
  do i=1,NN
     read(10,*,end=99)time(i),q(i)
     
!     if(time(i)>2E5.and.time(i)<5*2E5)write(100,*)time(i),0.01
!     if(time(i)>14*2E5.and.time(i)<(5+13)*2E5)write(101,*)time(i),0.01
!     if(time(i)>(2+2*13)*2E5.and.time(i)<(5+2*13)*2E5)write(101,*)time(i),0.01
!     if(time(i)>(2+3*13)*2E5.and.time(i)<(5+3*13)*2E5)write(101,*)time(i),0.01
!     if(time(i)>(8+3*13)*2E5.and.time(i)<(11+3*13)*2E5)write(101,*)time(i),0.01
     
     
     !if(time(i).ge.2753380.and.time(i).le.2753391)then
!    if (time(i).ge.2752478.and.time(i).le.2758091)then
!        nnp=nnp+1
        !write(74,*)time(i),nnp
        !              write(899,*)i,time2(i),iflag(i),irem(i),loop
!     endif
  enddo
        !        call flush(899)

        !pause

99 nev=i-1
  print*,nev,nnp
  
  do ix=-100,1000
     n(ix)=0
  enddo
  
  do i=1,nev
     ix=int(q(i)/bin)+1
     n(ix)=n(ix)+1
   enddo

  do ix=-100,1000
     if(n(ix)>0)write(201,*)ix*bin,n(ix)*1d0/(bin)!*nev)
  enddo

  stop
  end program
  
