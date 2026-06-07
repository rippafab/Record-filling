
program refill

  implicit none

  integer Ntot,NNtot
  double precision:: final_time,bin
  
  parameter(Ntot=1E4) !Total number of events
  parameter(NNtot=Ntot*100) !Total number of events
  parameter(final_time=100d0*Ntot)
  parameter(bin=0.01)
  
  double precision:: time(NNtot),time2(NNtot),trec(NNtot),dtrec(NNtot),prob,timee(NNtot),t,tt,kill_prob(0:NNtot),tblind(NNtot)
  double precision:: aver,dtmax,yy(NNtot),ave_nexp(NNtot)
    double precision:: time02(NNtot)
    real:: q02(NNtot)
    integer:: iflag0(NNtot)
    real:: q2(NNtot),nexp(NNtot),nexpt,nobst,dn,aa,q,qq,qq2(NNtot),qes(NNtot),q2s(NNtot),ran2,r,zz
  integer:: irem(NNtot),ind(NNtot),iflag(NNtot),nobs(NNtot),posiz(NNtot),ntrue(NNtot),jpos(NNtot),j2,jmin,nnn(NNtot),nmm,nloop
  integer:: i,j,n,nrec,nev,nev1,nevp,nrec1,k,nevpp,nr,loop,loop2,seed1,nnp,nobs0(NNtot)
  real:: cc,pp,max_rate,qmin,qth(NNtot),qth2(NNtot),qmax,qsup(NNtot),qparent(NNtot),qlastorig
  integer:: nm(100000),nc(0:100000),iq,iqmax,ith(NNtot),kp,nth2(100000)
  real:: nth(100000)
  integer:: pflag(NNtot),m,ploop
  
  open(70,file='catq2.dat',status='old')
  open(71,file='c3.dat',status='old')
  open(100,file='cat_filled.dat')!,status='new')

  seed1=-1701881
  qmin=huge(qmin)
  nobs=0
  nobs0=0
  iflag0=0
  kill_prob=0d0
  pflag=0

  do i=1,NNtot		!3000000
     read(70,*,end=99)t,q  !Legge dal file con la q trasformata
     q2(i)=q
     q02(i)=q
     qmin=min(q,qmin)
     time2(i)=t
     time02(i)=t
     iflag0(i)=0
     ind(i)=0
     read(71,*)tt,qq,k  !Legge dal file incompleta in cui nella terza colonna c'è la posizione giusta
     qq2(i)=qq
     timee(i)=tt
     posiz(i)=k
   enddo
  99 nev=i-1
  
   print*,'number of events in the catalog',nev
   do i=1,nev
      qes(i)=q2(i)
      q2s(i)=qq2(i)
   enddo
   call hpsortq(nev,qes)
   call hpsortq(nev,q2s)
   
   kill_prob(0)=0d0

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! SEARCH for RECORDS
  nrec=0
  nobst=0
  nexpt=0
  max_rate=10000.8 !THIS FIXES AN UPPER BOUND FOR THE RATE CAN BE ALSO SET A POSTERIORI IT CAN BE FIXED A POSTERIORI
  dtmax=3600*24*10  !MAximum age of a record
  nrec=0
  do i=1,nev
     do j=i+1,nev
        if(time2(j)-time2(i)>dtmax)exit
        if(q2(j)>q2(i))then
           nrec=nrec+1
           
           trec(nrec)=time2(i)
           qparent(nrec)=q2(i)
           tblind(i)=time2(j)!record blind time 
           dtrec(nrec)=time2(j)-time2(i)+1 !record age 
           ntrue(nrec)=posiz(j)-posiz(i) !true age
           nobs0(nrec)=j-i
!!$           rate(nrec)=exp(1*q2(i))/dtrec(nrec)
!!$           if(rate(nrec).gt.max_rate)then
!!$              nrec=nrec-1
!!$           else
!!$              iflag(i)=1
!!$           endif
           yy(nrec)=exp(-q2(i))
           !Nexp(nrec)=1E3*10!           rate(nrec)=exp(1*q2(i))/dtrec(nrec)
           if(nobs0(nrec)/dtrec(nrec).gt.max_rate)then
              nrec=nrec-1
           else
              iflag0(i)=1
           endif
           
           exit
        endif
     enddo
  enddo
  print*,'number of records=',nrec
  nevpp=0
  nmm=3*nev
  do nloop=1,100
     !pause
     nmm=nmm*.85

     do i=1,nrec
        nobs(i)=nobs0(i)
     enddo
     
     do i=1,nev
        q2(i)=q02(i)
        time2(i)=time02(i)
        iflag(i)=iflag0(i)
        ind(i)=0
        kill_prob(0)=-1
     enddo
     jmin=1
     max_rate=0
     do i=1,nrec
        k=1
        !        aver=1d0/rate(i)
        !nmm=int(6*nexp(i)+1)
        !if(nexp(i)>100)nmm=100000
        zz=1d0-yy(i)
        aver=dtrec(i)/(1d0-zz**nmm)/((1d0-zz**(nmm+1))/yy(i)-(nmm+1)*zz**nmm)
        do j=i-1,jmin,-1
           if(trec(j)<trec(i)-dtmax)then
              jmin=j+1
              exit
           endif
           if(trec(j)+dtrec(j)>trec(i))then
              k=k+1
              !nmm=int(6*nexp(j)+1)
              zz=1d0-yy(j)
              !rate(j)=((1d0-zz**(nn+1))-(nn+1)*zz**nn)/(yy(j)*(1d0-zz**nn))
              !rate(j)=rate(j)/dtrec(j)
              !aver=aver+1d0/rate(j)
              aver=aver+dtrec(j)*(1d0-zz**nmm)/((1d0-zz**(nmm+1))/yy(j)-(nmm+1)*zz**nmm)
           endif
        enddo
        do j=i+1,nrec
           if(trec(j)>trec(i)+dtrec(i))exit
           k=k+1
           !           aver=aver+1d0/rate(j)
           !nmm=int(6*nexp(j)+1)
           zz=1d0-yy(j)
           aver=aver+dtrec(j)*(1d0-zz**nmm)/((1d0-zz**(nmm+1))/yy(j)-(nmm+1)*zz**nmm)
        enddo
        max_rate=max(max_rate,dtrec(i)*k*1d0/aver)
        
        nexp(i)=max(0.,dtrec(i)*k*1d0/aver-0.99d0)
        
        !!nexp(i)=(1/0.7)*nexp(i)**(1./1.05)
        !     nexp(i)=ntrue(i)
        !     ave_nexp(ntrue(i))=ave_nexp(ntrue(i))+nexp(i)
        !     nnn(ntrue(i))=nnn(ntrue(i))+1
        !if(nexp(i)>0.and.k>0)write(100+nloop,*)ntrue(i),nexp(i),dtrec(i),k,time2(i),i
     
!     if(nobs(nrec).ge.0)nobst=nobst+ntrue(nrec)
!     if(nobs(nrec).ge.0)nexpt=nexpt+nexp(nrec) 
     enddo
!     print*,nexp(1),nexp(88),nexp(112)
     !pause
!  do i=1,1000
!     if(nnn(i)>0)write(15,*)i,ave_nexp(i)/nnn(i)
!  enddo
!  call flush(100+nloop)
!  print*,max_rate

!stop
!  max_rate=1E8
!$$  pause

  !Filling/Unfilling procedure
  nevp=nev
  aa=0.05
  do loop=1,1000
     n=0
     nobst=0
     do i=1,nrec
        dn=(nexp(i)-nobs(i))

!!$        
        if(nobs(i).gt.0)then
           kill_prob(i)=-dn*aa/nobs(i)*100
        else
           kill_prob(i)=0d0
        endif
        nobst=nobst+dn
        do j=1,int(dn*aa)
           n=n+1
           time2(nevp+n)=trec(i)+dtrec(i)*ran2(seed1)
           tblind(nevp+n)=time2(nevp+n)
           iflag(nevp+n)=-1
           ind(nevp+n)=i
           q2(nevp+n)=-1000
        enddo
!                do i=1,nrec
!        enddo

     enddo

        
     nev1=nevp+n
 !    print*,nev1,ind(10332)
     
!     print*,'prima',nevp,'dopo=',nev1,'aggiunti',n,nobst
     CALL HPSORT(nev1, time2,q2,ind,iflag,tblind)


        

!!!!     !!Unfilling ###############################
     k=0
     do i=1,nev1
 !       if(loop==2.and.i==10332)print*,'ppp',i,kill_prob(ind(i)),ind(i),q2(i)
        if(kill_prob(ind(i)).gt.rand())then
           irem(i)=1
           
 !          if(loop==2.and.i==10332)then
  !            write(17+nloop,*)i,nr,kill_prob(ind(i))!nexp(i),nobs(i)
  !            call flush(17+nloop)
  !         endif
        else
           irem(i)=0
        endif
     enddo
     
        
     nr=0
     do i=1,nev1
        if(irem(i).eq.1)then
           nr=nr+1
           
        else
           q2(i-nr)=q2(i)
           time2(i-nr)=time2(i)
           tblind(i-nr)=tblind(i)
           iflag(i-nr)=iflag(i)
           ind(i-nr)=ind(i)
           irem(i-nr)=irem(i)
        endif
  

     enddo
        nevp=nev1-nr
!        write(*,*)'tolti=',nr,nevp,'loop=',loop
        nnp=0
        

!###########################################################  
!!!!!Search for the new value of nobs
        nrec1=0
        do i=1,nevp
           if(iflag(i)==1)then
              do j=i+1,nevp
                 if(q2(j)>q2(i))then
                    nrec1=nrec1+1
                    nobs(nrec1)=j-i
                    exit
                 endif
              enddo
           endif
        enddo
     enddo



!!$     do i=1,nrec
!!$        write(289,*)nobs(i),nexp(i),ntrue(i),trec(i),trec(i)+dtrec(i)
!!$     enddo
!!$     call flush(89)
!!$     write(289,*)''

     do i=1,nevp
!        if(ind(i)==0)write(188,*)time2(i),tblind(i)
     enddo
!     if(kbest==1)exit
     if(nevp<nevpp*.99)then
        exit
     endif
     nevpp=nevp
     print*,'number of events loop',nloop, 'is',nevpp,nmm,nevpp*.99
  enddo
  
     qlastorig=0.
     do i=1,nevp
        if(ind(i)==0)then
           qsup(i)=0.
           qlastorig=q2(i)
        else
           qsup(i)=qlastorig
        endif
     enddo
     
     do i=1,nevp
        if(ind(i)==0)then
           !qsup(i)=max(qsup(i),qmin+2*bin) !Upper bound of magnitude for records
        endif
     enddo
     
     print*,'qmin=',qmin
     do iq=1,100000
        nm(iq)=0
        nth(iq)=0
        nth2(iq)=0
     enddo
     k=0
     qmax=-1
     do i=1,nevp
        if(ind(i)==0)then
           k=k+1
           iq=int((q2(i)-qmin)/bin)+1
           nm(iq)=nm(iq)+1   !Distribution of the incomplete catalog
           qmax=max(q2(i),qmax) !maximumm magnitude of original catalog
        endif
     enddo
     print*,qmin,qmax,k

!     do ploop=1,1
        print*,nevp,'prima'
        do i=1,nevp*10
           r=ran2(seed1)
           q=qmin-log(1.-r*(1-exp(-qmax)))
           iq=int((q-qmin)/bin)+1
           nth(iq)=nth(iq)+1./10    !theoric distribution of all events
        enddo
        iqmax=int((qmax-qmin)/bin)+1
!     do iq=1,iqmax
!        if(nm(iq)>0)write(11,*)iq*.1,nm(iq)*1./bin
!     enddo
!     write(11,*)
!     do iq=1,iqmax
!        if(nth(iq)>0)write(11,*)iq*.1,nth(iq)*1./bin
!     enddo
!     write(11,*)
        nc(0)=0
        do iq=1,iqmax
           nc(iq)=nc(iq-1)+max(0.,nth(iq)-nm(iq))  !Cumulative distribution of the difference between theoric and incomplete
        enddo
        do i=1,nevp
           if(ind(i)==0)cycle
33         r=ran2(seed1)*nc(iqmax)
           do iq=1,iqmax
              if(nc(iq)>=r)then
                 q2(i)=qmin+(iq-1)*bin+bin*ran2(seed1)  !Add magnitude according to the difference between theoric and incomplete
                 exit
              endif
           enddo
        enddo
        do i=1,nevp
           q2s(i)=q2(i)
           posiz(i)=i
        enddo
        call hpsort2(nevp,q2s,posiz)
        do j=1,nevp
           jpos(posiz(j))=j
           pflag(j)=0
        enddo
        m=0

        do while(m<nevp)
           i=ran2(seed1)*nevp+1
           !print*,m,nevp
           if(pflag(i)>0)cycle
           pflag(i)=1
           m=m+1
           if(ind(i)==0)cycle
           if(q2(i)>qsup(i))then !Added magnitude is larger than record!!! Faccio un search guardando alle magnitudo più piccolo (sono ordinate per grandezza) la prima che trovo che è minore del suo qsup(i) la scambio
              j=jpos(i)
              do k=j-1,1,-1
                 if(q2s(k)<=qsup(i))then
                    j2=posiz(k)
                    if(ind(j2)==0)cycle
                    if(q2(i)<=qsup(j2))then
                       if(q2(j2)<=qsup(i))then

                          qq=q2(i)
                          q2(i)=q2(j2)  !EXCHANGE
                          q2(j2)=qq     !EXCHANGE
                          pflag(i)=2
                          goto 35
                       endif
                    endif
                 endif
              enddo
!           q2(i)=qsup(i)
35            continue
           endif
        enddo
!!$        j=0
!!$        do i=1,nevp
!!$           if(ind(i)==0)then
!!$              q2(i-j)=q2(i)
!!$              ind(i-j)=ind(i)
!!$              time2(i-j)=time2(i)
!!$              goto 45
!!$           endif
!!$           if(q2(i)<=qsup(i))then
!!$              q2(i-j)=q2(i)
!!$              qsup(i-j)=qsup(i)
!!$              time2(i-j)=time2(i)
!!$              goto 45
!!$           endif
!!$        
!!$!           if(ran2(seed1)<.95)then
!!$!              j=j+1
!!$!           else
!!$!              q2(i-j)=q2(i)
!!$!              qsup(i-j)=qsup(i)
!!$!              time2(i-j)=time2(i)
!!$!           endif
!!$45         continue
!!$        enddo
 !       nevp=nevp-j
 !       print*,ploop,nevp,'zio'
 !    enddo
           
           
        
        do i=1,nevp
           if(ind(i)==0.or.q2(i)<=qsup(i))then
              write(100,*)time2(i),q2(i),qsup(i),ind(i),pflag(i)
!           else
!              q2(i)=qsup(i)!+ran2(seed1)*(qsup(i)-qmin)
!              write(100,*)time2(i),q2(i),qsup(i),ind(i),pflag(i)
           endif
!        if(ind(i)==0.or.q2(i)<qsup(i))then
           
!        endif
     enddo
     call flush(100)

     

  stop
  end program



SUBROUTINE HPSORT(N, R1, R2, R3, R4,R5)
  IMPLICIT NONE
  INTEGER :: N, L, IR, I, J,R3(N), R4(N)
  REAL :: R2(N)
  REAL*8 :: R1(N), R5(N)!, R3(N), R4(N),t1(N)
  REAL*8 :: RR1, RR2, RR3, RR4, RR5



  
  L = N / 2 + 1
  IR = N

10 CONTINUE
  IF (L > 1) THEN
     L = L - 1
     RR1 = R1(L)
     RR2 = R2(L)
     RR3 = R3(L)
     RR4 = R4(L)
     RR5 = R5(L)
  ELSE
     RR1 = R1(IR)
     RR2 = R2(IR)
     RR3 = R3(IR)
     RR4 = R4(IR)
     RR5 = R5(IR)

     R1(IR) = R1(1)
     R2(IR) = R2(1)
     R3(IR) = R3(1)
     R4(IR) = R4(1)
     R5(IR) = R5(1)

     IR = IR - 1
     IF (IR == 1) THEN
        R1(1) = RR1
        R2(1) = RR2
        R3(1) = RR3
        R4(1) = RR4
        R5(1) = RR5
        RETURN
     END IF
  END IF

  I = L
  J = L + L

20 IF (J <= IR) THEN
     IF (J < IR) THEN
        IF (R1(J) < R1(J + 1)) J = J + 1
     END IF
     IF (RR1 < R1(J)) THEN
        R1(I) = R1(J)
        R2(I) = R2(J)
        R3(I) = R3(J)
        R4(I) = R4(J)
        R5(I) = R5(J)
        I = J
        J = J + J
     ELSE
        J = IR + 1
     END IF
     GOTO 20
  END IF

  R1(I) = RR1
  R2(I) = RR2
  R3(I) = RR3
  R4(I) = RR4
  R5(I) = RR5
  GOTO 10

END SUBROUTINE HPSORT

SUBROUTINE HPSORT2(N, R1,R3)
  IMPLICIT NONE
  INTEGER :: N, L, IR, I, J,R3(N)
  REAL :: R1(N)
  REAL*8 :: RR1, RR3



  
  L = N / 2 + 1
  IR = N

10 CONTINUE
  IF (L > 1) THEN
     L = L - 1
     RR1 = R1(L)
     RR3 = R3(L)
  ELSE
     RR1 = R1(IR)
     RR3 = R3(IR)

     R1(IR) = R1(1)
     R3(IR) = R3(1)

     IR = IR - 1
     IF (IR == 1) THEN
        R1(1) = RR1
        R3(1) = RR3
        RETURN
     END IF
  END IF

  I = L
  J = L + L

20 IF (J <= IR) THEN
     IF (J < IR) THEN
        IF (R1(J) < R1(J + 1)) J = J + 1
     END IF
     IF (RR1 < R1(J)) THEN
        R1(I) = R1(J)
        R3(I) = R3(J)
        I = J
        J = J + J
     ELSE
        J = IR + 1
     END IF
     GOTO 20
  END IF

  R1(I) = RR1
  R3(I) = RR3
  GOTO 10

END SUBROUTINE HPSORT2


  

!*****************************************************
!*  Sorts an array RA of length N in ascending order *
!*                by the Heapsort method             *
!* ------------------------------------------------- *
!* INPUTS:                                           *
!*	    N	  size of table RA                       *
!*      RA	  table to be sorted                     *
!* OUTPUT:                                           *
!*	    RA    table sorted in ascending order        *
!*                                                   *
!* NOTE: The Heapsort method is a N Log2 N routine,  *
!*       and can be used for very large arrays.      *
!*****************************************************         
       SUBROUTINE HPSORTQ(N,RA)
	real RA(N)
	L=N/2+1
	IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR 
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
 10	continue
	if(L > 1)then
	   L=L-1
	   RRA=RA(L)
        else
	   RRA=RA(IR)
           RA(IR)=RA(1)
	   
	   IR=IR-1
	   if(IR.eq.1)then
	      RA(1)=RRA
	      
	      return
	   end if
	end if
	I=L
	J=L+L
 20	if(J.le.IR)then
	   if(J < IR)then
	      if(RA(J) < RA(J+1))  J=J+1
	   end if
	   if(RRA < RA(J))then
	      RA(I)=RA(J)
	      
	      I=J; J=J+J
	   else
	      J=IR+1
	   end if
	   goto 20
	end if
	RA(I)=RRA
	
	goto 10
	END


 
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,IA1=40014,IA2=40692)
      PARAMETER (IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!___________________________________________________________________________
