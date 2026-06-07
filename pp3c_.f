c************************************************************************
c    Programma per la simulazione del modello metas
c    La probabilità P(t,m) si calcola dalla relazione
c    P(t,m)=\sum_i alpha/{Exp[gamma(t-t_i)/10^(b(m_i-m))]-1+epsilon}+rate
c 
c*************************************************************************


	integer seed1,ZBQLPOI
	real*8 AAA,AAAf
	real*8 final_time,muaa
	real q(10000000),q2(10000000),bin0
        real*8 prob,probTOT(0:100000),probref,tt,t000,tt0
	real*8 time,cc,t00,itime(100000000),itime2(10000000),dt
	integer nrem(-100:1000,0:100)
	real*8 aveq,aveqp,aveqpp
	integer ave_nt(0:1000),ave_ntp(0:1000),ave_ntpp(0:1000)
	real*8 ave_b(0:1000),ave_bp(0:1000),ave_bpp(0:1000)
	real*8 ave_b2(0:1000),ave_bp2(0:1000),ave_bpp2(0:1000)
	integer nrealb(0:1000),irealiz,nrealbp(0:1000),nrealbpp(0:1000)
	integer ndt(0:2000),ndt2(0:2000),nt(10000000)
	real*8 bev,bev2,stev,lambda(10000000),ndt3(0:1000)
	real xq(10000000),phi,qe(1000000),dqp(1000000)
	integer nx(0:100000),ix,nqp(0:100000),posiz(10000000)
	real*8 itimee(10000000),rs(10000000)
	real qes(1000000),xqs(1000000),pdixc(0:10000000)
	
c*********************************************************
!	open(56,file='catq.dat')
	open(57,file='cat2.dat')
	open(58,file='cat.dat')
	


	!read(*,*)seed1
c	nmain=1200*10*1.9 !*10*realiz
	!n\\\\main=1400000*0.0005!15 per etasi !0.2
        nmain=1
	qmin0=0
	qsup=7.5
	final_time=4*30*3600*24.
	tau1=102.00000000001
	alpha2=0.75!/2.5!2.5/2.5*0.85
	tau2=20  !2.2!10**(-alpha*3.5*1.)*3600*24*0.000232!0.001
	stdv=.0001!sqrt(4*0.3)
	pp=1.2!1.08		!omori exponent
	alpha=0.9
	AAA=0.06*2*.5*1
	cc=400.47100!0.024*3600*24 !0.01
	bb=1!2.3!/log(10.)!0.97

	br=AAA*(1-10**((alpha-bb)*(qsup-qmin)))/((bb-alpha)*log(10.))
	br=br*bb*log(10.)/(1-10**((-bb)*(qsup-qmin))) !Normaliz. della P(m)
	print*,'branching ratio=',br
	bin0=2

	print*,final_time!,int(dlog(final_time-0.0001)/bin0+1)

	do loop=0,1000
	   ndt(loop)=0
	   ndt2(loop)=0
	   ndt3(loop)=0d0
	enddo
	
	
	nrealiz=1

	
	do irealiz=1,nrealiz
	   seed1=-7141231+irealiz	
	   do i=1,1000000
	      q(i)=-100
	   enddo
	   qmin=qmin0
	!print*,loop


	
	   nevent=1
	   q(1)=7.15
	   itime(1)=0.00001
		   
	   do j=1,1!  single loop
	      if(q(j).le.-10)goto 887
	      aa=AAA*10**(alpha*(q(j)-qmin))

	      muaa=AAA*10**(alpha*(q(j)-qmin))
	      nn=ZBQLPOI(muaa)

	      
	      tt0=itime(j)
	      qq0=q(j)
!	      nn=0.0003*final_time
	      
	      nafter=0
	      nevent=0
	      do i=1,nn
		 time=cc*((ran2(seed1))**(1d0/(1d0-pp)))-cc
!		 time=1*ran2(seed1)*final_time		 
		 if(tt0+time.le.final_time)then
		    nevent=nevent+1
		    itime(nevent)=tt0+time
		    
 34		    rr=ran2(seed1)
		    qq=qmin-(1./bb)*log(1-rr)
		    q(nevent)=qq !min(qq,qsup)

		    
		    
		    bb1=-0.5*bb
		    bb2=bb*2.5
		    xc=4.5

		    s1=sign(1.,bb1)
		    xm=30
		    a=exp(bb1*xc)*(bb1+bb2)/bb2
		    zk=s1/(-1+a-exp((bb2+bb1)*xc)*bb1*exp(-bb2*xm)/bb2)
		    rsc=zk*s1*(exp(bb1*xc)-1)

c$$$		    if(rr.le.rsc)then
c$$$		       q(nevent)=log(1+rr*s1/zk)/bb1
c$$$		    else
c$$$		       q(nevent)=(-log((a-1)*s1-rr/zk)
c$$$     &		    +(bb2+bb1)*xc+log(s1*bb1/bb2))/bb2
c$$$		    endif
!		    print*,q(nevent)
		    idt=int(dlog(1d0*itime(nevent))/0.1)+1
		    ndt(idt)=ndt(idt)+1
	         endif
	      enddo
	      
 79	   enddo

 887	   continue
	   print*,'nevent=',nevent
	   call hpsort(nevent,itime,q)

	   zeta=3.
	
	   do ix=1,1000
	      nx(ix)=0
	   enddo




	   
	   do i=1,nevent
	  
	      xq(i)=q(i)
	      write(58,*)itime(i),q(i)

	      
	      if(itime(i).gt.tau1)then
!		 ix=int(log(xq(i))/0.1)+1 !OR
		 ix=int((xq(i))/0.1)+1 !OR
		 nx(ix)=nx(ix)+1
	      endif
	   enddo
c$$$

c$$$	   do ix=1,1000
c$$$	      if(nx(ix).gt.0)write(15,*)
c$$$     &  	exp(ix*0.1),!*(10**(1.1*(-qmag(1)+qmin))),
c$$$     &	       nx(ix)*exp(-ix*0.1)/(1*(exp(0.1)-1))
c$$$           enddo
	   do ix=1,1000
	      if(nx(ix).gt.0)write(15,*)(ix-1)*0.1,nx(ix)/0.1
           enddo
	   
	   call flush(15)

	   do i=1,nevent
	      write(87,*)itime(i),i
	      write(187,*)itime(i),q(i)
	   enddo


	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	   
!!      THINNING 

	   
	   nevent2=nevent
	   do i=1,nevent
	      q2(i)=q(i)
	      itime2(i)=itime(i)
	      posiz(i)=i
	   enddo

	   print*,nevent,nevent2

	
	
	   do ix=1,1000
	      nx(ix)=0
	   enddo
	   do i=1,nevent2
	      xq(i)=q2(i)
	      write(57,*)itime2(i),q2(i),posiz(i)

	      
	      if(itime2(i).gt.tau1)then
	         ix=int((xq(i))/0.1)+1 !OR
	         nx(ix)=nx(ix)+1
	      endif
	   enddo
c$$$

	   do ix=1,1000
	      if(nx(ix).gt.0)write(16,*)(ix-1)*.1,nx(ix)/0.1
 !    &   	exp(ix*0.1),!*(10**(1.1*(-qmag(1)+qmin))),
 !    &        nx(ix)*exp(-ix*0.1)/(1*(exp(0.1)-1))
	   enddo

	   call flush(57)
	stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    SEARCH FOR BEST FIT	

	do i=1,nevent2
	   xqs(i)=q2(i)
	enddo
	call hpsortq(nevent2,itime2,xqs)
	bbp=1

	


	   
	do i=1,nevent2
	   rs(i)=i*1d0/nevent2!ncount*1d0/nevent2
	enddo
	


	xm=5*log(1.*nevent2)
	bb2=1
	sum_phim=1E12
	do loopx=1,60
	   xc=xm/(1.1**loopx)
	   
	   do loopb=-9,9
	      bb1=loopb*0.1-0.001
	      !call hpsortq(nevent2,itime2,q2)
	      
	      s1=sign(1.,bb1)

	      a=exp(bb1*xc)*(bb1+bb2)/bb2
	      zk=s1/(-1+a-exp((bb2+bb1)*xc)*bb1*exp(-bb2*xm)/bb2)
	      rsc=zk*s1*(exp(bb1*xc)-1)
	
	   
	
	      do i=1,nevent2-1
		 if(rs(i).le.rsc)then
		    qe(i)=log(1+rs(i)*s1/zk)/bb1
		 else
		    qe(i)=(-log((a-1)*s1-rs(i)/zk)
     &		 +(bb2+bb1)*xc+log(s1*bb1/bb2))/bb2
		 endif
		 itimee(i)=itime2(i)
		 ix=int(qe(i)/0.1)+1
		 nx(ix)=nx(ix)+1
	      enddo   
	      qe(nevent2)=qe(nevent2-1)+.1
	      itimee(nevent2)=itime2(nevent2-1)

	      do ix=1,1000
!		 if(nx(ix).ge.1)write(26,*)(ix-1)*0.1,nx(ix)
	      enddo

	      call hpsort(nevent2,itimee,qe)
	



	      avedq=0
	      ndq=0
	      dqth=0
	      do i=1,nevent2-100
		 do j=1,100
		    dq=qe(i+j)-qe(i)
		    if(dq.ge.0)then
		       if(dq.ge.dqth)then
			  avedq=avedq+(dq-dqth)
			  ndq=ndq+1
			  dqp(ndq)=dq
		       endif
		       exit
		    endif
		 enddo
	      enddo
	      bbp=ndq*1d0/avedq

	      call hpsortdq(ndq,dqp)

	      ncount=0
	      sum_phi=0.
	      do i=1,ndq-1
		 ncount=ncount+1
		 phi=ncount*1./ndq
		 z=-log(1-phi)/bbp
		 sum_phi=sum_phi+(z-dqp(i))**2
	      enddo
	      if(sum_phi.le.sum_phim)then
		 bbb=bb1
		 xcc=xc
		 print*,bbb,xcc,sum_phi/(ndq-1),bbp
		 sum_phim=sum_phi
	      endif
	   enddo
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	xm=5*log(1.*nevent2)
	bb2=1
	xc=xcc
	bb1=bbb
!	call hpsortq(nevent2,itime2,q2)
	      
	s1=sign(1.,bb1)

	a=exp(bb1*xc)*(bb1+bb2)/bb2
	zk=s1/(-1+a-exp((bb2+bb1)*xc)*bb1*exp(-bb2*xm)/bb2)
	rsc=zk*s1*(exp(bb1*xc)-1)
	print*,zk,bb1
	
	   
	do ix=1,1000
	   nx(ix)=0
	enddo


	
	do i=1,nevent2-1
	   !ncount=ncount+1
	   !rr=ncount*1d0/nevent2
	   if(rs(i).le.rsc)then
!	      print*,i,rs(i),1	   
	      qe(i)=log(1+rs(i)*s1/zk)/bb1
	      qes(i)=qe(i)  !111
!	      print*,1,rs(i),rsc,xc,qe(i)
	   else
!	      print*,i,rs(i),3
	      qe(i)=(-log((a-1)*s1-rs(i)/zk)
     &		 +(bb2+bb1)*xc+log(s1*bb1/bb2))/bb2
	      qes(i)=qe(i)
!	      itimee(i)=itime2(i)
	   endif
	   ix=int(qe(i)/0.1)+1
!       print*,ix
	   nx(ix)=nx(ix)+1
	enddo
	qe(nevent2)=qe(nevent2-1)+.1
	qes(nevent2)=qe(nevent2)
	      
	do ix=1,1000
	   if(nx(ix).ge.1)write(26,*)(ix-1)*0.1,nx(ix)
	enddo
	
	
	call hpsort(nevent2,itime2,qe)


	do ix=1,1000
	   nqp(ix)=0
	enddo 


	avedq=0
	ndq=0
	dqth=0
	do i=1,nevent2-100
!	   write(67,*)itime2(i),qe(i)
	   do j=1,100
	      dq=qe(i+j)-qe(i)
	      if(dq.ge.0)then
		 if(dq.ge.dqth)then
		    avedq=avedq+(dq-dqth)
		    ndq=ndq+1
		    dqp(ndq)=dq
		    idq=int(dq/0.1)+1
		    nqp(idq)=nqp(idq)+1
		 endif
		 exit
	      endif
	   enddo
	enddo
	bbp=ndq*1d0/avedq
	print*,'bbp3=',bbp,ndq
	call hpsortdq(ndq,dqp)
	do i=1,nevent2
	   write(56,*)itime2(i),bbp*qe(i),posiz(i)
	enddo
	call flush(56)
	call flush(57)
	call flush(58)
	
	ncount=0
	sum_phi=0.
	do i=1,ndq-1
	   ncount=ncount+1
	   phi=ncount*1./ndq
	   z=-log(1-phi)/bbp
	   write(66,*)dqp(i),z
	   sum_phi=sum_phi+(z-dqp(i))**2
	enddo
	write(*,*)'sum_phi=',bb1,xc,sum_phi/(ndq-1)
	
	
	do ix=1,1000
	   if(nqp(ix).ge.1)write(65,*)(ix-1)*0.1,2*nqp(ix)
	enddo 
	
	
	do i=1,nevent2
	   nt(i)=0
	   lambda(i)=-1
	enddo


	dm0=.0001
	dm1=1
	dm2=dm0
	!bbp=1
	
	do i=1,nevent2
!       lambda(i)=1/(itime2(i+1)-itime2(i)+0.000001)
!       flag=0
	   qmax=qe(i)
	   qth=qmax+dm0
	   do j=i+1,nevent2
	      
!	      if(q2(j).gt.qth)then
!       do k=i,j
	      if(qe(j).gt.qth)then
 		 pp=exp(bbp*(qth))
		 
		 lambda(i)=	!lambda(i)+
!     &		    (j-i)/(itime2(j)-itime2(i)+0.0000001)
     &		 pp/(1*(itime2(j)-itime2(i)+0.00001))
!		    write(66,*)ind(j)-ind(i),exp(bb*q2(i))		    
		 exit
	      endif
!	      qmax=max(qmax,q2(j))
!	      qth=max(qth,qmax+dm0)
	   enddo
	enddo


	do i=2,nevent2
!	   if(nt(i).ge.1)write(33,*)itime2(i),lambda(i)/nt(i)
	   idt=int(1d0*dlog(itime2(i))/0.1)+1
	   if(lambda(i).gt.0)then
	      ndt3(idt)=ndt3(idt)+1d0/(lambda(i)) !/nt(i))
	      ndt2(idt)=ndt2(idt)+1
	   endif
	enddo
!       call flush(33)
	end do
	
	
	do i=1,1000
	   if(ndt(i).ge.1)then
	   write(12,*)exp(i*0.1),
     &	   ndt(i)*exp(-(i-1)*0.1)/(exp(0.1)-1)/nrealiz
	endif
	   enddo
	   call flush(12)
	do i=1,1000
	   if(ndt2(i).ge.1)
     &	write(22,*)exp(i*0.1)
     &	   ,ndt2(i)*exp(-(i-1)*0.1)/(exp(0.1)-1)/nrealiz
	enddo
	call flush(22)
	znev=0
	do i=1,1000
	   if(ndt3(i).ge.0.0001)then
	      write(32,*)exp(i*0.1),1d0/(ndt3(i)/ndt2(i))
	      znev=znev+(exp(0.1)-1)*exp(i*0.1)*ndt2(i)/ndt3(i)
	   endif
	enddo
	call flush(32)
	print*,'refilled=',znev,'true=',nevent

	do i=51,nevent2
	   pdix=bbp*exp(-bbp*qes(i))*(qes(i)-qes(i-50))
	   pdix=pdix/(xqs(i)-xqs(i-50))
	   write(41,*)xqs(i),pdix*znev!(qes(i+100)-qes(i))/(qs(i+100)-qs(i))
	enddo



	stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	   

	
	end





c__________________________________________________________________________________________
	FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
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
c__________________________________________________________________________________
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
       SUBROUTINE HPSORTQ(N,q,RA)
	real RA(N)
	real*8 q(N)
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
	   qRA=q(L)
	else
	   RRA=RA(IR)
	   qRA=q(IR)
	   RA(IR)=RA(1)
	   q(IR)=q(1)
	   IR=IR-1
	   if(IR.eq.1)then
	      RA(1)=RRA
	      q(1)=qRA
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
	      q(I)=q(J)
	      I=J; J=J+J
	   else
	      J=IR+1
	   end if
	   goto 20
	end if
	RA(I)=RRA
	q(I)=qRA
	goto 10
	END
c___________________________________________________-__-------------------------

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
       SUBROUTINE HPSORT(N,RA,q)
	real*8 RA(N)
	real q(N)
	L=N/2+1
	IR=N
!The index L will be decremented from its initial value during the
!"hiring" (heap creation) phase. Once it reaches 1, the index IR 
!will be decremented from its initial value down to 1 during the
!"retirement-and-promotion" (heap selection) phase.
 110	continue
	if(L > 1)then
	   L=L-1
	   RRA=RA(L)
	   qRA=q(L)
	else
	   RRA=RA(IR)
	   qRA=q(IR)
	   RA(IR)=RA(1)
	   q(IR)=q(1)
	   IR=IR-1
	   if(IR.eq.1)then
	      RA(1)=RRA
	      q(1)=qRA
	      return
	   end if
	end if
	I=L
	J=L+L
 120	if(J.le.IR)then
	   if(J < IR)then
	      if(RA(J) < RA(J+1))  J=J+1
	   end if
	   if(RRA < RA(J))then
	      RA(I)=RA(J)
	      q(I)=q(J)
	      I=J; J=J+J
	   else
	      J=IR+1
	   end if
	   goto 120
	end if
	RA(I)=RRA
	q(I)=qRA
	goto 110
	END
c___________________________________________________-__-----------------------
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
       SUBROUTINE HPSORTDQ(N,RA)
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
c___________________________________________________-__-------------------------
