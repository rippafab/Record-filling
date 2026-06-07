program legge_q_fit_exp
  implicit none

  integer :: i, nevent2, ios
  real(8), allocatable :: q2(:),qe(:)
  real(16), allocatable :: itime2(:)
  integer, allocatable :: rank(:),ind1(:),ind2(:) 
  real :: qq
  real(16) :: time,avedq,sum_phi,sum_phim,phi,z
  real(8) :: r, kp, abest, kbest, ap, y,dq,dqth,dqp(1000000)
  real(8) :: bbp,bbpm
  integer :: j,nrec,kloop,aloop,loop
  
  ! 1. Pre-lettura per contare gli eventi o allocazione larga
  ! In questo esempio usiamo una dimensione fissa capiente o contiamo i record
  allocate(q2(1000000), itime2(1000000)) 

  open(71, file='cat2.dat', status='old')
  open(199, file='catq2.dat', status='new')
  
  i = 1
  do while (i <= 1000000)
      read(71, *, iostat=ios) time, qq
      if (ios /= 0) exit ! Fine file o errore
      q2(i) = dble(qq)
      itime2(i) = time
      i = i + 1
  end do
  nevent2 = i - 1
  close(71)
  allocate(qe(nevent2))
  
  print*, "Eventi letti:", nevent2

  ! CERCA RECORDS
  nrec=0
  dqth=0
  do i=1,nevent2-100
     do j=1,100
        dq=q2(i+j)-q2(i)
        if(dq.ge.0)then
           if(dq.ge.dqth)nrec=nrec+1
           exit
        endif
     enddo
  enddo
  allocate(ind1(nrec))
  allocate(ind2(nrec))
  nrec=0
  do i=1,nevent2-100
     do j=1,100
        dq=q2(i+j)-q2(i)
        if(dq.ge.0)then
           if(dq.ge.dqth)then
              nrec=nrec+1
              ind1(nrec)=i
              ind2(nrec)=i+j
           endif
           exit
        endif
     enddo
  enddo

  
  if (nevent2 <= 0) stop "Errore: file vuoto o non trovato."

  ! 2. Ri-allocazione precisa per il calcolo del rank
  allocate(rank(nevent2))
  
  call compute_rank(q2(1:nevent2), nevent2, rank)
    
  ! 3. Parametri per il metodo di Newton

  
!  k = 1!0.85931997139781713!92.d0
!  a = 0.01
  
  sum_phim = huge(0.0)
  do loop=1,100
     kp=0.01+0.1d0*loop
     do aloop=0,50
        ap=aloop*0.02
        !ap=a
        !     print*,aloop
!11      ap=abest*(1d0+0.1*sign(1.,rand()-0.5))
!         if(ap<0.or.ap>1)goto 11
         !     do kloop=1,1
         !        kp=kloop*0.1
         !        kp=k
         call convert_to_y(nevent2,rank,kp,ap,qe)
         call diff_to_exp(nevent2,qe,nrec,ind1,ind2,sum_phi,bbp)
         !print*,sum_phi,ap
         if(sum_phi<sum_phim)then
            abest=ap
            kbest=kp
            sum_phim=sum_phi
            bbpm=bbp
           ! print*,qe(14760)
            print*,sum_phi,'a=',ap,'k=',kp,'bbp=',bbp,aloop
         endif
         
      enddo
   enddo
      

      call convert_to_y(nevent2,rank,kbest,abest,qe)
      !print*,qe(14760)
      call diff_to_exp(nevent2,qe,nrec,ind1,ind2,sum_phi,bbp)
      !call update_best(sum_phi,be,xe,bbest,xbest)
      print*,'Best=',sum_phi,bbp,abest,kbest
      do i=1,nevent2
         write(199,*)itime2(i),bbp*qe(i),q2(i)
      enddo
      call flush(199)
  !call update_best(sum_phi,be,xe,bbest,xbest)

  
!  print*, "Risultato Newton (y):", y

contains

   subroutine compute_rank(q_in, n, rank_out)
     implicit none
     integer, intent(in) :: n
     real(8), intent(in) :: q_in(n)
     integer, intent(out) :: rank_out(n)
     
     real(8), allocatable :: val(:)
     integer, allocatable :: idx(:)
     integer :: i, curr_r

     allocate(val(n), idx(n))

     do i = 1, n
        val(i) = q_in(i)
        idx(i) = i
     end do

     call sort_with_index(val, idx, n)

     curr_r = 1
     rank_out(idx(1)) = curr_r

     do i = 2, n
        if (val(i) > val(i-1)) then
           curr_r = i
        end if
        rank_out(idx(i)) = curr_r
     end do

     deallocate(val, idx)
   end subroutine compute_rank

   subroutine sort_with_index(a, idx, n)
     integer, intent(in) :: n
     real(8), intent(inout) :: a(n)
     integer, intent(inout) :: idx(n)
     call quicksort(a, idx, 1, n)
   end subroutine sort_with_index

   recursive subroutine quicksort(a, idx, left, right)
     real(8), intent(inout) :: a(:)
     integer, intent(inout) :: idx(:)
     integer, intent(in) :: left, right
     integer :: i, j, itmp
     real(8) :: pivot, tmp

     i = left
     j = right
     pivot = a((left + right)/2)

     do
        do while (a(i) < pivot); i = i + 1; end do
        do while (a(j) > pivot); j = j - 1; end do
        if (i <= j) then
           tmp = a(i); a(i) = a(j); a(j) = tmp
           itmp = idx(i); idx(i) = idx(j); idx(j) = itmp
           i = i + 1
           j = j - 1
        end if
        if (i > j) exit
     end do

     if (left < j) call quicksort(a, idx, left, j)
     if (i < right) call quicksort(a, idx, i, right)
   end subroutine quicksort

    subroutine diff_to_exp(n,qe,nrec,ind1,ind2,sum_phi,bbp)
    implicit none
    integer, intent(in) :: n,nrec
    real(8), intent(in) :: qe(n)
    integer, intent(in) :: ind1(n),ind2(n)
    real(16), intent(out) :: sum_phi
    real(8) :: bbp
    real(16) :: dq,avedq,phi,z
    integer :: i

    avedq=0
    do i=1,nrec
       dq=qe(ind2(i))-qe(ind1(i))
       avedq=avedq+dq
       dqp(i)=dq
!       write(88,*)i,dq
    enddo
    bbp=nrec*1d0/avedq
!    call flush(88)
!      print*,bbp,ndq
!		 pause
    call hpsortdq(nrec,dqp)
		 
    sum_phi=0d0
    do i=1,nrec-1
       !phi = (i - 0.5d0)/nrec
       phi=i*1d0/nrec
       z=-log(1-phi)/bbp
       sum_phi=sum_phi+(z-dqp(i))**2
    enddo
!    print*,'bbp=',bbp,sum_phi
  end subroutine diff_to_exp
  
   
  subroutine convert_to_y(n,rank,kin,ain,qe_out)
     implicit none
     integer, intent(in) :: n
     integer, intent(in) :: rank(n)
     real(8), intent(in)  :: kin, ain
     real(8), intent(out) :: qe_out(n)
          
     real(8) :: AA, BB, NN, r
     integer :: i

     
!     AA = a + (1.d0 - a)/(1.d0 - exp(-k*xf)) 
!     BB = (1.d0 - a)/((k+1.d0)*(1.d0 - exp(-k*xf))) 
!     NN = a + ((1.d0-a)/(1.d0-exp(-k*xf)))*(k/(k+1.d0)) 
     AA=1.d0
     BB= (1.d0 - ain)/(kin+1.d0)
     NN= (kin+ain)/(kin+1.d0)
     
     do i=1,n
        !r=(rank(i) - 0.5d0)/n
        r=rank(i)*1d0/n-1d0/(4*n)  !maximum y=log(4*N)
!        call solve_y_newton(r, kin,AA, BB, NN, y)
!        qe(i)=y
!        print*,y,r,1
        call solve_y_bis(r, kin,AA, BB, NN, y)
        qe(i)=y
!        print*,y,2
!        pause
        

        !        write(199,*)itime2(i),y
     enddo
     !     call flush(199)
   end subroutine convert_to_y
   subroutine solve_y_newton(r,k,AA, BB, NN, y_res)
     real(8), intent(in)  :: r,k,AA, BB, NN
     real(8), intent(out) :: y_res
     real(8) :: CC, f, df, ynew
     integer :: iter
     real(8), parameter :: tol = 1.0d-12
     
!     AA = a + (1.d0 - a)/(1.d0 - exp(-k*xf)) 
!     BB = (1.d0 - a)/((k+1.d0)*(1.d0 - exp(-k*xf))) 
!     NN = a + ((1.d0-a)/(1.d0-exp(-k*xf)))*(k/(k+1.d0)) 
     CC = AA - BB - NN*r
     
     y_res = -log(1.d0 - r + 1.d-12) 
     do iter = 1, 50
        f  = AA*exp(-y_res) - BB*exp(-(k+1.d0)*y_res) - CC
        df = -AA*exp(-y_res) + (k+1.d0)*BB*exp(-(k+1.d0)*y_res)
        if (abs(df) < 1.d-14) exit
        ynew = y_res - f/df
        if (ynew < 0.d0) ynew = 0.5d0*y_res 
        if (abs(ynew - y_res) < tol) then
           y_res = ynew
           return
        end if
        y_res = ynew
     end do
     
   end subroutine solve_y_newton
    subroutine solve_y_bis(r,k, AA, BB, NN, y_res)
     real(8), intent(in)  :: r,k, AA, BB, NN
     real(8), intent(out) :: y_res
     real(8) :: CC, f, df, ynew
     integer :: iter
     real(8), parameter :: tol = 1.0d-12
     real(8) :: uL, uR, um
     real(8) :: fL, fR, fm
     
!     AA = a + (1.d0 - a)/(1.d0 - exp(-k*xf)) 
!     BB = (1.d0 - a)/((k+1.d0)*(1.d0 - exp(-k*xf))) 
!     NN = a + ((1.d0-a)/(1.d0-exp(-k*xf)))*(k/(k+1.d0)) 
     CC = AA - BB - NN*r
     

   ! --- dominio per u = exp(-y)
     uL = 0.d0
     uR = 1.d0

   ! --- funzione agli estremi
     fL = AA*uL - BB*uL**(k+1.d0) - CC   ! = -C
     fR = AA*uR - BB*uR**(k+1.d0) - CC   ! = A - B - C

   ! --- controllo di esistenza soluzione
     if (fL*fR > 0.d0) then
        print *, 'Errore: soluzione non bracketed'
        stop
     end if
     
     do iter = 1, 1000
        um = 0.5d0*(uL + uR)
        fm = AA*um - BB*um**(k+1.d0) - CC

        if (abs(fm) < tol) exit

        if (fL*fm > 0.d0) then
           uL = um
           fL = fm
        else
           uR = um
           fR = fm
        end if
     end do

   ! --- ritorno a y
     if (um <= 0.d0) then
        y_res = huge(1.d0)
     else
        y_res = -log(um)
     end if

   end subroutine solve_y_bis

     
   

   SUBROUTINE HPSORTDQ(N,RA)
     integer N,L,IR
     real*8 RA(N),RRA
     L=N/2+1
     IR=N
     !The index L will be decremented from its initial value during the
     !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
     !will be decremented from its initial value down to 1 during the
     !"retirement-and-promotion" (heap selection) phase.
10   continue
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
20   if(J.le.IR)then
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
   END SUBROUTINE HPSORTDQ

   
   
 end program legge_q_fit_exp
