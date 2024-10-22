program ti_fep_bar
implicit none

  integer*4 :: i,j,k,m,ist
  integer*4 :: nwindows 
  integer*4 :: frames
  integer*4 :: ti
  integer*4 :: fep
  integer*4 :: bar     
  character(len=160) :: mdoutfile
  character(len=160) :: tioutfile,fepoutfile,baroutfile
  character(len=160) :: filename
  character(len=160) :: fline
  character(len=160) :: infile 
  real*8, allocatable :: dvdl(:,:)
  real*8, allocatable :: energy(:,:,:)
  real*8 :: T
  real*8 :: RT
  real*8 :: cut_ti
  real*8 :: cut_fep
  real*8 :: cut_bar
 
  namelist /ti_fep_bar_cntrl/ nwindows,frames,mdoutfile,tioutfile,fepoutfile,baroutfile,T,cut_ti,cut_fep,cut_bar,ti,fep,bar
  T=300.0       !default: 300.0 
  cut_ti=3.0    !default: 3.0
  cut_fep=3.0   !default: 3.0
  cut_bar=999.0 !default: 999.0 
  ti=1
  fep=1
  bar=1

  i=iargc()
  do j=1,i
    call getarg(j,fline)
    if(fline.eq."-i")then
      call getarg(j+1,infile)
    endif
  enddo
 
  open(10,file=infile)
    read(10,nml=ti_fep_bar_cntrl)
  close(10)
   
  RT=8.314*T*0.001/4.1868
  allocate(dvdl(nwindows,frames*2))
  allocate(energy(nwindows,frames,nwindows)) 

  write(*,*)
  write(*,'(A53)')'----------' 
  write(*,'(A53)')'TI-FEP-BAR'
  write(*,'(A53)')'----------'
  write(*,*)
  write(*,'(A20)')'parameter settings'
  write(*,'(A20,I12,A20,F12.1)')'nwondow:',nwindows,'Temperature:',T
  write(*,'(A20,F12.1,A20,F12.1,A20,F12.1)')'cut_ti:',cut_ti, 'cut_fep:',cut_fep,'cut_bar:',cut_bar
  write(*,'(A20,I12,A20,I12,A20,I12)')'ti:',ti, 'fep:',fep,'bar:',bar
  write(*,*)
  write(*,'(A20)')'read mdoutfile ...'
  open(11,file=mdoutfile) 
    do i=1,nwindows
      j=0
      m=0
      read(11,'(A160)')filename
      open(12,file=trim(adjustl(filename)))
        do while(.true.)
        read(12,'(A79)',iostat=ist)fline
        if(ist/=0)exit
        if(fline(7:31).eq.'A V E R A G E S   O V E R')exit
          if(fline(2:6).eq.'DV/DL')then
            m=m+1 
            read(fline(10:24),'(F16.4)')dvdl(i,m)
          endif
          if(fline(1:4).eq.'MBAR')then
            j=j+1
            do k=1,nwindows
              read(12,'(A80)')fline
              if(fline(20:20).eq.'*')then
                energy(i,j,k)=99999999
              else
                read(fline(19:31),'(F16.3)')energy(i,j,k)
              endif
            enddo
          endif
        enddo
      close(12)
      write(*,'(A20,I12,A20,F12.1,A20,F12.1)')'window:',i,'Frame(TI):',real(m)/2.0,'Frame(FEP/BAR):',real(j)
    enddo
  close(11)

  if(ti.eq.1)then
    write(*,*)
    write(*,'(A20)')'calculate TI ...'
    call tical(dvdl,cut_ti,frames,nwindows,tioutfile)
  endif

  if(fep.eq.1)then
    write(*,*)
    write(*,'(A20)')'calculate FEP ...'
    call fepcal(RT,energy,cut_fep,frames,nwindows,fepoutfile)
  endif

  if(bar.eq.1)then
    write(*,*)
    write(*,'(A20)')'calculate BAR ...'
    call barcal(RT,energy,cut_bar,frames,nwindows,baroutfile)
  endif

  write(*,'(A3)')'End'
end program ti_fep_bar

subroutine tical(dvdl,cut_ti,frames,nwindows,tioutfile)
implicit none

  integer*4 :: i,m,counter
  integer*4 :: nwindows
  integer*4 :: frames
  integer*4 :: user(frames)
  real*8 :: cut_ti
  real*8 :: sum_ti
  real*8 :: sem_ti
  real*8 :: sum_squ
  real*8 :: sum_dvdl
  real*8 :: sem(nwindows)
  real*8 :: std(nwindows)
  real*8 :: ave_dvdl(nwindows)
  real*8 :: dvdl(nwindows,frames*2)
  character(len=160) :: tioutfile

  open(21,file=tioutfile)
  write(21,'(4A16)')"Nwindows","dvdl","SEM_dvdl","(kcal/mol)"
  sum_ti=0.0
  sem_ti=0.0
  do i=1,nwindows
    sum_dvdl=0.0
    sum_squ=0.0
    do m=1,frames
      sum_dvdl=sum_dvdl+dvdl(i,2*m)
    enddo
    ave_dvdl(i)=sum_dvdl/frames
    do m=1,frames
      sum_squ=sum_squ+(dvdl(i,2*m)-ave_dvdl(i))**2
    enddo
    std(i)=sqrt(sum_squ/frames)
    counter=0
    do m=1,frames
      if(dvdl(i,2*m).gt.(ave_dvdl(i)-cut_ti*std(i)) .and. dvdl(i,2*m).lt.(ave_dvdl(i)+cut_ti*std(i)))then
        user(m)=1
        counter=counter+1
      else
        user(m)=0
      endif
    enddo
    sum_dvdl=0.0
    do m=1,frames
      if(user(m).eq.1)then
        sum_dvdl=sum_dvdl+dvdl(i,2*m)
      endif
    enddo
    ave_dvdl(i)=sum_dvdl/counter
    sum_ti=sum_ti+ave_dvdl(i)
    sum_squ=0.0
    do m=1,frames
      if(user(m).eq.1)then
        sum_squ=sum_squ+(dvdl(i,2*m)-ave_dvdl(i))**2
      endif
    enddo
    sem(i)=sqrt(sum_squ)/counter
    sem_ti=sem_ti+sem(i)**2
    write(21,'(I16,2F16.4)')i,ave_dvdl(i),sem(i)
    write(*,'(A20,I12,A20,F12.1)')'window:',i,'Frame(TI):',real(counter)
  enddo
  write(21,'(A16,2F16.4)')"Sum TI:",(sum_ti-ave_dvdl(1)/2-ave_dvdl(nwindows)/2)/(nwindows-1),&
                                    sqrt(sem_ti-(3.0/4.0)*(sem(1)**2)-(3.0/4.0)*(sem(nwindows)**2))/(nwindows-1)
  close(21)
end subroutine tical

subroutine fepcal(RT,energy,cut_fep,frames,nwindows,fepoutfile)
implicit none

  integer*4 :: i,j,k,counter 
  integer*4 :: nwindows !the number of windows in MD simulation
  integer*4 :: frames
  integer*4 :: user(frames)
  real*8 :: RT
  real*8 :: cut_fep
  real*8 :: energy(nwindows,frames,nwindows)
  real*8 :: delta_u(frames)
  real*8 :: std(nwindows)
  real*8 :: sum_delta_u
  real*8 :: ave_delta_u
  real*8 :: delta_A(nwindows-1)
  real*8 :: sem_delta_A(nwindows-1)
  real*8 :: sum_squ
  real*8 :: sum_exp_u
  real*8 :: ave_exp_u
  real*8 :: sum_delta_A
  real*8 :: sum_sem_delta_A
  character(len=160) :: fepoutfile

  open(21,file=fepoutfile)
  write(21,'(4A16)')"Nwindows","Delta_A","SEM_A","(kcal/mol)"
  sum_delta_A=0.0
  sum_sem_delta_A=0.0
  do i=1,nwindows-1
    sum_delta_u=0.0
    do j=1,frames
      delta_u(j)=energy(i,j,i+1)-energy(i,j,i)
      sum_delta_u=sum_delta_u+delta_u(j)
    enddo
    ave_delta_u=sum_delta_u/frames
    sum_squ=0.0
    do j=1,frames
      sum_squ=sum_squ+(delta_u(j)-ave_delta_u)**2
    enddo
    std(i)=sqrt(sum_squ/frames)
    counter=0
    do j=1,frames
      if(delta_u(j).gt.(ave_delta_u-cut_fep*std(i)) .and. delta_u(j).lt.(ave_delta_u+cut_fep*std(i)))then
        user(j)=1
        counter=counter+1
      else
        user(j)=0
      endif
    enddo 
    sum_exp_u=0.0
    do j=1,frames
      if(user(j).eq.1)then
        sum_exp_u=sum_exp_u+exp(-delta_u(j)/RT)
      endif
    enddo
    ave_exp_u=sum_exp_u/counter
    sum_squ=0.0
    do j=1,frames
      if(user(j).eq.1)then
        sum_squ=sum_squ+(exp(-delta_u(j)/RT)-ave_exp_u)**2
      endif
    enddo
    delta_A(i)=-RT*log(ave_exp_u)
    sem_delta_A(i)=abs((RT/ave_exp_u)*sqrt(sum_squ)/frames)
    sum_delta_A=sum_delta_A+delta_A(i)
    sum_sem_delta_A=sum_sem_delta_A+(sem_delta_A(i))**2
    write(21,'(I16,2F16.4)')i,delta_A(i),sem_delta_A(i)
    write(*,'(A20,I12,A20,F12.1)')'window:',i,'Frame(FEP):',real(counter)
  enddo
  write(21,'(A16,2F16.4)')"Sum FEP:",sum_delta_A,sqrt(sum_sem_delta_A)
  close(21)
end subroutine fepcal 

subroutine barcal(RT,energy,cut_bar,frames,nwindows,baroutfile)
implicit none

   integer*4 :: i,j,counter
   integer*4 :: nwindows
   integer*4 :: frames
   integer*4 :: user(frames)
   real*8 :: cut_bar
   real*8 :: RT
   real*8 :: C
   real*8 :: sum_f0_1
   real*8 :: sum_f1_0
   real*8 :: sem_f0_1
   real*8 :: sem_f1_0
   real*8 :: sum_squ_f0_1
   real*8 :: sum_squ_f1_0
   real*8 :: sum_sem_f0_1
   real*8 :: sum_sem_f1_0
   real*8 :: sum_delta_A
   real*8 :: sum_sem_delta_A
   real*8 :: std_f0_1(nwindows)
   real*8 :: std_f1_0(nwindows)
   real*8 :: delta_A(nwindows)
   real*8 :: sem_delta_A(nwindows)
   real*8 :: u0_1(frames)
   real*8 :: u1_0(frames)
   real*8 :: f0_1(frames)
   real*8 :: f1_0(frames)
   real*8 :: ave_f0_1(nwindows)
   real*8 :: ave_f1_0(nwindows)
   real*8 :: energy(nwindows,frames,nwindows)
   character(len=160) :: baroutfile

   open(21,file=baroutfile)
   write(21,'(4A16)')"Nwindows","Delta_A","SEM_A","(kcal/mol)"
   sum_delta_A=0.0
   sum_sem_delta_A=0.0
   do i=1,nwindows-1
     do j=1,frames
       u0_1(j)=energy(i+1,j,i)-energy(i+1,j,i+1)
       u1_0(j)=energy(i,j,i+1)-energy(i,j,i)
     enddo
     C=0.0
     do while(.true.)
       sum_f0_1=0.0
       sum_f1_0=0.0
       do j=1,frames
         f0_1(j)=1.0/(1.0+exp((u0_1(j)+C)/RT))
         f1_0(j)=1.0/(1.0+exp((u1_0(j)-C)/RT))
         sum_f0_1=sum_f0_1+f0_1(j)
         sum_f1_0=sum_f1_0+f1_0(j)
       enddo
       ave_f0_1(i)=sum_f0_1/frames
       ave_f1_0(i)=sum_f1_0/frames
       sum_squ_f0_1=0.0
       sum_squ_f1_0=0.0
       do j=1,frames
         sum_squ_f0_1=sum_squ_f0_1+(f0_1(j)-ave_f0_1(i))**2
         sum_squ_f1_0=sum_squ_f1_0+(f1_0(j)-ave_f1_0(i))**2
       enddo
       std_f0_1(i)=sqrt(sum_squ_f0_1/frames)
       std_f1_0(i)=sqrt(sum_squ_f1_0/frames)
       counter=0
       do j=1,frames
         if(f0_1(j).gt.(ave_f0_1(i)-cut_bar*std_f0_1(i)) .and. f0_1(j).lt.(ave_f0_1(i)+cut_bar*std_f0_1(i))&
      .and. f1_0(j).gt.(ave_f1_0(i)-cut_bar*std_f1_0(i)) .and. f1_0(j).lt.(ave_f1_0(i)+cut_bar*std_f1_0(i)))then
           user(j)=1
           counter=counter+1
         else
           user(j)=0
         endif
       enddo
       sum_f0_1=0.0
       sum_f1_0=0.0
       do j=1,frames
         if(user(j).eq.1)then
           sum_f0_1=sum_f0_1+f0_1(j)
           sum_f1_0=sum_f1_0+f1_0(j)
         endif
       enddo
       ave_f0_1(i)=sum_f0_1/counter
       ave_f1_0(i)=sum_f1_0/counter
       sum_sem_f0_1=0.0
       sum_sem_f1_0=0.0
       do j=1,frames
         if(user(j).eq.1)then
           sum_sem_f0_1=sum_sem_f0_1+(f0_1(j)-ave_f0_1(i))**2
           sum_sem_f1_0=sum_sem_f1_0+(f1_0(j)-ave_f1_0(i))**2
         endif
       enddo
       sem_f0_1=sqrt(sum_sem_f0_1)/counter
       sem_f1_0=sqrt(sum_sem_f0_1)/counter
       sem_delta_A(i)=sqrt((sem_f0_1*RT/ave_f0_1(i))**2+(sem_f1_0*RT/ave_f1_0(i))**2)
       if(abs(ave_f0_1(i)-ave_f1_0(i)).lt.1E-3)then
         delta_A(i)=C
         exit
       else
         delta_A(i)=RT*log(ave_f0_1(i)/ave_f1_0(i))+C
         C=delta_A(i)
       endif
     enddo
     sum_delta_A=sum_delta_A+delta_A(i)
     sum_sem_delta_A=sum_sem_delta_A+sem_delta_A(i)**2
     write(21,'(I16,2F16.4)')i,delta_A(i),sem_delta_A(i)
     write(*,'(A20,I12,A20,F12.1)')'window:',i,'Frame(BAR):',real(counter)
   enddo
   write(21,'(A16,2F16.4)')"Sum BAR:",sum_delta_A,sqrt(sum_sem_delta_A)
   close(21)
end subroutine barcal

