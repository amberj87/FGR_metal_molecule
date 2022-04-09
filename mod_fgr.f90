Module mod_fgr
!! JPCA_117_6196
implicit none
real*8, parameter :: clight=2.99792458D10,av=6.0221367D23,hbar=1.05457266D-34
real*8, parameter :: mass_h=1.007825d0,kb=1.3806503d-23
real*8, parameter :: pascal_to_atm=9.86923267d-6,kcal_to_J=4184.d0
real*8, parameter :: amu2kg=1.66053892d-27
real*8, parameter :: au2kg=9.10938291d-31,au2J=4.35974417d-18,au2m=5.2917721092d-11
complex*16,parameter :: iota = (0.d0,1.d0)
real*8 pi,wave_to_J

!! Common Variables
real*8 V_coup,V_exothermicity,gama_coup,V_ex
real*8 omg1,V_reorg1,g_coup1
real*8 omg2,V_reorg2,g_coup2
real*8 gamma
real*8 temperature,mass,x_cr,beta
real*8 VER_rate

!! integration
real*8 wmax,tmax,Emax,dE
integer nclass,n_t,nE
real*8,allocatable:: omg(:),ck(:)

!! diag
integer nold
real*8,allocatable:: work(:)
integer,allocatable:: iwork(:),isuppz(:)

contains
!---------------------------------------------------------- 
!---------------------------------------------------------- 

subroutine setup
  implicit none
  character st_ch
  real*8 c_0,c_e

  pi=dacos(-1.d0)
  wave_to_J=2*pi*clight*hbar
  open(10,file="fgr.inp")
  read(10,*) mass
  read(10,*) gama_coup
  read(10,*) V_exothermicity
  read(10,*) omg1
  read(10,*) V_reorg1
  read(10,*) omg2
  read(10,*) gamma
  read(10,*) VER_rate
  read(10,*) temperature
  read(10,*) nclass
  read(10,*) wmax
  read(10,*) n_t
  read(10,*) tmax
  read(10,*) nE
  read(10,*) Emax
  read(10,*) st_ch
  close(10)
  !----------------------------------------------------------

  if(st_ch.ne.'x') then
    write(6,*) "problem in reading input file"
    stop
  endif

  !---------------------------------------------------------- 

  allocate(omg(nclass),ck(nclass))

  mass=mass*au2kg
  gama_coup=gama_coup*au2J
  omg1=omg1*(2*pi*clight)
  omg2=omg2*(2*pi*clight)
  gamma=gamma*(2*pi*clight)
  wmax=wmax*(2*pi*clight)
  V_exothermicity=V_exothermicity*wave_to_J
  V_coup=V_coup*wave_to_J
  V_reorg1=V_reorg1*wave_to_J
  Emax=Emax*au2J
  dE=2*Emax/real(nE)
  beta=1.d0/(kb*temperature)

  g_coup1=dsqrt(V_reorg1*mass*omg1**2/2.d0)
  call calculate_lambda

  !c_0=2*V_barrier/(mass*omg**2)
  !c_e=V_exothermicity/(2*mass*omg**2)

  !s01=-0.5d0*(dsqrt(c_0)+dsqrt(c_0+4*c_e))
!  s01=-dsqrt(V_reorg/(2*mass*omg**2))
  !s02=-s01
  !x_cr=V_exothermicity/(2*mass*omg**2*s01)

  !V_reorg=2*mass*omg**2*s01**2

!  s01=dsqrt(V_reorg/(2*mass*omg**2))

end subroutine setup
!---------------------------------------------------------- 

subroutine main
  implicit none
  integer i,j,k
  real*8 kM_f,kM_b
  real*8 kF_f,kF_b
  real*8 k_FGR,E


!write(6,*) (V_exothermicity-V_reorg)**2/(4*V_reorg)/wave_to_J
!do i=1,200
!  V_reorg=10.d0*i*wave_to_J
!write(6,*) i,(V_exothermicity-V_reorg)**2/(4*V_reorg)/wave_to_J
!enddo
!stop

!call draw_pes

  do k=1,4
    temperature=650+(k-1)*300/3.d0
    beta=1.d0/(kb*temperature)


    kM_f=0.d0
    kM_b=0.d0
    kF_f=0.d0
    kF_b=0.d0
    call convert_spectra(omg,ck)
    call compute_k_FGR(kF_f,kF_b)
    do i=1,nE
      E=-Emax + i*dE
      V_ex=V_exothermicity+E
      kM_f=kM_f+gama_coup*(1-fermi(E))*F_reorg(-V_ex)/hbar

      V_ex=-V_exothermicity-E
      kM_b=kM_b+gama_coup*(fermi(E))*F_reorg(-V_ex)/hbar
    enddo
    kM_f=kM_f*dE
    kM_b=kM_b*dE

    write(6,*) "Marcus, ps-1",kM_f/1.d12,kM_b/1.d12,(kM_f+kM_b)/1.d12
    write(6,*) 1/(1+kM_b/kM_f),1/(1+dexp(V_exothermicity/(kb*temperature)))
    write(6,*) "FGR, ps-1",kF_f/1.d12,kF_b/1.d12,(kF_f+kF_b)/1.d12
    write(6,*) 1/(1+kF_b/kF_f),1/(1+dexp(V_exothermicity/(kb*temperature)))
    write(6,*) 

    write(30,*) "Marcus, ps-1",kM_f/1.d12,kM_b/1.d12,(kM_f+kM_b)/1.d12
    write(30,*) 1/(1+kM_b/kM_f),1/(1+dexp(V_exothermicity/(kb*temperature)))
    write(30,*) "FGR, ps-1",kF_f/1.d12,kF_b/1.d12,(kF_f+kF_b)/1.d12
    write(30,*) 1/(1+kF_b/kF_f),1/(1+dexp(V_exothermicity/(kb*temperature)))
    write(30,*) 
    write(31,*)  temperature,(kM_f+kM_b)/1.d12,(kF_f+kF_b)/1.d12,1/(1+kF_b/kF_f),1/(1+dexp(V_exothermicity/(kb*temperature)))
  enddo


end subroutine main
!---------------------------------------------------------- 

subroutine compute_k_FGR(k_FGR_f,k_FGR_b)
  implicit none
  real*8,intent(out) :: k_FGR_f,k_FGR_B
  integer i,j
  real*8 t,tmin,dt
  real*8 wt(nclass)
  complex*16 su1,su2f,su2b,fac(nclass),su3,gg1,gg2
  real*8 t1,t2
  real*8 kk
  real*8 V_ex,E,g1(nE),en1(nE),g2(nE),en2(nE)

  call cpu_time(t1)

  !call convert_spectra(omg,ck)

  do j=1,nE
    E=-Emax + j*dE
    V_ex=V_exothermicity+E
    g1(j)=gama_coup*(1-fermi(E))
    en1(j)=V_ex
    g2(j)=gama_coup*fermi(E)
    en2(j)=-V_ex
  enddo


  tmin=0.d0
  dt=(tmax-tmin)/dfloat(n_t-1)
  
  su2f=0.d0
  su2b=0.d0
  !call omp_set_num_threads(8)
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,t,wt,fac,su1)
  !$OMP DO SCHEDULE(STATIC) REDUCTION(+:su2f) REDUCTION(+:su2b)
  do i=1,n_t
    t=tmin+(i-1)*dt
    wt=omg*t
    fac=(1-dcos(wt))/dtanh(0.5d0*beta*hbar*omg)-iota*dsin(wt)
    su1=2/mass*sum(ck**2/omg**3*fac)/hbar
    gg1 = sum(g1*cdexp(iota*en1*t/hbar))*dE/hbar
    gg2 = sum(g2*cdexp(iota*en2*t/hbar))*dE/hbar
write(30,*) t*1.d15,real(gg1),real(cdexp(-su1))
    su2f=su2f+cdexp(-su1)*gg1
    su2b=su2b+cdexp(-su1)*gg2
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  su2f=su2f*dt
  su2b=su2b*dt

  !k_FGR=2*V_coup**2*dble(su2)/hbar**2
  k_FGR_f=dble(su2f)/(2*pi*hbar)
  k_FGR_b=dble(su2b)/(2*pi*hbar)

  call cpu_time(t2)
  !write(6,*) t2-t1

end subroutine compute_k_FGR
!-----------------------------------------------------------------  

function density(w)
  implicit none
  real*8 density,w

  density=0.5*V_reorg2*omg2*omg2 * gamma*w/((w*w-omg2*omg2)**2+(gamma*w)**2)
!  density=mass*gamma*w

end function density
!-----------------------------------------------------------------  

function spec_brown(w)
  implicit none
  real*8 spec_brown,w

  !spec_brown=0.5*V_reorg1*omg1*omg1 * gamma*w/((w*w-omg1*omg1)**2+(gamma*w)**2)
  spec_brown=0.5*V_reorg1*omg1*omg1 * 10*omg1*w/((w*w-omg1*omg1)**2+(10*omg1*w)**2)
!  density=mass*omg1*w

end function spec_brown
!-----------------------------------------------------------------  


subroutine convert_spectra(omg_n,ck)
  implicit none
  real*8,intent(out):: omg_n(nclass),ck(nclass)
  integer i,j
  real*8 mat(nclass,nclass),en(nclass),vect(nclass,nclass)
  real*8 ck_brown(nclass),omg(nclass)
  real*8 delw

  delw=wmax/real(nclass)

  do i=1,nclass
    omg(i)=i*wmax/real(nclass)
    !ck(i)=dsqrt(2*spec_ohm(omg(i))*mass(i)*omg(i)*delw/pi)
    ck(i)=dsqrt(2*density(omg(i))*mass*omg(i)*delw/pi)
    !ck(i)=dsqrt(2*spec_brown(omg(i))*mass*omg(i)*delw/pi)
    !omg(i)=-omg_c*log((i-0.5)/real(nclass))
    !ck(i)=mass(i)*omg(i)*dsqrt(2*gamma_B*omg_c/(nclass*pi))
    write(20,*)omg(i)/(2*pi*clight),density(omg(i))
  enddo
  ck(1)=0.d0

  write(6,*)"reorg1", 2*sum(ck**2/(mass*omg**2)),V_reorg2

  mat=0.d0
  mat(1,1)=0.5*mass*omg1**2+sum(ck**2/(2*mass*omg**2))

  do i=2,nclass
    mat(i,i)=0.5*mass*omg(i)**2
  enddo
  do i=2,nclass
    mat(1,i)=ck(i)/2.d0
    mat(i,1)=ck(i)/2.d0
  enddo

  call diag(mat,nclass,en,vect,nclass)

  omg_n=dsqrt(2*en/mass)
  write(21,*) en

  ck=vect(1,:)*g_coup1
  write(6,*)"reorg2",2*sum(ck**2/(mass*omg_n**2))/wavE_to_J,V_reorg1/wave_to_J

  do i=1,nclass-1
    !ck_brown(i)=dsqrt(2*spec_brown(omg_n(i))*mass(i)*omg_n(i)*(omg_n(i+1)-omg_n(i))/pi)
    write(22,*) omg_n(i)/(2*pi*clight),pi/2.d0*ck(i)**2/(mass*omg_n(i)*(omg_n(i+1)-omg_n(i)))/wave_to_J,spec_brown(omg_n(i))/wave_to_J
  enddo

end subroutine convert_spectra
!-----------------------------------------------------------------  

subroutine calculate_lambda
  implicit none
  real*8 w

  w=omg1

  !V_reorg2=VER_rate*(1-dexp(-beta*hbar*omg1)) * mass*omg1
  V_reorg2=VER_rate*mass*omg1
  V_reorg2=V_reorg2 * 2 * ((w*w-omg2*omg2)**2+(gamma*w)**2)/(omg2**2*gamma*w)
  !V_reorg2=423.77084238981416

  g_coup2=dsqrt(V_reorg2*mass*omg2**2/2.d0)

end subroutine calculate_lambda
!-----------------------------------------------------------------  

subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0 

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
    allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!---------------------------------------------------------- 

function fermi(E)
  implicit none
  real*8 fermi,E

  fermi = 1.d0/(1.d0+dexp(E/(kb*temperature)))

end function fermi
!-----------------------------------------------------------------  

function F_reorg(x)
  implicit none
  real*8 F_reorg,x

  F_reorg = 1/dsqrt(4*pi*V_reorg1*kb*temperature) * dexp(-(x-V_reorg1)**2/(4*V_reorg1*kb*temperature))

end function F_reorg
!-----------------------------------------------------------------  


End Module mod_fgr
