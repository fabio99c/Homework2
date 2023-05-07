module isingf

  implicit none
  public :: init,metropolis,DeltaEPBC,DeltaEOBC
  integer, parameter, public :: dp=selected_real_kind(13)
  integer, public, dimension(8) :: seed
  integer, public :: L,N,M, Nequi, Nmcs, ic
  real (kind = dp), public, dimension(-8:8) :: w
  integer, public, dimension(:,:), allocatable :: spin
  real (kind = dp), public, dimension(:), allocatable :: datae, datac
  real(kind=dp), public :: T,E,k,J
  character(len = 3),public :: bconditions
  character(len = 6),public :: mcsselection
  integer, public :: accept

contains

  subroutine init()
	integer :: dE
	real :: p
	ic=1
	seed=[1,2,3,4,5,6,7,8]
    call random_seed(put=seed)
	do dE = -8,8,4
		p=-J*real(dE)/real(k*T)
		w(dE) = exp(p)
	end do
  end subroutine init

  subroutine alloc()
	allocate(spin(L,L))
	allocate(datae((Nequi+Nmcs)*N))
    allocate(datac((Nequi+Nmcs)*N))
  end subroutine alloc

  subroutine metropolis()
    !  one Monte Carlo step per spin
    integer :: ispin,x,y,dE
    real(kind=dp) :: rnd,p
    if (mcsselection=="random") then
        do ispin = 1,N
            call random_number(rnd)
            x = int(L*rnd) + 1
            call random_number(rnd)
            y = int(L*rnd) + 1
            call random_number(rnd)
            if (bconditions=="PBC") then
                dE = DeltaEPBC(x,y)
                p=w(dE)
            else if (bconditions=="OBC") then
                dE = DeltaEOBC(x,y)
                p=-J*dE/(T*k)
                rnd = log(rnd)
            end if
            if (rnd <= p) then
                spin(x,y) = -spin(x,y)
                accept = accept + 1
                M = M + 2*spin(x,y)  ! factor 2 is to account for the variation:
                E = E + dE       ! (-(-)+(+))
            end if
			datae(ic)=E
			datac(ic)=M
			ic=ic+1
        end do
    else if (mcsselection =="ordered") then
        do x = 1,N
            do y = 1,N
                call random_number(rnd)
                if (bconditions=="PBC") then
                    dE = DeltaEPBC(x,y)
                    p=w(dE)
                else if (bconditions=="OBC") then
                    dE = DeltaEOBC(x,y)
                    p=-J*dE/(T*k)
                    rnd = log(rnd)
                end if
                if (rnd <= p) then
                    spin(x,y) = -spin(x,y)
                    accept = accept + 1
                    M = M + 2*spin(x,y)  ! factor 2 is to account for the variation:
                    E = E + dE           ! (-(-)+(+))
                end if
				datae(ic)=E
				datac(ic)=M
				ic=ic+1
            end do
        end do     
    end if
	
  end subroutine metropolis

  function DeltaEPBC(x,y) result (DeltaE_result)
    !  periodic boundary conditions
    integer, intent (in) :: x,y
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down
    if (x == 1) then
       left = spin(L,y)
       right = spin(2,y)
    else if (x == L) then
       left = spin(L-1,y)
       right = spin(1,y)
    else
       left = spin(x-1,y)
       right = spin(x+1,y)
    end if
    if (y == 1) then
       up = spin(x,2)
       down = spin(x,L)
    else if (y == L) then
       up = spin(x,1)
       down = spin(x,L-1)
    else
       up = spin(x,y+1)
       down = spin(x,y-1)
    end if
    DeltaE_result = 2*spin(x,y)*(left + right + up + down)
! also here the factor 2 is to account for the variation
  end function DeltaEPBC
  
    
  function DeltaEOBC(x,y) result (DeltaE_result)
    !  open boundary conditions
    integer, intent (in) :: x,y
    integer :: DeltaE_result
    integer :: left
    integer :: right
    integer :: up
    integer :: down
    if (x == 1) then
       left = 0
       right = spin(2,y)
    else if (x == L) then
       left = spin(L-1,y)
       right = 0
    else
       left = spin(x-1,y)
       right = spin(x+1,y)
    end if
    if (y == 1) then
       up = spin(x,2)
       down = 0
    else if (y == L) then
       up = 0
       down = spin(x,L-1)
    else
       up = spin(x,y+1)
       down = spin(x,y-1)
    end if
    DeltaE_result = 2*spin(x,y)*(left + right + up + down)
! also here the factor 2 is to account for the variation
  end function DeltaEOBC
  
end module isingf

