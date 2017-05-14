!ifort lapack.f90 MainMatlab.f90 -lmkl_lapack95_lp64 -mkl ###### For first compile - need to have lapack.f90 in directory
!ifort -o MainMatlab MainMatlab.f90 -lmkl_lapack95_lp64 -mkl ###### After first compile

!Debugging flags ifort -o noness MainNonEssential.f90 -g -traceback -check all -fp-stack-check
 
      module int_stack8
!      * * F E A P * * A Finite Element Analysis Program

!....  Copyright (c) 1984-2016: Regents of the University of California
!                               All rights reserved

!-----[--.----+----.----+----.-----------------------------------------]
!     Modification log                                Date (dd/mm/year)
!       Original version                                    28/12/2016
!-----[--.----+----.----+----.-----------------------------------------]
!      Purpose:  Stack class for integers

!      Data:
!         data(:) - stack data    (public)
!         size    - size of stack (public)

!      Methods:
!         push  - push items onto stack
!         empty - check if stack is empty
!         
!-----[--.----+----.----+----.-----------------------------------------]
        implicit none

        public
 
        ! Type def for the stack data
        type int_stack_obj8
           integer(kind=8), allocatable :: data(:)
           integer              :: size = 0
        end type int_stack_obj8
 
        ! Block size
        integer, parameter, private :: block_size = 1024 !It shouldnt be too big here
 
 
      contains
 
        subroutine push8(s, e)
          ! Add element to stack, allocate if new
          ! Inputs: s - stack, e - integer to put on stack
          integer(kind=8), intent(in)                :: e
          integer(kind=8), allocatable               :: wk(:)
          type(int_stack_obj8), intent(inout) :: s

          if (.not. allocated(s%data)) then
             ! Allocate space if not yet done
             allocate(s%data(block_size))
       
          elseif (s%size == size(s%data)) then
             ! Grow the allocated space if last entry is full
             allocate(wk(size(s%data)+block_size))
             wk(1:s%size) = s%data
             call move_alloc(wk,s%data)   ! reset pointer and free wk
       
          end if
       
          ! Store the data in the stack
          s%size = s%size + 1
          s%data(s%size) = e
        end subroutine push8
 
        logical function empty8(s)
          ! Check if the stack is empty
          type(int_stack_obj8), intent(inout) :: s
          empty8 = (s%size == 0 .or. .not. allocated(s%data))
        end function empty8
 
      end module int_stack8

program NonEssentialLattice
      use int_stack8
      implicit none

      !!Integer (kind=8) variables
      integer (kind=8) :: v0(4), i0, i10, i20, j0, k0 !The indices yielding the minimum value
      integer (kind=8) :: sizeAL, sizeL  !Sizes of Permutations and set L for adding multiples of lattice vectors to shifts
      integer (kind=8) :: sizeSL1, sizeSL2, sizeSLinv, sizeSLinvIndex, sizeSLint1, sizeSLint2 !Sizes for loading binaries for SL sets
      integer (kind=8) :: sizeSL10, sizeSL20, sizeSLinv0, sizeSLinvIndex0, sizeSLint10, sizeSLint20 !Sized for loading binaries for preliminary search    
      integer (kind=8) :: LengthOfMinParameters, LengthOfMinParameters0 !Length of the int_stack8 objects
      integer (kind=8) :: i, j, m !Counter variables
      integer (kind=8) :: Factorial !Counter variables
      
      !!Integer (kind=2) variables
      integer(kind=2) :: detvalue, Det


      !!Integer (kind=1) variables
      integer (kind=1) :: N, N1, N2, NS, NS1, NS2, Kfactor !Lattice numbers and number if shifts, N-lattice and NS - number of shifts
      integer (kind=1) :: gcd0, gcd
      integer (kind=1) :: SLkmax1, SLkmax2, LkMax !Maximal entries of parent and product SL set
      integer (kind=1) :: j1, j2 !determinant of parent and product lattice transformation
      integer (kind=1) :: Mu(9), Mu0(9) !The running Mu for full and preliminary search
      integer (kind=1) :: Mu10(3,3), Mu20(3,3) !The final minimal Mu's giving rise to the minimal deformation
      integer (kind=1), dimension(:),   allocatable :: muint1, muint2, SL1, SLint1, SL2, SLint2, SL10, SLint10, SL20, SLint20 !The sets for storage in RAM
      integer (kind=1), dimension(:),   allocatable :: L, AL            !The sets L and AL to iterate over
      integer (kind=1), dimension(:,:), allocatable :: Perm, Lelement   !Singe elements in Al and L

      !!Integer variables
      integer :: pos(2), pos0(2)    !The positions from SLinvIndex/SLinvIndex0

      !Int_stack objects
      type(int_stack_obj8) :: MinParameters, MinParameters0 !Each quadruple gives the index of SLproduct, SLparent, AL and L (in this order)

      !Real variables
      real :: m0,m0tot !The minimal value with and without the shifts 
      real :: m1,m1tot !The potential minimal value with and without the shifts 
      real :: d        !The determinant of the SLinv matrices
      real :: start, finish, finish0 !Timing variables
      real :: F(3,3), G(3,3), Fi(3,3), Hmu(3,3), HmuMin(3,3), sing(3), normsq, snormsq, norm,snorm
      real :: F1(3,3), G1(3,3)  !The supercell lattice vectors
      real :: Mu10i(3,3)        !Mu10 inverse
      real, dimension(:,:), allocatable :: P, Q, P1, Q1 !The original shift vectors/super call shift vectors

      !Real parameters
      real, parameter :: tol=10.**(-5.) !Tolerance to find other minima
      real, parameter :: r=2./2.      !norm for d_r distance - divided by 2 (it's because we are using the square singvals )
      real, parameter :: alpha=.4     !Weight for shifts
      real, parameter :: qpower=2     !Shift norm exponent

      !Logical variables
      logical :: FoundInPreliminary=.true. !True if the minimum was already found in the preliminary search

      !Chracters
      character*1  :: NString, LkString
      character*1  :: kString1, kString2
      character*1  :: iString1, iString2
      character*2  :: kiString1, kiString2
      character*7  :: LString
      character*10 :: AlphaString
      character*8  :: SLkiString1, SLkiString2, SLkiString10, SLkiString20
      character*11 :: SLkiIntString1, SLkiIntString2, SLkiIntString10, SLkiIntString20
      character*14 :: SLinvString, SLinvString0
      character*19 :: SLinvStringindex, SLinvStringindex0


      !Initialisation of counters
      start=0
      finish=0
      LengthOfMinParameters=0
      LengthOfMinParameters0=0
      call cpu_time(start)

      ! If Kfactor>1, then Kafcctor times bigger sublattice than necessary are considered
      Kfactor=1


     
! -------------------- Test with Bain Strain ------------------------------------------------------!

! -------------------- The input lattices  --------------------------------------------------------!

      !Parent
      F=real(reshape((/ 1,0,0,0,1,0,0,0,1 /), shape(F)))
      write(*,*) F
      NS1=3 !Number of Shifts for the parent lattice
      allocate(P(3,NS1))
      P=.5*real(reshape((/ 0, 1, 1, 1, 0, 1, 1, 1, 0 /), shape(P)))
      write(*,*) P  
      N1=NS1+1 !Thus the parent lattice is an N1 lattice

      !Product
      G=2.**(-1./3.)*F
      write(*,*) G
      NS2=1 !Number of Shifts for the product lattice
      allocate(Q(3,NS2))
      Q=.5*2.**(-1./3.)*real(reshape((/1,1,1/), shape(Q)))
      N2=NS2+1 !Thus the product lattice is an N2 lattice

      !Determine the correct lattice transformation
      gcd0=gcd(N1,N2)
      N=Kfactor*N1*N2/gcd0 ! Thus it is a transformation between N lattices
      j1=Kfactor*N2/gcd0
      j2=Kfactor*N1/gcd0

      write(*,*) j1
      write(*,*) j2
      print '("The transformations is between",I2," lattices.")', N


      d=real(1./real(j1)) !The parent transformations needs this factor for their inverse

      allocate(P1(3,j1*n1-1)) !They have the same size
      allocate(Q1(3,j2*n2-1)) !They have the same size

      allocate(MuInt1(2*(j1-1)))
      allocate(MuInt2(2*(j2-1)))

      NS=N-1 !Number of Shifts for both lattices
      allocate(Perm(NS,NS))
      allocate(Lelement(3,NS))

      !The size of the search sets

      SLkmax1=2 !The absolut maximal entries for the parent lattice 
      SLkmax2=1 !The absolut maximal entries for the product lattice 
      Lkmax = 1 !The absolut maximal entries for L - this should probably stay at 1 

      write(NString, 45) N
45    format (I1) 
      write(LkString, 48) Lkmax
48    format (I1)

!!!!!!!!!!!!!!! Open the SL files and SLint files and read to RAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'Opening the necessary sets'
      !Set Parent lattice SL strings
      write(kString1, 15) SLkmax1
15    format (I1) 
      write(iString1, 16) j1
16    format (I1) 
      SLkiString1 = trim('SL') // trim(kString1) // trim(iString1)// trim('.bin')
      write(*,*) SLkiString1

      ! Open Parent lattice SL 
      open(file=SLkiString1,unit=30,form='unformatted', action='read',status='old',access='stream')
      read(30) sizeSL1
      write(*,*) sizeSL1
      allocate(SL1(9*sizeSL1))
      read(30) SL1
      write(*,*) SL1(1:9)
      close(30)

      ! Open Parent lattice SLint
      SLkiIntString1 = trim('SLint') // trim(kString1) // trim(iString1)// trim('.bin')
      write(*,*) SLkiIntString1
      open(file=SLkiIntString1,unit=100,form='unformatted', action='read',status='old',access='stream')
      read(100) sizeSLint1

      if (j1==1) then !If the determinant=1 then there are no new Interior vectors
        allocate(SLint1(3*sizeSL1)) !Need to check this convention
        SLint1=0
      else
        allocate(SLint1(3*(j1-1)*sizeSL1))
        read(100) SLint1
        write(*,*) SLint1(1:3)
      end if
      close(100)

      ! Open Parent lattice SL 
      write(kString2, 17) SLkmax2 
17    format (I1) 
      write(iString2, 18) j2
18    format (I1) 
      SLkiString2 = trim('SL') // trim(kString2) // trim(iString2)// trim('.bin')
      write(*,*) SLkiString2

      open(file=SLkiString2,unit=31,form='unformatted', action='read',status='old',access='stream')
      read(31) sizeSL2
      write(*,*) sizeSL2
      allocate(SL2(9*sizeSL2))
      read(31) SL2
      write(*,*) SL2(1:9)
      close(31)

      ! Open Product lattice SLint
      SLkiIntString2 = trim('SLint') // trim(kString2) // trim(iString2)// trim('.bin')
      write(*,*) SLkiIntString2
      open(file=SLkiIntString2,unit=101,form='unformatted', action='read',status='old',access='stream')
      read(101) sizeSLint2


      if (j2==1) then !If the determinant=1 then there are no new interior vectors
        allocate(SLint2(3*sizeSL2)) !Need to check this convention
        SLint2=0
      else
        allocate(SLint2(3*(j2-1)*sizeSL2))
        read(101) SLint2
        write(*,*) SLint2(1:3)
      end if
      close(101)


      !!!Open the SLinv sets (incl index)
      SLinvString = trim('SLinv') // trim(kString2) // trim(iString2)// trim('o') // trim(kString1) &
                  // trim(iString1)// trim('.bin')
      write(*,*) SLinvString
      !Open the SLinv set
      open(file=SLinvString,unit=20,form='unformatted', action='read',status='old',access='stream')
      read(20) sizeSLinv !Number of matrices (not elements)
      write(*,*) sizeSLinv
      read(20) Mu
      write(*,*) Mu
      rewind(20)
      read(20) sizeSLinv 
      
      !Set the SLinvindex string
      SLinvStringindex = trim('SLinv') // trim(kString2) // trim(iString2)// trim('o') // trim(kString1) &
                    // trim(iString1)// trim('index.bin')
      write(*,*) SLinvStringindex
      !Open the SLinvindex set
      open(file=SLinvStringindex,unit=21,form='unformatted', action='read',status='old',access='stream')
      read(21) sizeSLinvIndex !Number of tuples (not elements)
      write(*,*) sizeSLinvIndex
      read(21) pos
      write(*,*) pos
      rewind(21)
      read(21) sizeSLinvIndex 


      !Open Alphaset and store in RAM
      AlphaString=trim('Alpha') // trim(NString) // trim('.bin')
      write(*,*) AlphaString
      open(file=AlphaString,unit=10,form='unformatted', action='read',status='old',access='stream')
      read(10) sizeAL !This is the number of elements which is (k+1)!*k^2 (they are (k+1)! kxk matrices)
      allocate(AL(sizeAL)) !Allocate array
      read(10) AL         !Write everything to allocated array
      sizeAL=sizeAL/((NS)**2) !This is the number of matrices in the set
      write(*,*) sizeAL
      write(*,*) AL(1:NS**2)
      close(10) 

      !Open LSet and store in RAM
      LString=trim('L') // trim(NString) // trim(LkString) // trim('.bin')
      write(*,*) LString
      open(file=LString,unit=11,form='unformatted', action='read',status='old',access='stream')
      read(11) sizeL !This is the number of elements which is (k+1)!*k^2 (they are (k+1)! kxk matrices)
      allocate(L(sizeL)) !Allocate array
      read(11) L         !Write everything to allocated array
      sizeL=sizeL/(3*NS) !This is the number of matrices in the set
      write(*,*) sizeL
      write(*,*) L(1:3*NS)
      close(11)

      !Post a message to use the other routine if it is a transformationen between essential N lattices; proceed otherwise
      if (j1 ==1 .and. j2==1) then
        write(*,*) 'This is a transformation between two essential N lattices - use the corresponding routine'
      else
      !Open the analogous sets for preliminary search
      write(*,*) 'Opening the necessary sets for a preliminary search'

      !Open the SLinv0 set
      SLinvString0 = trim('SLinv') // trim('1') // trim(iString2)// trim('o') // trim('1') &
                    // trim(iString1)// trim('.bin')
      write(*,*) SLinvString0
      open(file=SLinvString0,unit=70,form='unformatted', action='read',status='old',access='stream')
      read(70) sizeSLinv0 !Number of matrices (not elements)
      write(*,*) sizeSLinv0
      read(70) Mu0
      write(*,*) Mu0
      rewind(70)
      read(70) sizeSLinv0 

      !Open the SLinvindex0 set
      SLinvStringindex0 = trim('SLinv') // trim('1') // trim(iString2)// trim('o') // trim('1') &
                    // trim(iString1)// trim('index.bin')
      write(*,*) SLinvStringindex0
      open(file=SLinvStringindex0,unit=71,form='unformatted', action='read',status='old',access='stream')
      read(71) sizeSLinvIndex0 !Number of tuples (not elements)

      ! Open Parent lattice SL0
      SLkiString10 = trim('SL') // trim('1') // trim(iString1)// trim('.bin')
      write(*,*) SLkiString10
      open(file=SLkiString10,unit=130,form='unformatted', action='read',status='old',access='stream')
      read(130) sizeSL10
      allocate(SL10(9*sizeSL10))
      read(130) SL10
      close(130)

      ! Open Parent lattice SLint0
      SLkiIntString10 = trim('SLint') // trim('1') // trim(iString1)// trim('.bin')
      write(*,*) SLkiIntString10
      open(file=SLkiIntString10,unit=140,form='unformatted', action='read',status='old',access='stream')
      read(140) sizeSLint10

      if (j1==1) then !If the determinant=1 then there are no new interior vectors
        allocate(SLint10(3*sizeSL10)) !Need to check this convention
        SLint10=0
      else
        allocate(SLint10(3*(j1-1)*sizeSL10))
        read(140) SLint10
      end if
      close(140)


      ! Open Product lattice SL0
      SLkiString20 = trim('SL') // trim('1') // trim(iString2)// trim('.bin')
      write(*,*) SLkiString20
      open(file=SLkiString20,unit=131,form='unformatted', action='read',status='old',access='stream')
      read(131) sizeSL20
      allocate(SL20(9*sizeSL20))
      read(131) SL20
      close(131)

      !Open Product lattice SLint0
      SLkiIntString20 = trim('SLint') // trim('1') // trim(iString2)// trim('.bin')
      write(*,*) SLkiIntString20
      open(file=SLkiIntString20,unit=141,form='unformatted', action='read',status='old',access='stream')

      if (j2==1) then   !If the determinant=1 then there are no new interior vectors
        allocate(SLint20(3*sizeSL20)) !Need to check this convention
        SLint20=0
      else
        allocate(SLint20(3*(j2-1)*sizeSL20))
        read(141) SLint20
      end if
      write(*,*) sizeSL20
      close(141)


      !Compute the inverse of F
      call inverse(F,Fi)
      ! write(*,*) Fi

! --------------------- We first do a preliminary minimisation step to find a m_0  -----------------------------------------!

      !Set some arbitrary first value - the sets SLinv are ordered such that the row with the least absolute sum comes first
      read(70) Mu0
      read(71) pos0
      Hmu=d*MATMUL(MATMUL(G, reshape(real(Mu0), shape(G))),Fi)
      call SingvalsSq(Hmu, sing) 

      !Calculate (only) the skeletal deformation distance
      m0=norm(sing**r-[1,1,1]) 

      !Calculate the Parent supercell
      if (j1==1) then
        call SuperCell(F,P,SL1(1:9),0,j1,N1,F1,P1) !In this case there are no new  no SLint vector
      else
        call SuperCell(F,P,SL1(1:9),SLint10(1:3*(j1-1)),j1,N1,F1,P1) 
      end if

      !Calculate the Product SuperCell
      if (j2==1) then
        call SuperCell(G,Q,SL2(1:9),0,j2,N2,G1,Q1) !In this case there are no new  no SLint vector
      else
        call SuperCell(G,Q,SL2(1:9),SLint20(1:3*(j2-1)),j2,N2,G1,Q1)
      end if

      !Calculate the (full) deformation distance
      m0tot=m0+alpha*snorm(Q1-MATMUL(Hmu,P1),qpower)

      !Cover the unlikely case that this is already the optimal
      I10=pos0(2)
      I20=pos0(1)
      J0=int8(1)
      K0=int8(1) !this is the 0 matrix 
      call push8(MinParameters0,int8(pos0(2))) !Product
      call push8(MinParameters0,int8(pos0(1))) !Parent
      call push8(MinParameters0,int8(j0)) !This is the Alpha element identity
      call push8(MinParameters0,int8(k0)) !This is the 0 matrix in L 
      LengthOfMinParameters0=LengthOfMinParameters0+1

      !Do the preliminary minimisation with only SLinv0 and Alpha  
      ! --> here we are only looking for the absolute minimum (without tolerance)
      do i=2,sizeSLinv0 !Starts at 2 since we have already used the first element
        read(70) Mu0    !Read in mu0 from SLinv0
        read(71) pos0   !Read in the corresponding positions from SLinvindex0
        Hmu=d*MATMUL(MATMUL(G, reshape(real(Mu0), shape(G))),Fi)
        
        call SingvalsSq(Hmu, sing) 
        m1=norm(sing**r-[1,1,1]) !Set the current skeletal deformation distance

        if (m1<m0tot) then !Only enter the loop if the skeletal part is smaller than m0tot
             !Parent supercell
             call SuperCell(F,P,SL10(pos0(2)*9-8+1:pos0(2)*9),SLint10((pos0(2)-1)*3*(j1-1)+1:pos0(2)*3*(j1-1)),j1,N1,F1,P1)
             !Product SuperCell
             call SuperCell(G,Q,SL20(pos0(1)*9-8+1:pos0(1)*9),SLint20((pos0(1)-1)*3*(j2-1)+1:pos0(1)*3*(j2-1)),j2,N2,G1,Q1)

             do j=1,sizeAL !Try all possible lattice permutations
              m1tot=m1+alpha*snorm(Q1-MATMUL(MATMUL(Hmu,P1),real(reshape(AL((j-1)*NS**2+1:j*NS**2),shape(Perm)))),qpower)
              if (m1tot<m0tot) then !If the total deformation distance is smaller than the previos one - add to the list
                call push8(MinParameters0,int8(pos0(2))) !Product
                call push8(MinParameters0,int8(pos0(1))) !Parent
                call push8(MinParameters0,int8(j)) !This is the Alpha element
                call push8(MinParameters0,int8(1)) !This is the 0 matrix in L 
                LengthOfMinParameters0=LengthOfMinParameters0+1
                m0tot=MIN(m1tot,m0tot)
              end if
            end do !Alphaset
        end if
      end do !SLinv0
      !Display time taken for preliminary search
      ! write(*,*) 'Preliminary search completed: time'
      call cpu_time(finish0)
      print '("Preliminary search completed in ",e12.5," seconds.")',finish0-start

!------------- Do the full minimisation (this additionally involves the set L and potentially bigger ki's) -------------------------------!
      do i=1,sizeSLinv
        read(20) Mu
        read(21) pos
        Hmu=d*MATMUL(MATMUL(G, reshape(real(Mu), shape(G))),Fi)
        call SingvalsSq(Hmu, sing) 
        m1=norm(sing**r-[1,1,1]) !Just the skeletal part
        ! print '("m1 ",e12.5," .")',m1
        ! print '("tol ",e12.5," .")',tol
        ! print '("m0tot ",e12.5," .")',m0tot+tol
        if (m1<m0tot+tol) then !If just the skeletal part is already bigger than m0tot then do no enter ther loop
          !Parent supercell
          call SuperCell(F,P,SL1(pos(2)*9-8:pos(2)*9),SLint1((pos(2)-1)*3*(j1-1)+1:pos(2)*3*(j1-1)),j1,N1,F1,P1)
          !Product SuperCell
          call SuperCell(G,Q,SL2(pos(1)*9-8:pos(1)*9),SLint2((pos(1)-1)*3*(j2-1)+1:pos(1)*3*(j2-1)),j2,N2,G1,Q1)
          do j=1,sizeAL !Try all possible lattice permutations
            do m=1,sizeL !Add skeletal vectors to the shift vectors
              m1tot=m1+alpha*snorm(Q1+MATMUL(G1,real(reshape(L((m-1)*3*NS+1:m*3*NS),shape(Lelement)))) & !???? Is it G or G1??? - probably G1 since Hmu maps F1 to G1
                    -MATMUL(MATMUL(Hmu,P1),real(reshape(AL((j-1)*NS**2+1:j*NS**2),shape(Perm)))),qpower)
              if (m1tot<m0tot+tol) then 
                !If the total deformation distance is smaller than the previos one (+ tolerance)
                !-> add to the list
              call push8(MinParameters,int8(pos(2))) !Product
              call push8(MinParameters,int8(pos(1))) !Parent
              call push8(MinParameters,int8(j)) !This is the Alpha element
              call push8(MinParameters,int8(m)) !This is the L element
              LengthOfMinParameters=LengthOfMinParameters+1
              m0tot=MIN(m1tot,m0tot)
              end if
            end do !Lset
          end do !AlphaSet
        end if !(m1<m0tot)
      end do !Slinv
      end if !non-essential case


!------------------------------------------ Output results and stats -------------------------------------------------!
      call cpu_time(finish)
      print '("Full search completed in ",e12.5," seconds.")',finish-start

      !Display the results

      !Cover the case that it was already found in the preliminary search - this can only happen if tol<eps
      FoundInPreliminary=empty8(MinParameters) 
      if (FoundInPreliminary==.true.) then
        print '("Minimal deformation was found in Preliminary Search")'
        ! write(*,*) LengthOfMinParameters0
        ! write(*,*) MinParameters0%data(Max(1,(LengthOfMinParameters0-10)*4+1):LengthOfMinParameters0*4)
        print '("Minimal total deformation distance is ",e12.5)', m0tot
        v0=MinParameters0%data(LengthOfMinParameters0*4-3:LengthOfMinParameters0*4)
        i10=v0(1)
        i20=v0(2)
        j0=v0(3)
        k0=v0(4)
        write(*,*) i10
        write(*,*) i20
        write(*,*) j0
        write(*,*) k0
        Mu10=Transpose(reshape(SL10(9*(i10-1)+1:9*i10),shape(F))) !These transposes shouldnt be here but oh well
        call inverse(real(Mu10),Mu10i)
        Mu20=reshape(SL20(9*(i20-1)+1:9*i20),shape(F))
        Perm=reshape(AL((j0-1)*NS**2+1:j0*NS**2),shape(Perm))
        Lelement=reshape(L((k0-1)*3*NS+1:k0*3*NS),shape(Lelement))
        Hmu=Transpose(MATMUL(MATMUL(MATMUL(G, real(Mu20)),Mu10i),Fi)) !These transposes shouldnt be here but oh well
      else
        ! write(*,*) LengthOfMinParameters
        ! write(*,*) MinParameters%data(Max(1,(LengthOfMinParameters-10)*4+1):LengthOfMinParameters*4)
        print '("Minimal total deformation distance is ",e12.5)', m0tot
        v0=MinParameters%data(LengthOfMinParameters*4-3:LengthOfMinParameters*4)
        i10=v0(1)
        i20=v0(2)
        j0=v0(3)
        k0=v0(4)
        write(*,*) i10
        write(*,*) i20
        write(*,*) j0
        write(*,*) k0
        Mu10=Transpose(reshape(SL1(9*(i10-1)+1:9*i10),shape(F))) !These transposes shouldnt be here but oh well
        call inverse(real(Mu10),Mu10i)
        Mu20=reshape(SL2(9*(i20-1)+1:9*i20),shape(F))
        Perm=reshape(AL((j0-1)*NS**2+1:j0*NS**2),shape(Perm))
        Lelement=reshape(L((k0-1)*3*NS+1:k0*3*NS),shape(Lelement))
        Hmu=Transpose(MATMUL(MATMUL(MATMUL(G, real(Mu20)),Mu10i),Fi)) !These transposes shouldnt be here but oh well
      end if

      !Show the final minimal matrices
      write(*,*) Mu10
      write(*,*) Mu20
      write(*,*) int1(Perm)
      write(*,*) int1(Lelement)
      write(*,*) Hmu
end program NonEssentialLattice

! -------------------------------------------------------------------------------------------------------------------------------------------------- !
!                                                       SUBROUTINES AND FUNCTIONS                                                                    !
! -------------------------------------------------------------------------------------------------------------------------------------------------- !

subroutine SuperCell(F,P,mu,muint,j,n,G,Q)
      ! Returns the new lattice vectors for the non-essential lattice 
      !   F is the skeletal matrix
      !   P are the orginal shifts
      !   mu is the integer transformation 
      !   j=det(mu)
      !   muint are the skeletal shifts
      !   [~,n-1]=size(P) and so it is an n lattice
      IMPLICIT NONE
      real,             intent(in)  :: F(3,3), P(3,n-1)
      integer(kind=1),  intent(in)  :: j,n,mu(9), muint(3*(j-1))
      real,             intent(out) :: G(3,3), Q(3,j*n-1)
      real                          :: Psk(3,j-1), Psh(3,n-1,j) !skeletal and Interior shift vectors
      integer(kind=1)               :: i

      !Cover the case of not having any shift vectors
      if (n==1 .and. (j==1)) then 
        G=MatMul(F,real(reshape(mu,shape(F))))
      else if (n==1) then
        G=MatMul(F,real(reshape(mu,shape(F)))) !This computes the new skelateal vectors and outputs them as a matrix
        Q=MatMul(F,real(reshape(muint, shape(Q)))) ! This computes the skeletal shift vectors
      
      !Cover the case of not having any additional shift vectors, i.e. det(mu)=1
      else if (j==1) then
        Q=P
        G=MatMul(F,real(reshape(mu,shape(F))))

      !Cover the generic case
      else
        G=MatMul(F,real(reshape(mu,shape(F))))  ! This computes the new skelateal vectors and outputs them as a matrix
        Psk=MatMul(F,real(reshape(muint,shape(Psk))))       ! This computes the skeletal shift vectors
        Psh(:,:,1)=P !The first block are the orignal shifts
        do i=2,j !This adds the skeletal shift vectors to the orignal shifts
          Psh(1,:,i) = P(1,:) + Psk(1,i-1)
          Psh(2,:,i) = P(2,:) + Psk(2,i-1)
          Psh(3,:,i) = P(3,:) + Psk(3,i-1)
        end do
        Q(1:3,1:3)=Psk 
        Q(1:3,4:j*n-1)=reshape(Psh,(/3, (n-1)*j/) ) ! this puts them all into one vector
      end if
end subroutine supercell



subroutine inverse(A,B)
      IMPLICIT NONE
      real,intent(in) :: A(3,3)
      real, intent(out)::B(3,3)

      B=Transpose(reshape( &
      [-(A(2,3)*A(3,2)) + A(2,2)*A(3,3),  &
      A(1,3)*A(3,2) - A(1,2)*A(3,3),      &
      -(A(1,3)*A(2,2)) + A(1,2)*A(2,3),   &
      A(2,3)*A(3,1) - A(2,1)*A(3,3),      &
      -(A(1,3)*A(3,1)) + A(1,1)*A(3,3),   &
      A(1,3)*A(2,1) - A(1,1)*A(2,3),      &
      -(A(2,2)*A(3,1)) + A(2,1)*A(3,2),   &
      A(1,2)*A(3,1) - A(1,1)*A(3,2),      &
      -(A(1,2)*A(2,1)) + A(1,1)*A(2,2)],shape(B))) &
      / (-(A(1,3)*A(2,2)*A(3,1)) + A(1,2)*A(2,3)*A(3,1) + A(1,3)*A(2,1)*A(3,2) &
      -  A(1,1)*A(2,3)*A(3,2) - A(1,2)*A(2,1)*A(3,3) + A(1,1)*A(2,2)*A(3,3))
end subroutine inverse

real function norm(a) !Calculates the 2-norm of a
      IMPLICIT NONE
      real :: a(3)

      norm=SQRT(DOT_PRODUCT(a,a))
end function

real function normsq(a) !Calculates the square of the 2-norm of a
      IMPLICIT NONE
      real :: a(3)

      normsq=DOT_PRODUCT(a,a)
end function

real function snorm(a,s) !Calculates the the s-norm of a
      IMPLICIT NONE
      real :: a(3), s

      if (s==2) then
        snorm=SQRT(DOT_PRODUCT(a,a))
      end if
      snorm=DOT_PRODUCT(a,a)**(s/2)
end function

real function snormsq(a,s) !Calculates the square of the s-norm of a
      IMPLICIT NONE
      real :: a(3), s

      if (s==2) then
        snormsq=DOT_PRODUCT(a,a)
      end if
      snormsq=(DOT_PRODUCT(a,a))**s
end function



subroutine SingvalsSq(B,s) !returns the singular values squared
      IMPLICIT NONE
      real,intent(in) :: B(3,3)
      real, intent(out)::s(3)
      real :: A(3,3)

      A=MATMUL(Transpose(B),B)
      call cubroot(-a(1,1) - a(2,2) - a(3,3),&
      (-a(1,2)**2 - a(1,3)**2 + a(1,1)*a(2,2) - a(2,3)**2 + a(1,1)*a(3,3) + &
      a(2,2)*a(3,3)),&
      a(1,3)**2*a(2,2) - 2*a(1,2)*a(1,3)*a(2,3) + a(1,1)*a(2,3)**2 + &
      a(1,2)**2*a(3,3) - a(1,1)*a(2,2)*a(3,3),s)
end subroutine SingvalsSq

subroutine cubroot(b,c,d,v) !Solves a cubic equation x^3+bx^2+cx+d=0 with three real solutions and returns them as the vector v
      IMPLICIT NONE
      real, intent(in) :: b,c,d
      real, intent(out) :: v(3)
      real :: d0, d1, Cr, root
      complex, parameter :: xi1=(-.5,0.8660254037844386), xi2=(-.5,-0.8660254037844386), xi3=(1.,0.)
      complex :: Cc, d0c,d1c, bc

      bc=cmplx(b)
      d0c=cmplx(b**2-3*c)
      d1c=cmplx(2*b**3-9*b*c+27*d)
      CC=(1./2.*(d1c+Sqrt(d1c**2-4.*d0c**3)))**(1./3.)
      v(1)=real(-cmplx(1./3.)*(bc+xi1*CC+d0c/(xi1*cc))) 
      v(2)=real(-cmplx(1./3.)*(bc+xi2*CC+d0c/(xi2*cc))) 
      v(3)=real(-cmplx(1./3.)*(bc+xi3*CC+d0c/(xi3*cc))) 
end subroutine cubroot

RECURSIVE FUNCTION Factorial(n)  RESULT(Fact)
      IMPLICIT NONE
      INTEGER(kind=8) :: Fact
      INTEGER(kind=8), INTENT(IN) :: n

      IF (n == 0) THEN
      Fact = 1
      ELSE
      Fact = n * Factorial(n-1)
      END IF
END FUNCTION Factorial

Integer(kind=1) function  GCD(ain,bin) !Computes the gcd of two positive integers using Euclid's method
      IMPLICIT  NONE
      INTEGER(kind=1)   :: a, b, c, ain, bin
      a=ain
      b=bin
      IF (a < b) THEN       ! since a >= b must be true, they
        c = a              ! are swapped if a < b
        a = b
        b = c
      END IF

      DO                    ! now we have a <= b
        c = MOD(a, b)      !    compute c, the reminder
        IF (c == 0) EXIT   !    if c is zero, we are done.  GCD = b
        a = b              !    otherwise, b becomes a
        b = c              !    and c becomes b
      END DO  
      gcd=b              !    go back
END function  GCD