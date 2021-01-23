! 
!     Copyright (C) 2021  Chee Kwan Gan (ihpcganck@gmail.com)
!     Copyright (C) 2020  Chee Kwan Gan (ihpcganck@gmail.com)
! 
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
! 
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
! 
! 

module commod
  use sgsym
  use derived_constants
  implicit none

  interface fargv
    module procedure getarg_str,getarg_int,getarg_double
  end interface

contains

  subroutine getarg_str(n,str)
    integer :: n
    character(len=*) :: str
    call get_command_argument(n,str)
  end subroutine getarg_str

  subroutine getarg_int(n,k)
    integer :: n, k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(I10)') k
  end subroutine getarg_int

  subroutine getarg_double(n,k)
    integer :: n
    real(double) :: k
    character(len=100) :: tmpstr
    call get_command_argument(n,tmpstr)
    read(unit=tmpstr,FMT='(D22.15)') k
  end subroutine getarg_double


  function str2N(str) result(n)
    character(len=*) :: str
    integer :: n
    read(str,'(I10)') n
  end function str2N

  function N2str(i) result (str)
    integer :: i,j,tmpv,intarr(20),nd,digit
    character(len=20) :: str

    if(i < 0) then
      tmpv = -i
    else
      tmpv = I
    endif

    nd = 0
    do
      digit = modulo(tmpv,10)
      nd = nd + 1
      intarr(nd) = digit
      tmpv = tmpv/10
      if(tmpv == 0) exit
    enddo
    str = ' '
    do j = 1, nd
      str(j:j) = digitarr(intarr(nd-j+1))
    enddo
    if(i < 0) then
      nd = nd + 1
      do j = nd,2,-1
        str(j:j) = str(j-1:j-1)
      enddo
      str(1:1) = '-'
    endif
  end function N2str

  subroutine fargn(n)
    integer :: n
    integer :: m

    m = command_argument_count()
    if(n /= m) then
      write(*,'(A,I4,A)') 'we need ',n,' arguments'
      write(*,'(A)') 'err: check number of arguments'
      stop 1
    endif
  end subroutine fargn

  function num_of_strings(tmp)
    character(len=*) :: tmp
    integer :: num_of_strings
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100

    character(len=100) :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)

      if(io /= 0) then
        exit
      endif
    enddo
    num_of_strings = i-1
    if(num_of_strings == max_num_cols) then

    endif
  end function num_of_strings

  function num_of_integers(tmp)
    character(len=*) :: tmp
    integer :: num_of_integers
    character(len=DSL) :: line
    integer,parameter :: max_num_cols = 100
    integer :: test_array(max_num_cols)
    integer :: i,io

    line = trim(tmp)
    write(*,'(A)') 'line is '//trim(line)
    do i = 1, max_num_cols

      read(line,*,iostat=io) test_array(1:i)
      if(io /= 0) then
        exit
      endif
    enddo
    num_of_integers = i-1
    if(num_of_integers == max_num_cols) then
      write(*,*) 'count_int_number= ',num_of_integers
      write(*,*) 'err: we may need to increase max_num_cols.'
      stop 1
    endif
  end function num_of_integers

  subroutine getlenang(a,L)
    real(double) :: a(3,3),L(6)
    real(double) :: sintheta,costheta

    L(1) = vecmag3(a(1:3,1))
    L(2) = vecmag3(a(1:3,2))
    L(3) = vecmag3(a(1:3,3))

    costheta = dotprod3(a(1:3,2),a(1:3,3))/(L(2)*L(3))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(4) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,3),a(1:3,1))/(L(3)*L(1))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(5) = atan2(sintheta,costheta)*rad2deg

    costheta = dotprod3(a(1:3,1),a(1:3,2))/(L(1)*L(2))
    sintheta = sqrt((1d0-costheta)*(1+costheta))
    L(6) = atan2(sintheta,costheta)*rad2deg
  end subroutine getlenang

  subroutine real2recip(reallatt,reciplatt)
    real(double) :: reallatt(3,3),reciplatt(3,3)
    real(double) :: vol,tmp2(3,3)

    vol = det3(reallatt(1,1))
    if(vol < 0) then
      write(*,*) 'vol = ',vol
      write(*,'(A)') 'err: vol is negative. Check reallatt.'
      stop 1
    endif

    tmp2 = inv3x3(reallatt(1,1))

    reciplatt(1:3,1) = tmp2(1,1:3)
    reciplatt(1:3,2) = tmp2(2,1:3)
    reciplatt(1:3,3) = tmp2(3,1:3)
  end subroutine real2recip

  function dotprod3(A,B)
    real(double) :: A(1:3),B(1:3),dotprod3
    dotprod3 = dotprodn(3,a(1),b(1))
  end function dotprod3

  function dotprodn(n,A,B)
    integer :: n,i
    real(double) :: A(n),B(n),dotprodn
    dotprodn = zero
    do i = 1, n
      dotprodn = dotprodn + A(i)*B(i)
    enddo
  end function dotprodn

  function tripprod(A)
    real(double) :: A(3,3),tripprod
    tripprod = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - &
               A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + &
               A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end function tripprod

  function det3(a)
    real(double) :: det3,a(3,3)
    det3 = tripprod(a)
  end function det3

  subroutine frac2abs(reallatt,FracC,AbsC)
    real(double) :: FracC(1:3),AbsC(1:3),reallatt(3,3)
    integer :: i,j
    do i = 1, 3
      AbsC(i) = 0
      do j = 1, 3
        AbsC(i) = AbsC(i) + reallatt(i,j)*FracC(j)
      enddo
    enddo
  end subroutine frac2abs

  subroutine abs2frac(reciplatt,AbsC,FracC)
    real(double) :: AbsC(3),FracC(3),reciplatt(3,3)
    integer :: i
    do i = 1, 3
      FracC(i) = dotprod3(AbsC(1:3),reciplatt(1:3,i))
    enddo
  end subroutine abs2frac

  function vecmag3(x)
    real(double) :: vecmag3,x(3)
    vecmag3 = sqrt(dotprod3(x,x))
  end function vecmag3

  function inv3x3(a) result(inv)
    real(double) :: a(3,3),inv(3,3)
    real(double) :: det,absdet

    inv(1,1) = a(2,2)*a(3,3)-a(2,3)*a(3,2)
    inv(2,1) = -(a(2,1)*a(3,3) - a(2,3)*a(3,1))
    inv(3,1) = a(2,1)*a(3,2) - a(2,2)*a(3,1)
    inv(1,2) = -(a(1,2)*a(3,3) - a(1,3)*a(3,2))
    inv(2,2) = a(1,1)*a(3,3) - a(1,3)*a(3,1)
    inv(3,2) = -(a(1,1)*a(3,2) - a(1,2)*a(3,1))
    inv(1,3) = a(1,2)*a(2,3) - a(1,3)*a(2,2)
    inv(2,3) = -(a(1,1)*a(2,3) - a(1,3)*a(2,1))
    inv(3,3) = a(1,1)*a(2,2) - a(1,2)*a(2,1)

    det = a(1,1)*inv(1,1)+a(1,2)*inv(2,1)+a(1,3)*inv(3,1)
    absdet = abs(det)

    if(absdet < 1.0d-13) then
      write(*,*) 'absdet = ',absdet
      write(*,'(A)') 'err: det is too small in inv3x3.'
      stop 1
    endif
    inv = inv/det
  end function inv3x3

  subroutine onespace(ns,comment1,comment)
    integer :: ns
    character(len=ns) :: comment1,comment
    integer :: i,j

    comment = 's'

    i = 1
    j = 1

    loop2: do
      loop1: do
        if(comment1(i:i) /= ' ') exit loop1
        i = i + 1
        if(i == ns) exit loop2
      enddo loop1

      loop3: do
        if(comment1(i:i) == ' ') exit loop3
        comment(j:j) = comment1(i:i)
        j = j + 1
        i = i + 1
      enddo loop3
      comment(j:j) = ' '
      j = j + 1
    enddo loop2
  end subroutine onespace

  subroutine assign_str_fr_struct_out(inu,s)
    integer :: inu
    type(supercell) :: s
    integer :: i

    do i = 1, 3
      read(inu,*) s%a(1:3,i)
    enddo
    read(inu,*) s%n
    do i = 1, s%n
      read(inu,*) s%at(i)%gr, s%at(i)%z,s%at(i)%f(1:3)
    enddo
    call real2recip(s%a,s%b)
    close(inu)
  end subroutine assign_str_fr_struct_out

  subroutine assign_str_fr_xsf(controlu,sc)
    integer :: i,controlu
    type(supercell) :: sc
    character(LEN=DSL) :: tmp

    read(controlu,*) tmp
    if(trim(tmp) /= 'CRYSTAL') then
      write(*,'(A)') 'err: Keyword CRYSTAL not found.'
      stop 1
    endif
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMVEC') then
      write(*,'(A)') 'err: Keyword PRIMVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    read(controlu,*) tmp
    if(trim(tmp) /= 'CONVVEC') then
      write(*,'(A)') 'err: Keyword CONVVEC not found.'
      stop 1
    endif

    read(controlu,*) sc%a(1:3,1)
    read(controlu,*) sc%a(1:3,2)
    read(controlu,*) sc%a(1:3,3)
    call real2recip(sc%a,sc%b)
    read(controlu,*) tmp
    if(trim(tmp) /= 'PRIMCOORD') then
      write(*,'(A)') 'err: Keyword PRIMCOORD not found.'
      stop 1
    endif
    read(controlu,*) sc%n
    do i = 1, sc%n
      read(controlu,*) sc%at(i)%z, sc%at(i)%ac
      call abs2frac(sc%b,sc%at(i)%ac,sc%at(i)%f)
    enddo
    sc%sg_label = '1-P1'
  end subroutine assign_str_fr_xsf

  subroutine assign_str_fr_poscar(controlu,sc)
    integer :: controlu
    type(supercell) :: sc
    real(double) :: totalscaling,vol
    integer :: i,j,s,nsp,v,z,ind
    character(LEN=DSL) :: comment1,comment2,tmp1,tmp2,valuestr,targetlabel
    real(double) :: pos(3)
    integer,allocatable :: numarr(:)
    integer :: totn,k,labellen,SGnumber,loc
    logical :: done,found,puredigits
    integer :: nvalid
    real(double) :: density

    read(controlu,DSF) comment1

    call onespace(DSL,comment1,comment2)

    sc%commentline = trim(comment2)
    write(*,'(A,A,A)') 'POSCAR comment is [',trim(comment2),'].'
    s = scan(comment2,' ')

    if(comment2(1:s-1) == 'ATMS:') then
      tmp1= comment2(s+1:DSL)

      read(tmp1,*) v
      sc%n = v
      if(sc%n < 1) then
        write(*,*) 'sc%n = ',sc%n
        write(*,'(A)') 'err: check sc%n in assign_str_fr_poscar'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)

      read(tmp1,*) v
      sc%nsp = v
      if(sc%nsp < 1) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,'(A)') 'err: check sc%nsp in assign_str_fr_poscar.'
        stop 1
      elseif(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp,zordmax=',sc%nsp,zordmax
        write(*,'(A)') 'err: check sc%nsp, i.e., the number of species/elements'
        stop 1
      endif
      s = scan(tmp1,' ')
      tmp1 = tmp1(s+1:DSL)
      nsp = v
      if(sc%nsp > zordmax) then
        write(*,*) 'sc%nsp = ',sc%nsp
        write(*,*) 'zordmax = ',zordmax
        write(*,'(A)') 'err: out of bound.'
        stop 1
      endif

      do i = 1, nsp

        read(tmp1,*) v

        sc%n_per_species(i) = v

        s = scan(tmp1,' ')
        tmp1 = tmp1(s+1:DSL)

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        found = .false.
        do j = 1, zordmax
          if(found) exit
          if(trim(valuestr) == trim(Ats(j)) .or. trim(valuestr) == trim(CAts(j)) ) then
            z = j
            found = .true.
          endif
        enddo
        if(.not. found) then
          write(*,*) 'valuestr = ',trim(valuestr)
          write(*,'(A)') 'err: this valuestr is not found in the periodic table.'
          stop 1
        endif
        sc%zarr(i) = z
        tmp1 = tmp1(s+1:DSL)
      enddo

      ind = 0
      do j = 1, sc%nsp
        ind = ind + sc%n_per_species(j)
      enddo

      if(ind /= sc%n) then
        write(*,*) 'ind,sc%n=',ind,sc%n
        write(*,'(A)') 'err: missing atoms in the assign_str_fr_poscar.'
        stop 1
      endif

      sc%sg_label = ''
      if(trim(tmp1) /= '') then

        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        if(trim(valuestr) /= 'SG') then
          write(*,*) 'valuestr is '//trim(valuestr)
          write(*,*) 'The first optional parameter must be SG.'
          write(*,'(A)') 'err: Check poscar.'
          stop 1
        endif
        tmp1 = tmp1(s+1:DSL)
        s = scan(tmp1,' ')
        valuestr = tmp1(1:s-1)
        sc%sg_label = trim(valuestr)
      else
        write(*,*) 'No SG label'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      allocate(numarr(sc%nsp))

      read(controlu,DSF) tmp1
      write(*,*) 'tmp1 is .'//trim(tmp1)//'.'
      call onespace(DSL,tmp1,tmp2)
      s = scan(tmp2,'0123456789')

      if(s == 0) then

        write(*,'(A)') 'We are going to process the species line |'//trim(tmp2)//'|'

        nsp = num_of_strings(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        do i = 1, sc%nsp
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)
          write(*,*) 'valuestr is '//trim(valuestr)
          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              z = j
              if(z /= sc%zarr(i)) then
                write(*,*) 'z, sc%zarr(i) = ',z,sc%zarr(i)
                write(*,'(A)') 'err: Inconsistency.'
                stop 1
              endif
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
        enddo

        read(controlu,DSF) tmp2
        nsp = num_of_integers(tmp2)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif

        read(tmp2,*) numarr(1:sc%nsp)
      else

        nsp = num_of_integers(tmp1)
        if(nsp /= sc%nsp) then
          write(*,*) 'err: nsp,sc%nsp = ',nsp,sc%nsp
          stop 1
        endif
        read(tmp1,*) numarr(1:sc%nsp)
      endif

      do i = 1, sc%nsp
        if(numarr(i) /= sc%n_per_species(i)) then
          write(*,*) 'species: i=',i
          write(*,*) 'numarr(i) = '//trim(N2str(numarr(i)))//', sc%n_per_species(i) = '//trim(N2str(sc%n_per_species(i)))
          write(*,'(A)') 'err: numarr(i) and sc%n_per_species(i) are not consistent'
          stop 1
        endif
      enddo
      deallocate(numarr)

    else
      if(comment2(1:s-1) == 'SG') then

        tmp2 = comment2(s+1:DSL)
        targetlabel = trim(tmp2)

        write(*,*) 'targetlabel = '//trim(targetlabel)

        puredigits = .true.
        labellen = len(trim(targetlabel))
        do i = 1, labellen
          if(.not. puredigits) then
            exit
          endif
         loc = scan(targetlabel(i:i),'0123456789')
         if(loc /= 1) then
           puredigits = .false.
         endif
        enddo

        if(puredigits)  then
          SGnumber = str2N(trim(targetlabel))
          write(*,*) 'SGnumber is ',SGnumber
          if(SGnumber >= 1 .and. SGnumber <= 230) then
            targetlabel = trim(SGbase(1,SGnumber)%sglabel)
            write(*,*) 'Fully expand the default SG label to the full SG label:'
            write(*,*) 'SGnumber= ',SGnumber
            write(*,*) 'Full targetlabel= ',trim(targetlabel)
          else
            write(*,*) 'err: must be between 1 and 230 inclusive. impossible. '
            stop 1
          endif
        endif
        sc%sg_label = trim(targetlabel)
        write(*,'(A)') 'Possibly adjusted sg_label is '//trim(sc%sg_label)
      else
        write(*,'(A)') 'WARNING WARNING: No ATMS:, no SG, hence we assume SG 1'
        sc%sg_label = '1-P1'
      endif

      read(controlu,*) totalscaling
      read(controlu,*) sc%a(1:3,1)
      read(controlu,*) sc%a(1:3,2)
      read(controlu,*) sc%a(1:3,3)

      read(controlu,DSF) tmp1
      call onespace(DSL,tmp1,tmp2)

      s = scan(tmp2,'0123456789')

      if(s == 0) then

        nsp = 0
        done = .false.
        do
          if(done) exit
          s = scan(tmp2,' ')
          valuestr = tmp2(1:s-1)

          found = .false.
          do j = 1, zordmax
            if(found) exit
            if(trim(valuestr) == trim(CATS(j))) then
              nsp = nsp + 1
              sc%zarr(nsp) = j
              found = .true.
            endif
          enddo
          if( .not. found) then
            write(*,'(A)') 'err: not found.'
            stop 1
          endif
          tmp2 = tmp2(s+1:DSL)
          if(len(trim(tmp2)) == 0) then
            write(*,'(A)') 'nsp = '//trim(N2str(nsp))
            sc%nsp = nsp
            done = .true.
          endif
        enddo

        read(controlu,DSF) tmp1

        nvalid = num_of_integers(tmp1)
        if(nsp /= nvalid) then
          write(*,*) 'nsp,nvalid=',nsp,nvalid
          write(*,*) 'err: inconsistency in species line and species line.'
          stop 1
        endif
        read(tmp1,*) sc%n_per_species(1:nsp)

      else
        write(*,'(A)') 'err: Since this is without ATMs, we must have species line.'
        stop 1
      endif
    endif

    totn = 0
    do j = 1, sc%nsp
      v = sc%n_per_species(j)
      do k = 1, v
        totn = totn + 1
        sc%at(totn)%z = sc%zarr(j)
        sc%at(totn)%mass = massofa(sc%zarr(j))
      enddo
    enddo
    write(*,'(A)') 'totn = '//trim(N2str(totn))
    sc%n = totn

    sc%a = sc%a*totalscaling

    call real2recip(sc%a,sc%b)

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Real latt:', sc%a(1:3,1),'    | ',vecmag3(sc%a(1:3,1)),sc%a(1:3,2),'    | ',vecmag3(sc%a(1:3,2)),sc%a(1:3,3),'    | ',vecmag3(sc%a(1:3,3))
    vol = det3(sc%a(1,1))
    if(vol .le. 0) then
      write(*,'(A)') 'err: vol is not positive.'
      stop 1
    endif

    sc%vol = vol

    write(*,'(A,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10,/,3F15.10,A,F15.10)') &
      'Reciprocal latt:', sc%b(1:3,1),'    | ',vecmag3(sc%b(1:3,1)),sc%b(1:3,2),'    | ',vecmag3(sc%b(1:3,2)),sc%b(1:3,3),'    | ',vecmag3(sc%b(1:3,3))

    read(controlu,*) tmp1
    if(trim(tmp1) == 'Selective') then
      read(controlu,*) tmp1
    endif
    if(trim(tmp1) /= 'Direct' .and. trim(tmp1) /= 'direct' .and. trim(tmp1) /= 'Cartesian' ) then
      write(*,'(A)') 'err: dummy should be Direct, direct, or Cartesian'
      stop 1
    endif

    do i = 1, sc%n
      sc%at(i)%force(1:3) = (/zero,zero,zero/)
    enddo
    do i = 1, sc%n
      if(trim(tmp1) == 'Direct' .or. trim(tmp1) == 'direct') then
        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) sc%at(i)%f(1:3)

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      else if(trim(tmp1) == 'Cartesian') then

        read(controlu,DSF) comment1
        call onespace(DSL,comment1,comment2)
        read(comment2,*) pos(1:3)

        call abs2frac(sc%b(1,1),pos(1),sc%at(i)%f(1))

        call frac2abs(sc%a(1,1),sc%at(i)%f(1),sc%at(i)%ac(1))
        s = scan(comment2,'|')
        if(s > 0) then
          tmp2= comment2(s+1:DSL)

          read(tmp2,*) sc%at(i)%force(1:3)
          write(*,'(A,3F20.10)') 'sc%at(i)%force(1:3) = ',sc%at(i)%force(1:3)
        endif
      endif
    enddo
    call getlenang(sc%a,sc%la(1))
    write(*,'(A)') 'la(1:6)='
    do i = 1, 6
      write(*,'(3F15.10)') sc%la(i)
    enddo
    write(*,*) 'Volume of the supercell= ',vol
    density = crys_density(sc)
    sc%density = density
    write(*,*) 'Density is ',density,' g/cm^3'

  end subroutine assign_str_fr_poscar

  subroutine assign_str_fr_fdf(inu,s)
    integer :: inu,natom,nspecies,sp,i,j
    type(supercell) :: s
    character(len=DSL) :: str1,str2,str3
    real(double) :: x(3)
    character(len=DSL) :: icf

    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfAtoms') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofAtoms'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') natom
      write(*,*) 'natom = ',natom
      s%n = natom
    endif
    read(inu,*) str1,str2
    if(trim(str1) /= 'NumberOfSpecies') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be NumberofSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    else
      read(str2,'(I8)') nspecies
      write(*,*) 'nspecies = ',nspecies
      s%nsp = nspecies
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, nspecies
      read(inu,*) j,sp
      if(j/=i) then
        write(*,*) 'j,i=',j,i
        write(*,*) 'sp = ',sp
        write(*,'(A)') 'err: something wrong with the species numbering ?'
        stop 1
      else
        s%zarr(i) = sp
        write(*,*) 's%zarr(',i,')=',sp
      endif
    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'ChemicalSpeciesLabel') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be ChemicalSpeciesLabel'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2,str3
    if(trim(str1) /= 'LatticeConstant') then
      write(*,*) 'str1 = ',trim(str1)
      write(*,*) 'but str1 must be LatticeConstant'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str2) /= '1.0') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be 1.0'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    if(trim(str3) /= 'Ang') then
      write(*,*) 'str3 = ',trim(str3)
      write(*,*) 'but str3 must be Ang'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be  LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) s%a(1:3,1)
    read(inu,*) s%a(1:3,2)
    read(inu,*) s%a(1:3,3)
    write(*,*) 's%a is '
    write(*,*) s%a(1:3,1)
    write(*,*) s%a(1:3,2)
    write(*,*) s%a(1:3,3)
    call real2recip(s%a,s%b)
    read(inu,*) str1,str2
    if(trim(str2) /= 'LatticeVectors') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be LatticeVectors'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) == 'Fractional') then
      icf = 'frac'
    else if(trim(str2) == 'NotScaledCartesianAng') then
      icf = 'abs'
    else
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be either Fractional or NotScaledCartesianAng'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
    do i = 1, natom
      read(inu,*) x(1:3),sp
      if(trim(icf) == 'frac') then
        s%at(i)%f(1:3) = x(1:3)
      else if(trim(icf) == 'abs') then
        call abs2frac(s%b,x,s%at(i)%f)
      else
        write(*,'(A)') 'err: frac and abs problem.'
        stop 1
      endif
      s%at(i)%z = s%zarr(sp)

    enddo
    read(inu,*) str1,str2
    if(trim(str2) /= 'AtomicCoordinatesAndAtomicSpecies') then
      write(*,*) 'str2 = ',trim(str2)
      write(*,*) 'but str2 must be AtomicCoordinatesAndAtomicSpecies'
      write(*,'(A)') 'err: check fdf file.'
      stop 1
    endif
  end subroutine assign_str_fr_fdf

  function massofa(z)
    integer :: z
    real(double) :: massofa
    if(z < 0 .or. z > zordmax) then
      write(*,*) 'z = ',z
      write(*,'(A)') 'massofa not defined yet.'
      write(*,*) 'err: massofa problem'
      stop 1
    endif
    massofa = MASSOFATOM(z)
  end function massofa

  function crys_density(s)
    type(supercell) :: s
    real(double) :: mass,crys_density
    integer :: i

    mass = 0.0d0
    do i = 1, s%n
      mass = mass + massofatom(s%at(i)%z)
    enddo
    crys_density = (mass*AMU/(det3(s%a(1:3,1:3))*1.d-30))*(1.0d3/1.0d6)
  end function crys_density

  subroutine read_struc(inputu,inputfile,filestem,inputformat,s)
    integer :: inputu,length,ind
    character(len=*) :: inputfile,filestem,inputformat
    type(supercell) :: s

    length = len_trim(inputfile)
    ind = index(trim(inputfile),'.',BACK=.TRUE.)
    filestem=inputfile(1:ind-1)

    inputformat=inputfile(ind+1:length)

    open(inputu,file=trim(inputfile),status='old',action='read')

    if(trim(inputformat)=='vasp' .or. trim(inputformat) == 'VASP') then
      call assign_str_fr_poscar(inputu,s)

    else if(trim(inputformat)=='STRUCT_OUT') then
      call assign_str_fr_struct_out(inputu,s)
    else if(trim(inputformat)=='fdf') then
      call assign_str_fr_fdf(inputu,s)
    else if(trim(inputformat)=='xsf') then
      call assign_str_fr_xsf(inputu,s)
    else
      write(*,*)
      write(*,*) 'WARNING: accepted format are vasp,arc/car,STRUCT_OUT,gjf,fdf'
      write(*,*) 'but inputformat is ',trim(inputformat)
      write(*,*) 'unrecognized input file format.'
      write(*,'(A)') 'err: check input file format.'
      stop 1
    endif
    close(inputu)

  end subroutine read_struc
end module commod
