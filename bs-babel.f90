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

module extra
  use commod
  implicit none
  integer,parameter :: inputu=10
  integer,parameter :: kptu=11
  integer,parameter :: eigenvu=12
  integer,parameter :: bsu=13
  integer,parameter :: anchoru=14
  integer,parameter :: zerou=15

  logical :: vasp,siesta,qphonon

  type bst
    integer :: nspin
    integer :: ndy
    integer :: nsegments
    integer :: maxn
    real(double) :: elowerb(2),eupperb(2)
    integer :: commdline_nbands_print
    integer :: actual_nbands_print

    integer,allocatable :: npoints(:)
    real(double),allocatable :: bvec(:,:),evec(:,:)
    real(double),allocatable :: q(:,:,:)
    real(double),allocatable :: xdis(:,:)
    real(double),allocatable :: anchorxdis(:)
    real(double),allocatable :: freq(:,:,:,:)
    real(double),allocatable :: r(:)
    integer,allocatable :: ind(:)

  end type bst

contains

  subroutine get_branches_information(bs,dir,s)
    type(bst) :: bs
    type(supercell) :: s
    character(len=*) :: dir
    character(len=DSL) :: kptfile,line
    integer :: nsegments,countn
    real(double) :: f(3),startv(3),endv(3)
    integer :: n,i,maxn,nintervals,j
    integer :: ios1,ios2
    real(double) :: p,x(3),curpos,branchdis,delta

    if(vasp) then
      kptfile=trim(dir)//'/KPOINTS'
    else if(siesta) then
      kptfile=trim(dir)//'/writebands.fdf'
    else if(qphonon) then
      kptfile=trim(dir)//'/vasp_line_mode_kpoints.dat'
    else
      write(*,*) 'You must provide a kpoint file.'
      stop
    endif

    open(unit=kptu,file=trim(kptfile),action='read')
    write(*,'(A)') 'kptfile= '//trim(kptfile)//' is opened'
    if(vasp .or. qphonon) then
      read(kptu,*) line
      read(kptu,*) nintervals
      read(kptu,*) line
      if(trim(line) /= 'Line-mode') then
        write(*,*) 'line is '//trim(line)
        stop 'err: must be Line-mode'
      endif
      read(kptu,*) line
      if(trim(line) /= 'rec') then
        write(*,*) 'line is '//trim(line)
        stop 'err: must be rec.'
      endif

      i = 0
      do
        read(kptu,*,iostat=ios1) startv(1:3)
        if(ios1 /= 0) then
          write(*,*) 'end of pair reached'
          nsegments=i
          write(*,*) 'the number of branches is ',nsegments
          exit
        endif
        read(kptu,*,iostat=ios2) endv(1:3)
        if(ios2 /= 0) then
          write(*,*) 'endv not defined.'
          stop 'err: check KPOINTS'
        endif
        i =  i + 1
      enddo

      close(kptu)
      bs%nsegments = nsegments
      allocate(bs%npoints(nsegments))
      bs%npoints(1:nsegments) = nintervals
      allocate(bs%bvec(3,nsegments),bs%evec(3,nsegments))

      open(unit=kptu,file=trim(kptfile),action='read')
      do i = 1, 4
        read(kptu,*) line
      enddo
      do i = 1, nsegments
        read(kptu,*) startv(1:3)
        read(kptu,*) endv(1:3)
        bs%bvec(1:3,i) = startv(1:3)
        bs%evec(1:3,i) = endv(1:3)
      enddo
    else if(siesta) then
      read(kptu,'(A17)') line
      if(line(1:17) /= 'WriteBands .true.') then
        write(*,'(A17)') line(1:17)
        stop 'err: must be WriteBands .true.'
      endif
      read(kptu,'(A39)') line
      if(line(1:39) /= 'BandLinesScale ReciprocalLatticeVectors') then
      write(*,'(A39)') line(1:39)
        stop 'err: must be BandLinesScale ReciprocalLatticeVectors'
      endif
      read(kptu,*) line
      if(trim(line) /= '%block') then
        write(*,'(A)') trim(line)
        stop 'err: must be %block'
      endif

      countn=0
      do
        read(kptu,*) line
        if(trim(line)=='%endblock') then
          exit
        endif
        countn = countn + 1
      enddo
      close(kptu)
      nsegments = countn-1
      write(*,*) 'nsegments = ',nsegments
      bs%nsegments = nsegments
      allocate(bs%npoints(nsegments))
      allocate(bs%bvec(3,nsegments),bs%evec(3,nsegments))

      open(unit=kptu,file=trim(kptfile),action='read')
      do i = 1, 3
        read(kptu,*) line
      enddo

      read(kptu,*) n, f(1:3)
      if(n /= 1) stop 'err: n must be 1'
      bs%bvec(1:3,1) = f(1:3)
      do i = 1, nsegments
        read(kptu,*) n,f(1:3)
        bs%npoints(i) = n+1
        bs%evec(1:3,i) = f(1:3)
        if(i /= nsegments) then
          bs%bvec(1:3,i+1) = f(1:3)
        endif
      enddo
    else
       write(*,*)
       stop 'err: later.'
    endif
    close(kptu)

    do i = 1, nsegments
      write(*,'(I6,6F10.5)') bs%npoints(i),bs%bvec(1:3,i),bs%evec(1:3,i)
    enddo
    maxn = -10000
    do i = 1, nsegments
      if(bs%npoints(i) > maxn) then
        maxn = bs%npoints(i)
      endif
    enddo
    if(maxn < 1) then
      write(*,*) 'err: maxn = ',maxn
      stop 1
    endif
    write(*,*) 'maxn = ',maxn
    bs%maxn=maxn

    allocate(bs%q(3,maxn,nsegments))
    allocate(bs%xdis(maxn,nsegments))
    allocate(bs%anchorxdis(nsegments+1))

    curpos = zero
    bs%anchorxdis(1) = curpos
    do i = 1, bs%nsegments
      call frac2abs(s%b,bs%bvec(1:3,i),startv)
      call frac2abs(s%b,bs%evec(1:3,i),endv)
      x = endv - startv
      branchdis = vecmag3(x)
      delta = branchdis/(bs%npoints(i)-one)
      do j = 1, bs%npoints(i)
        p = curpos + (j-one)*delta
        bs%xdis(j,i) = p
      enddo
      curpos = curpos + branchdis
      bs%anchorxdis(i+1) = curpos
      do j = 1, bs%npoints(i)
        bs%q(1:3,j,i) = bs%bvec(1:3,i) + (j-one)*(bs%evec(1:3,i)-bs%bvec(1:3,i))/(bs%npoints(i)-one)
      enddo
    enddo
  end subroutine get_branches_information

  subroutine energy_for_all_branches(bs,dir,SystemLabel,efermi)
    type(bst) :: bs
    character(len=*) :: dir,SystemLabel
    character(len=300) :: eigenvaluefile,line
    integer :: countn,a,b,nbands,i,nspin,j,k,p,nkpt
    real(double) :: v,x(3),ener_freq,efermi

    if(vasp) then
      eigenvaluefile=trim(dir)//'/EIGENVAL'
    else if(siesta) then

      eigenvaluefile=trim(dir)//'/'//trim(SystemLabel)//'.bands'
    else if(qphonon) then
      eigenvaluefile=trim(dir)//'/tbEIGENVAL'
    else
    endif
    open(unit=eigenvu,file=trim(eigenvaluefile),action='read')
    write(*,'(A)') 'eigenvaluefile '//trim(eigenvaluefile)//' is opened'
    bs%elowerb(1:2) = 1d100
    bs%eupperb(1:2) = -1d100

    if(vasp .or. qphonon) then
      nspin=1
      bs%nspin = 1
      do i = 1, 5
        read(eigenvu,*) line
      enddo

      read(eigenvu,*) a,b,nbands
      if(b /= bs%nsegments*bs%npoints(1)) then
        write(*,*) 'b,nsegments,nintervals=',b,bs%nsegments,bs%npoints(1)
        write(*,*) 'b must be equal to nsegments*nintervals, now: b,nsegments*nintervals=',b,bs%nsegments*bs%npoints(1)
        stop 'err: check number of kpoints.'
      endif

      bs%ndy = nbands
      allocate(bs%freq(nbands,bs%maxn,bs%nsegments,nspin))
      do i = 1, bs%nsegments
        do j = 1, bs%npoints(i)
          read(eigenvu,*) x(1:3)
          if(vecmag3(bs%q(1:3,j,i)-x(1:3)) > 1d-7) then
            write(*,*) 'bs%q(1:3,j,i),x(1:3)=',bs%q(1:3,j,i),x(1:3)
            stop 'err: bs%q and x are not the same.'
          endif
          do k = 1, nbands
            read(eigenvu,*) a,ener_freq
            if(a /= k) then
              write(*,*) 'a,k=',a,k
              stop 'err: a must be equal to k'
            endif
            ener_freq = ener_freq - efermi
            bs%freq(k,j,i,nspin) = ener_freq
            if(ener_freq < bs%elowerb(nspin)) then
              bs%elowerb(nspin) = ener_freq
            endif
            if(ener_freq > bs%eupperb(nspin)) then
              bs%eupperb(nspin) = ener_freq
            endif
          enddo
        enddo
      enddo

    else if(siesta) then

      read(eigenvu,*) line
      write(*,*) 'line is ',trim(line)
      read(eigenvu,*) line
      write(*,*) 'line is ',trim(line)
      read(eigenvu,*) line
      write(*,*) 'line is ',trim(line)
      read(eigenvu,*) nbands,nspin,nkpt
      write(*,*) 'nbands,nspin,nkpt=',nbands,nspin,nkpt
      bs%ndy = nbands
      bs%nspin= nspin
      write(*,*) 'bs%maxn,bs%nsegments=',bs%maxn,bs%nsegments
      allocate(bs%freq(nbands,bs%maxn,bs%nsegments,nspin))
      countn=0
      do i = 1, bs%nsegments
        if(i == 1) then
          do j = 1, bs%npoints(i)
            countn = countn + 1
            read(eigenvu,*) v,bs%freq(1:nbands,j,i,1:nspin)
          enddo
        else
          if(nspin > 1) then

          endif

          bs%freq(1:nbands,1,i,1:nspin) = bs%freq(1:nbands,bs%npoints(i-1),i-1,1:nspin)
          do j = 2, bs%npoints(i)
            countn = countn + 1
            read(eigenvu,*) v,bs%freq(1:nbands,j,i,1:nspin)
          enddo
        endif
      enddo
      if(countn /= nkpt) then
        write(*,*) 'countn,nkpt=',countn,nkpt
      endif
      do i = 1, bs%nsegments
        do j = 1, bs%npoints(i)
          do k = 1, bs%ndy
            do p = 1, nspin
              v = bs%freq(k,j,i,p)

              v = v-efermi
              bs%freq(k,j,i,p) = v
              if(v < bs%elowerb(p)) then
                bs%elowerb(p) = v
              endif
              if(v > bs%eupperb(p)) then
                bs%eupperb(p) = v
              endif
            enddo
          enddo
        enddo
      enddo
      write(*,*) 'bs%elowerb(1:nspin),bs%eupperb(1:nspin)=',bs%elowerb(1:nspin),bs%eupperb(1:nspin)
    endif
  end subroutine energy_for_all_branches

  subroutine output_bandstructure(bs)
    type(bst) :: bs
    integer :: i,j,k
    character(len=DSL) :: filename

    if(bs%commdline_nbands_print == 0) then
      bs%actual_nbands_print = bs%ndy
    else
      bs%actual_nbands_print = bs%commdline_nbands_print
    endif
    write(*,*) 'to print ',bs%actual_nbands_print , ' bands.'
    do k = 1, bs%nspin
      filename=trim(N2str(k))//'-bs.dat'
      open(unit=bsu,file=trim(filename),status='replace')
      write(bsu,'(A)') '# '//trim(N2str(bs%nsegments))//' '//trim(N2str(bs%npoints(1)))//'   '//trim(N2str(bs%actual_nbands_print))
      do i = 1, bs%nsegments
        do j = 1, bs%npoints(i)
          write(bsu,'(F15.8,1000F20.10)') bs%xdis(j,i),bs%freq(1:bs%actual_nbands_print,j,i,k)
        enddo
      enddo
      close(bsu)
      write(*,*) trim(filename)//' is created.'
    enddo
  end subroutine output_bandstructure

  subroutine output_anchor_points(bs)
    type(bst) :: bs
    real(double) maxv,minv
    integer :: i,j
    character(len=DSL) :: filename

    real(double),parameter :: shift=50.0d0

    do j = 1, bs%nspin
      maxv = bs%eupperb(j)+shift
      minv = bs%elowerb(j)-shift
      filename=trim(N2str(j))//'-anchor.dat'
      open(unit=anchoru,file=trim(filename),status='replace')
      do i = 1, bs%nsegments+1
        if(mod(i,2) == 0) then
          write(anchoru,*) bs%anchorxdis(i),maxv,minv
          write(anchoru,*) bs%anchorxdis(i),minv,maxv
        else
          write(anchoru,*) bs%anchorxdis(i),minv,maxv
          write(anchoru,*) bs%anchorxdis(i),maxv,minv
        endif
      enddo
      close(anchoru)
      write(*,*) trim(filename)//' is created'
      open(unit=zerou,file='zero.dat',status='replace')
      write(zerou,*) 0.d0, 0.0d0
      write(zerou,*) bs%anchorxdis(bs%nsegments+1),0.0d0
      close(zerou)
      write(*,*) 'plotbs or xmgrace -nxy '//trim(N2str(j))//'-bs.dat '//trim(N2str(j))//'-anchor.dat zero.dat'
    end do
    if(bs%nspin == 2) then
      write(*,'(A)') 'You may want to do xmgrace -nxy 1-bs.dat 1-anchor.dat zero.dat -nxy 2-bs.dat 2-anchor.dat zero.dat'
    endif
  end subroutine output_anchor_points
end module extra

program sample
  use extra
  implicit none
  character(len=300) :: cmd,strucfile,dir,filestem
  character(len=DSL) :: codename,inputformat,SystemLabel
  type(supercell) :: s
  real(double) :: efermi
  type(bst) :: bs

  write(*,'(///A)') 'b.out (1) vasp,siesta,... (2) working directory (3) SystemLabel (any dummy label when vasp bs is to be plotted) (4) efermi (5) commdline_nbands_print, 0 means all bands are printed'
  write(*,*) 'version 2010-04-05'

  call fargn(5)
  call fargv(1,codename)
  call fargv(2,dir)
  call fargv(3,SystemLabel)
  call fargv(4,efermi)
  call fargv(5,bs%commdline_nbands_print)

  if(bs%commdline_nbands_print == 0) then
    write(*,*) 'Going to print all bands eigenvalues'
  else
    write(*,*) 'Going to print ', bs%commdline_nbands_print, ' bands only.'
  endif

  siesta=.FALSE.
  vasp=.FALSE.
  qphonon=.FALSE.
  if(trim(codename) == 'vasp') then
    vasp=.TRUE.
  else if(trim(codename) == 'siesta') then
    siesta=.TRUE.
  else if(trim(codename) == 'qphonon') then
    qphonon=.TRUE.
  else
    write(*,'(A)') 'codename is '//trim(codename)//' not supported.'
  endif

  if(vasp) then
    cmd='cp -f '//trim(dir)//'/POSCAR '//trim(dir)//'/bs.vasp'
    write(*,'(A)') 'cmd is '//trim(cmd)
    call system(trim(cmd))
    strucfile=trim(dir)//'/bs.vasp'
  else if(siesta) then
    strucfile=trim(dir)//'/poscar.fdf'
  else if(qphonon) then
    strucfile=trim(dir)//'/p.vasp'
  else
    write(*,*) 'missing structure file ?'
    stop 'err: check codename'
  endif
  call read_struc(inputu,trim(strucfile),filestem,inputformat,s)
  write(*,*) 'reciprocal lattice vector information:'
  write(*,'(3F10.5)') s%b(1:3,1)
  write(*,'(3F10.5)') s%b(1:3,2)
  write(*,'(3F10.5)') s%b(1:3,2)
  call get_branches_information(bs,dir,s)
  call energy_for_all_branches(bs,dir,SystemLabel,efermi)
  call output_bandstructure(bs)
  call output_anchor_points(bs)
end program sample
