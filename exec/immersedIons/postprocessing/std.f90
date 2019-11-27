program std

  implicit none

  character(len=80) filename

  integer i, j, npts, nfiles
  
  double precision, allocatable :: file(:)

  double precision :: temp1, temp2
  
  print*,'enter filename'
  read(*,*) filename

  print*,'how many points'
  read(*,*) npts

  allocate(file(npts))

  open(unit=999,file=trim(filename))
  do i=1,npts
     read(unit=999, FMT=*) file(i)
  end do

  temp1 = 0.
  do i=1,npts
     temp1 = temp1 + file(i) 
  end do
  temp1 = temp1 / npts

  temp2 = 0.
  do i=1,npts
     temp2 = temp2 + ( file(i) - temp1 )**2
  end do
  temp2 = sqrt(temp2/npts)

  print*,'mean',temp1
  print*,'standard deviation',temp2
  print*,'standard error',temp2/sqrt(dble(npts))

end program std
