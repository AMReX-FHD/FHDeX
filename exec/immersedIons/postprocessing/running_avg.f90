program running_avg

  implicit none

  character(len=80) filename

  integer i, npts
  
  double precision sum

  double precision, allocatable :: data(:)
  
  print*,'enter filename'
  read(*,*) filename

  print*,'how many total points'
  read(*,*) npts
  
  allocate(data(npts))

  open(unit=999,file=trim(filename))
  do i=1,npts
     read(unit=999, FMT=*) data(i)
  end do

  sum = 0
  
  do i=1,npts
     sum = sum + data(i)
     print*,sum/i
  end do

end program running_avg
