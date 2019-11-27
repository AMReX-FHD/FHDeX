program read_curr

  implicit none

  character(len=80) filename

  integer i, j, npts, interval, nmpi, nintervals, n_intervals_skip
  
  double precision x, y, curr_temp, average

  double precision, allocatable :: curr(:)
  double precision :: a, b, c, d, e
  
  print*,'enter filename'
  read(*,*) filename

  print*,'how many total points'
  read(*,*) npts

  print*,'stats reset interval'
  read(*,*) interval

  print*,'how many sample to throw out?'
  read(*,*) n_intervals_skip
  
  allocate(curr(npts))

  open(unit=999,file=trim(filename))
  do i=1,npts
     read(unit=999, FMT=*) curr(i), a, b, c, d, e
  end do

  nintervals = npts / (interval)
  average = 0.
  
  do j=n_intervals_skip+1,nintervals
     curr_temp = curr(j*interval-1)
     print*,curr_temp
     average = average + curr_temp
  end do

  average = average / (nintervals-n_intervals_skip)
  print*,'overall average',average

end program read_curr
