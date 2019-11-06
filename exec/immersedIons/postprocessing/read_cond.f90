program read_cond

  implicit none

  character(len=80) filename

  integer i, j, npts, interval, nmpi, nintervals, n_intervals_skip
  
  double precision x, y, cond_temp, average

  double precision, allocatable :: cond(:)
  
  print*,'enter filename'
  read(*,*) filename

  print*,'how many total points'
  read(*,*) npts

  print*,'stats reset interval'
  read(*,*) interval

  print*,'how many sample to throw out?'
  read(*,*) n_intervals_skip
  
  allocate(cond(npts))

  open(unit=999,file=trim(filename))
  do i=1,npts
     read(unit=999, FMT=*) cond(i)
  end do

  nintervals = npts / (interval)
  average = 0.
  
  do j=1,nintervals
     cond_temp = cond(j*interval)
     print*,cond_temp
     if (j > n_intervals_skip) then
        average = average + cond_temp
     end if
  end do

  average = average / (nintervals-n_intervals_skip)
  print*,'overall average',average

end program read_cond
