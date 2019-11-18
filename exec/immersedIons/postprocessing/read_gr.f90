program read_gr

  implicit none

  character(len=80) filename

  integer i, j, npts, nfiles
  
  double precision, allocatable :: gr(:,:), gr_avg(:,:)
  
  print*,'enter filename'
  read(*,*) filename

  print*,'how many points per sub-file'
  read(*,*) npts

  print*,'how many files'
  read(*,*) nfiles

  allocate(gr(npts*nfiles,5))
  allocate(gr_avg(npts,5))

  open(unit=999,file=trim(filename))
  do i=1,npts*nfiles
     read(unit=999, FMT=*) gr(i,1), gr(i,2), gr(i,3), gr(i,4), gr(i,5)
  end do

  gr_avg(:,:) = 0.
  do i=1,nfiles
     do j=1,npts
        gr_avg(j,1) = gr_avg(j,1) + gr( (i-1)*npts + j, 1 )
        gr_avg(j,2) = gr_avg(j,2) + gr( (i-1)*npts + j, 2 )
        gr_avg(j,3) = gr_avg(j,3) + gr( (i-1)*npts + j, 3 )
        gr_avg(j,4) = gr_avg(j,4) + gr( (i-1)*npts + j, 4 )
        gr_avg(j,5) = gr_avg(j,5) + gr( (i-1)*npts + j, 5 )
     end do
  end do

  gr_avg = gr_avg / nfiles


  do j=1,npts
     print*,gr_avg(j,1:5)
  end do

end program read_gr
