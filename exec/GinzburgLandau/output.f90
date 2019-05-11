       subroutine output(ist,iend,jst,jend,phi,time,ndim)

       implicit none

       integer ist,iend,jst,jend,ndim
       double precision time

       double precision phi(0:ndim+1,0:ndim+1)

       integer i,j

       write(6,*)" solution at  time = ", time

       do j=jst,jend


       write(6,*)j
       write(6,1000) (phi(i,j),i=ist,iend)
1000   format(5e15.7)

       enddo

       return
       end
