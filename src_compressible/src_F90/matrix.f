      subroutine decomp(neq,ndim,a,ip)
      implicit none
      integer neq,ndim
      real*8  a(ndim,ndim)
      integer ip(ndim)
      real*8 eps,eps2,amult,apm,apn
      integer n,nm1,i,ip1,k,j,l

      eps = 1.d-30
      eps2 = 1.d-60
      
      nm1=neq-1
      do n=1,neq
         ip(n)=n
      end do
      do i=1,nm1
         ip1=i+1  
         k=i
         apm    =abs(a(ip(i),i))
         do j=ip1,neq
            apn=abs(a(ip(j),i))
            if(apm.ge.apn) cycle
            apm=apn  
            k=j
         end do
         j=ip(k)
         ip(k)=ip(i)
         ip(i)=j
         if(apm.lt.eps) then
            write(6,*)a(ip(neq),neq)
            stop"singular matrix in decomp"
            return
         end if
         do l=ip1,neq
            n=ip(l)
            amult=a(n,i)/a(j,i)
            a(n,i)=amult
            if(abs(amult).lt.eps2) cycle
            do k=ip1,neq
               a(n,k)=a(n,k)-amult*a(j,k)
            end do
         end do
      end do
      if(abs(a(ip(neq),neq)).lt.eps) then
         write(6,*)a(ip(neq),neq)
         stop"singular matrix in decomp"
      end if
      
      return
      end

      subroutine solve(neq,ndim,a,b,ip)
      implicit none
      real*8 a(ndim,ndim),b(ndim)
      integer ip(ndim)
      integer neq,ndim
      integer nm1,l,lm1,k,j,jp1,n
      real*8 scr
      
      nm1 = neq-1
      do l=2,neq
         n=ip(l)
         lm1=l-1
         do k=1,lm1
            b(n)=b(n)-a(n,k)*b(ip(k))
         end do
      end do
      b(ip(neq))=b(ip(neq))/a(ip(neq),neq)
      do l=1,nm1
         j=neq-l
         jp1=j+1
         n=ip(j)
         do k=jp1,neq
            b(n)=b(n)-a(n,k  )*b(ip(k))
         end do
         b(n)=b(n)/a(n,j)
      end do
      do n=1,neq
         do while (ip(n) .ne. n)
            j=ip(n)
            scr=b(j)
            ip(n)=ip(j)
            b(j)=b(ip(j))
            b(ip(j))=scr
            ip(j)=j
         end do
      end do
         
      return
      end
