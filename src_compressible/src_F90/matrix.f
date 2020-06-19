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

c$$$      ! this is a zero-based array routine for easier conversion to c++
c$$$      subroutine decomp(neq,ndim,a,ip)
c$$$      implicit none
c$$$      integer neq,ndim
c$$$      real*8  a(0:ndim-1,0:ndim-1)
c$$$      integer ip(0:ndim-1)
c$$$      real*8 eps,eps2,amult,apm,apn
c$$$      integer n,nm1,i,ip1,k,j,l
c$$$
c$$$      eps = 1.d-30
c$$$      eps2 = 1.d-60
c$$$      
c$$$      nm1=neq-1
c$$$      do n=0,neq-1
c$$$         ip(n)=n+1
c$$$      end do
c$$$      do i=0,nm1-1
c$$$         ip1=i+2
c$$$         k=i+1
c$$$         apm    =abs(a(ip(i)-1,i))
c$$$         do j=ip1-1,neq-1
c$$$            apn=abs(a(ip(j)-1,i))
c$$$            if(apm.ge.apn) cycle
c$$$            apm=apn  
c$$$            k=j+1
c$$$         end do
c$$$         j=ip(k-1)
c$$$         ip(k-1)=ip(i)
c$$$         ip(i)=j
c$$$         if(apm.lt.eps) then
c$$$            write(6,*)a(ip(neq-1)-1,neq-1)
c$$$            stop"singular matrix in decomp"
c$$$            return
c$$$         end if
c$$$         do l=ip1-1,neq-1
c$$$            n=ip(l)
c$$$            amult=a(n-1,i)/a(j-1,i)
c$$$            a(n-1,i)=amult
c$$$            if(abs(amult).lt.eps2) cycle
c$$$            do k=ip1-1,neq-1
c$$$               a(n-1,k)=a(n-1,k)-amult*a(j-1,k)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      if(abs(a(ip(neq-1)-1,neq-1)).lt.eps) then
c$$$         write(6,*)a(ip(neq-1)-1,neq-1)
c$$$         stop"singular matrix in decomp"
c$$$      end if
c$$$      
c$$$      return
c$$$      end


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


c$$$      ! this is a zero-based array routine for easier conversion to c++
c$$$      subroutine solve(neq,ndim,a,b,ip)
c$$$      implicit none
c$$$      real*8 a(0:ndim-1,0:ndim-1),b(0:ndim-1)
c$$$      integer ip(0:ndim-1)
c$$$      integer neq,ndim
c$$$      integer nm1,l,lm1,k,j,jp1,n
c$$$      real*8 scr
c$$$      
c$$$      nm1 = neq-1
c$$$      do l=1,neq-1
c$$$         n=ip(l)
c$$$         lm1=l
c$$$         do k=0,lm1-1
c$$$            b(n-1)=b(n-1)-a(n-1,k)*b(ip(k)-1)
c$$$         end do
c$$$      end do
c$$$      b(ip(neq-1)-1)=b(ip(neq-1)-1)/a(ip(neq-1)-1,neq-1)
c$$$      do l=0,nm1-1
c$$$         j=neq-(l+1)
c$$$         jp1=j+1
c$$$         n=ip(j-1)
c$$$         do k=jp1-1,neq-1
c$$$            b(n-1)=b(n-1)-a(n-1,k)*b(ip(k)-1)
c$$$         end do
c$$$         b(n-1)=b(n-1)/a(n-1,j-1)
c$$$      end do
c$$$      do n=0,neq-1
c$$$         do while (ip(n) .ne. n+1)
c$$$            j=ip(n)
c$$$            scr=b(j-1)
c$$$            ip(n)=ip(j-1)
c$$$            b(j-1)=b(ip(j-1)-1)
c$$$            b(ip(j-1)-1)=scr
c$$$            ip(j-1)=j
c$$$         end do
c$$$      end do
c$$$         
c$$$      return
c$$$      end
