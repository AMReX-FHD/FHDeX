      subroutine decompnp(neq,ndim,a)
      implicit none
      integer neq,ndim
      real*8  a(ndim,ndim)
      real*8 eps,eps2,amult,apm,apn
      integer n,nm1,i,ip1,k,j,l

      eps = 1.d-30
      eps2 = 1.d-60
      
      nm1=neq-1
      do i=1,nm1
         ip1 = i+1
         do l=ip1,neq
            amult=a(l,i)/a(i,i)
            a(l,i)=amult
            if(abs(amult).lt.eps2) cycle
            do k=ip1,neq
               a(l,k)=a(l,k)-amult*a(i,k)
            end do
         end do
      end do
      if(abs(a(neq,neq)).lt.eps) then
         write(6,*)a(neq,neq)
         stop"singular matrix in decompnp"
      end if
      
      return
      end

      subroutine solvenp(neq,ndim,a,b)
      implicit none
      real*8 a(ndim,ndim),b(ndim)
      integer neq,ndim
      integer nm1,l,lm1,k,j,jp1,n
      real*8 scr

      nm1 = neq-1
      do l=2,neq
         lm1=l-1
         do k=1,lm1
            b(l)=b(l)-a(l,k)*b(k)
         end do
      end do
      b(neq)=b(neq)/a(neq,neq)
      do l=1,nm1
         j=neq-l
         jp1=j+1
         do k=jp1,neq
            b(j)=b(j)-a(j,k  )*b(k)
         end do
         b(j)=b(j)/a(j,j)
      end do
      
      return
      end
