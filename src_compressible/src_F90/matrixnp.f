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

! zero-based version for C++ conversion      
c$$$      subroutine decompnp(neq,ndim,a)
c$$$      implicit none
c$$$      integer neq,ndim
c$$$      real*8  a(0:ndim-1,0:ndim-1)
c$$$      real*8 eps,eps2,amult
c$$$      integer n,nm1,i,ip1,k,j,l
c$$$
c$$$      eps = 1.d-30
c$$$      eps2 = 1.d-60
c$$$      
c$$$      nm1=neq-1
c$$$      do i=0,nm1-1
c$$$         ip1 = i+2
c$$$         do l=ip1-1,neq-1
c$$$            amult=a(l,i)/a(i,i)
c$$$            a(l,i)=amult
c$$$            if(abs(amult).lt.eps2) cycle
c$$$            do k=ip1-1,neq-1
c$$$               a(l,k)=a(l,k)-amult*a(i,k)
c$$$            end do
c$$$         end do
c$$$      end do
c$$$      if(abs(a(neq-1,neq-1)).lt.eps) then
c$$$         write(6,*)a(neq,neq)
c$$$         stop"singular matrix in decompnp"
c$$$      end if
c$$$      
c$$$      return
c$$$      end

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


!     zero-based version for C++ conversion
c$$$      subroutine solvenp(neq,ndim,a,b)
c$$$      implicit none
c$$$      real*8 a(0:ndim-1,0:ndim-1),b(0:ndim-1)
c$$$      integer neq,ndim
c$$$      integer nm1,l,lm1,k,j,jp1,n
c$$$      real*8 scr
c$$$
c$$$      nm1 = neq-1
c$$$      do l=1,neq-1
c$$$         lm1=l
c$$$         do k=0,lm1-1
c$$$            b(l)=b(l)-a(l,k)*b(k)
c$$$         end do
c$$$      end do
c$$$      b(neq-1)=b(neq-1)/a(neq-1,neq-1)
c$$$      do l=0,nm1-1
c$$$         j=neq-(l+1)
c$$$         jp1=j+1
c$$$         do k=jp1-1,neq-1
c$$$            b(j-1)=b(j-1)-a(j-1,k)*b(k)
c$$$         end do
c$$$         b(j-1)=b(j-1)/a(j-1,j-1)
c$$$      end do
c$$$      
c$$$      return
c$$$      end
