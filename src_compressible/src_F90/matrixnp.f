      subroutine decompnp(neq,ndim,a)
      implicit none
      integer neq,ndim
      real*8  a(ndim,ndim)
      real*8 eps,eps2,amult,apm,apn
      integer n,nm1,i,ip1,k,j,l
      data eps,eps2 /1.d-30,1.d-60/
      nm1=neq-1
      do 100 i=1,nm1
      ip1 = i+1
      do 120 l=ip1,neq
      amult=a(l,i)/a(i,i)
      a(l,i)=amult
      if(abs(amult).lt.eps2) go to 120
      do 130 k=ip1,neq
      a(l,k)=a(l,k)-amult*a(i,k)
  130 continue
  120 continue
  100 continue
      if(abs(a(neq,neq)).lt.eps) go to 1000
      return
 1000 continue
      write(6,*)a(neq,neq)
      stop"singular matrix in decompnp"
      end

      subroutine solvenp(neq,ndim,a,b)
      implicit none
      real*8 a(ndim,ndim),b(ndim)
      integer neq,ndim
      integer nm1,l,lm1,k,j,jp1,n
      real*8 scr
      nm1 = neq-1
      do 10 l=2,neq
      lm1=l-1
      do 20 k=1,lm1
      b(l)=b(l)-a(l,k)*b(k)
  20  continue
   10 continue
      b(neq)=b(neq)/a(neq,neq)
      do 30 l=1,nm1
      j=neq-l
      jp1=j+1
      do 40 k=jp1,neq
      b(j)=b(j)-a(j,k  )*b(k)
   40 continue
      b(j)=b(j)/a(j,j)
   30 continue
      return
      end
