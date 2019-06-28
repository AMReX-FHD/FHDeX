      subroutine decomp(neq,ndim,a,ip)
      implicit none
      integer neq,ndim
      real*8  a(ndim,ndim)
      integer ip(ndim)
      real*8 eps,eps2,amult,apm,apn
      integer n,nm1,i,ip1,k,j,l
      data eps,eps2 /1.d-30,1.d-60/
      nm1=neq-1
      do 10 n=1,neq
  10  ip(n)=n
      do 100 i=1,nm1
      ip1=i+1  
       k=i
      apm    =abs(a(ip(i),i))
      do 110 j=ip1,neq
      apn=abs(a(ip(j),i))
      if(apm.ge.apn) go to 110
      apm=apn  
       k=j
  110 continue
      j=ip(k)
      ip(k)=ip(i)
      ip(i)=j
  115 if(apm.lt.eps) go to 1000
      do 120 l=ip1,neq
      n=ip(l)
      amult=a(n,i)/a(j,i)
      a(n,i)=amult
      if(abs(amult).lt.eps2) go to 120
      do 130 k=ip1,neq
      a(n,k)=a(n,k)-amult*a(j,k)
  130 continue
  120 continue
  100 continue
      if(abs(a(ip(neq),neq)).lt.eps) go to 1000
      return
 1000 continue
      write(6,*)a(ip(neq),neq)
      stop"singular matrix in decomp"
      end

      subroutine solve(neq,ndim,a,b,ip)
      implicit none
      real*8 a(ndim,ndim),b(ndim)
      integer ip(ndim)
      integer neq,ndim
      integer nm1,l,lm1,k,j,jp1,n
      real*8 scr
      nm1 = neq-1
      do 10 l=2,neq
      n=ip(l)
      lm1=l-1
      do 20 k=1,lm1
      b(n)=b(n)-a(n,k)*b(ip(k))
  20  continue
   10 continue
      b(ip(neq))=b(ip(neq))/a(ip(neq),neq)
      do 30 l=1,nm1
      j=neq-l
      jp1=j+1
      n=ip(j)
      do 40 k=jp1,neq
      b(n)=b(n)-a(n,k  )*b(ip(k))
   40 continue
      b(n)=b(n)/a(n,j)
   30 continue
      do 50 n=1,neq
   60 if(ip(n) .eq. n) go to 50
      j=ip(n)
      scr=b(j)
      ip(n)=ip(j)
      b(j)=b(ip(j))
      b(ip(j))=scr
      ip(j)=j
      go to 60
   50 continue
      return
      end
