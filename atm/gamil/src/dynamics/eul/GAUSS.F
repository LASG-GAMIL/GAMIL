c
        subroutine gauss(a,b,n)
c
        implicit none
        integer n,i,k
        real*8 a(n-1,5),b(n+1),y
c
        do k=1,n-2
           i=k+1
           y=a(i,1)/a(k,2)
           a(i,1)=0.0
           a(i,2)=a(i,2)-a(k,3)*y
           a(i,4)=a(i,4)-a(k,4)*y
           a(i,5)=a(i,5)-a(k,5)*y
c
           y=b(k)/a(k,2)
           b(i)=b(i)-a(k,3)*y
           b(n)=b(n)-a(k,4)*y
           b(n+1)=b(n+1)-a(k,5)*y
        end do
        k=n-1
        i=n
        y=b(k)/a(k,2)
        b(n)=b(n)-a(k,4)*y
        b(n+1)=b(n+1)-a(k,5)*y
c
        b(n+1)=b(n+1)/b(n)
        a(n-1,5)=(a(n-1,5)-a(n-1,4)*b(n+1))/a(n-1,2)
        do k=n-2,1,-1
           y=a(k,4)*b(n+1)+a(k,3)*a(k+1,5)
           a(k,5)=(a(k,5)-y)/a(k,2)
        end do
c
        return
        end
