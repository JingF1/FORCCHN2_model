      subroutine  zeros3(m,n,k,b)!ÈıÎ¬0Õó
      integer m,n,k,i1,i2,i3
      real b(m,n,k)
      do i1=1,k
          do i2=1,m
             do i3=1,n
                 b(i2,i3,i1)=0.0
             end do
          end do
      end do
      end
      subroutine  zeros2(m,n,b)!¶şÎ¬0Õó
      integer m,n,i2,i3
      real b(m,n)
      do i2=1,m
          do i3=1,n
              b(i2,i3)=0.0
          end do
      end do   
      end
      subroutine  zeros1(m,b)!Ò»Î¬0Õó
      integer m,i2
      real b(m)
      do i2=1,m
          b(i2)=0.0     
      end do   
      end
       