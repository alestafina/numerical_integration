      PROGRAM main
      COMMON /dim/ n
      COMMON /gap/start,final
      REAL*8 Grid(1000000)
      REAL*8 res, Gaus5, Gaus2, trapz
      
      CALL input
      CALL makinGrid(Grid(1),n,start,final)

      res = trapz(Grid(1))
      PRINT*,'Trapz:', res
      res = Gaus2(Grid(1))
      PRINT*,'Gaus 2:', res
      res = Gaus5(Grid(1))
      PRINT*,'Gaus 5:', res
      
      END
     
      REAL*8 FUNCTION fun(x)
      REAL*8 x
      fun=10*x*sin(40*x)
      END
     
      SUBROUTINE input
      COMMON /dim/n
      COMMON /gap/start,final
      OPEN(1, FILE='input.txt', STATUS='old')
      READ(1,*)n,start,final
      CLOSE(1)
      END
     
      SUBROUTINE makinGrid(grid,m,beg,en)
      COMMON /dim/n
      COMMON /gap/start,final
      REAL*8 grid(*)
      h=(en - beg)/m
      DO i=0,m,1
        grid(i+1)=beg+i*h  
      END DO
      END
     
      REAL*8 FUNCTION trapz(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      REAL*8 grid(*)
      REAL*8 res, fun
      res = (fun(grid(1))+fun(grid(n+1)))/2D0
      DO i=2,n,1
          res = res+fun(grid(i)) 
      END DO
      res=res*(final-start)/n
      trapz=res
      END

      REAL*8 FUNCTION Gaus2(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      REAL*8 grid(*)
      REAL*8 q(2),x(2)
      REAL*8 sum1, sum2, fun, h
      h=(final-start)/n
      sum2=0.0
      
      x(1)=-1.0D0/SQRT(3.0D0)
      x(2)=1.0D0/SQRT(3.0D0)
      
      q(1)=1.0D0
      q(2)=1.0D0
      DO j=1,2
      sum1=0.0
      DO i=1,n
      sum1=sum1+h*fun(((grid(i))+grid(i+1))/2.0D0+x(j)*h/2.0D0)
      END DO
      sum2=sum2+q(j)*sum1
      END DO
      Gaus2=sum2/2.0D0
      END

      REAL*8 FUNCTION Gaus5(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      REAL*8 grid(*)
      REAL*8 q(5),x(5)
      REAL*8 sum1, sum2, fun, h
      h=(final-start)/n
      sum2=0.0
      
      x(1)=-1.0D0/3.0D0*(sqrt(5.0D0+2.0*sqrt(10D0/7.0D0)))
      x(2)=-1.0D0/3.0D0*(sqrt(5.0D0-2.0D0*sqrt(10D0/7.0D0)))
      x(3)=0.0D0
      x(4)=1.0D0/3.0D0*(sqrt(5.0D0-2.0D0*sqrt(10D0/7.0D0)))
      x(5)=1.0D0/3.0D0*(sqrt(5.0D0+2.0D0*sqrt(10D0/7.0D0)))
      
      q(1)=(322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
      q(2)=(322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
      q(3)=128.0D0/225.0D0
      q(4)=(322.0D0+13.0D0*sqrt(70.0D0))/900.0D0
      q(5)=(322.0D0-13.0D0*sqrt(70.0D0))/900.0D0
      DO j=1,5
      sum1=0.0
      DO i=1,n
      sum1=sum1+h*fun(((grid(i))+grid(i+1))/2.0D0+x(j)*h/2.0D0)
      END DO
      sum2=sum2+q(j)*sum1
      END DO
      Gaus5=sum2/2.0D0
      END