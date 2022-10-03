      PROGRAM main
      COMMON /dim/ n
      COMMON /gap/start,final
      DIMENSION Grid(1000000)
      REAL*8 res
      CALL input
      CALL makinGrid(Grid(1),n,start,final)
  
      res = trapz(Grid(1))
      PRINT*,'Trapz:', res
      res = Gaus2(Grid(1))
      PRINT*,'Gaus 2:', res
      res = Gaus5(Grid(1))
      PRINT*,'Gaus 5:', res
      
      END
     
      FUNCTION fun(x)
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
      DIMENSION grid(*)
      h=(en - beg)/m
      DO i=0,m,1
        grid(i+1)=beg+i*h  
      END DO
      END
     
      FUNCTION trapz(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      DIMENSION grid(*)

      res = (fun(grid(1))+fun(grid(n+1)))/2
      DO i=2,n,1
          res = res+fun(grid(i)) 
      END DO
      res=res*(final-start)/n
      trapz=res
      END
  
      FUNCTION Gaus2(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      DIMENSION grid(*), q(2),x(2)
        
      h=(final-start)/n
      sum2=0.0
      
      x(1)=-1.0/SQRT(3.0)
      x(2)=1.0/SQRT(3.0)
      
      q(1)=1.0
      q(2)=1.0

      DO j=1,2
      sum1=0.0
      DO i=1,n
      sum1=sum1+h*fun(((grid(i))+grid(i+1))/2.0+x(j)*h/2.0)
      END DO
      sum2=sum2+q(j)*sum1
      END DO
      Gaus2=sum2/2.0
      END
  
      FUNCTION Gaus5(grid)
      COMMON /dim/ n
      COMMON /gap/start,final
      DIMENSION grid(*), q(5),x(5)
      h=(final-start)/n
      sum2=0.0
      
      x(1)=-1.0/3.0*(sqrt(5.0+2.0*sqrt(10/7.0)))
      x(2)=-1.0/3.0*(sqrt(5.0-2.0*sqrt(10/7.0)))
      x(3)=0.0
      x(4)=1.0/3.0*(sqrt(5.0-2.0*sqrt(10/7.0)))
      x(5)=1.0/3.0*(sqrt(5.0+2.0*sqrt(10/7.0)))
      
      q(1)=(322.0-13.0*sqrt(70.0))/900.0
      q(2)=(322.0+13.0*sqrt(70.0))/900.0
      q(3)=128.0/225.0
      q(4)=(322.0+13.0*sqrt(70.0))/900.0
      q(5)=(322.0-13.0*sqrt(70.0))/900.0
      DO j=1,5
      sum1=0.0
      DO i=1,n
      sum1=sum1+h*fun(((grid(i))+grid(i+1))/2.0+x(j)*h/2.0)
      END DO
      sum2=sum2+q(j)*sum1
      END DO
      Gaus5=sum2/2.0
      END