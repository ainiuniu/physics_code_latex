      program main
      implicit none
      integer:: i,j,k,l,m,n,h
      real(kind=8) ::psi(100),p(10),w(10),sum(100),t(400)
      real(kind=8)::pi=3.141592653589793d0,eta=1.0d0/4.0d0*1.06d0**2
      real(kind=8)::mid,mid1,mid2,mid3 
      complex(kind=8)::sumt,thta,thta1
      real(kind=8)::sumt1(400),dsumt(400)
      
      do l=187,187
        t(l)=1.0d-3*l
        sumt=(0.d0,0.d0)

      do k=1,100
         psi(k)=k*2.0d0*pi/100.0d0
         sum(k)=0.0d0 
      
      mid=1.06d0/(4.0d0*pi*t(l))-psi(k)/(2.0d0*pi)
      mid1=-1.06d0/(4.0d0*pi*t(l))-psi(k)/(2.0d0*pi)
      !write(*,*) mid1,mid
      n=int(mid1)
      !write(*,*)n
      if(mid>=0.d0) then
      
      m=int(mid)
      
      do i=n,m
      
      
      mid2=eta-(2.d0*pi*t(l)*(i+psi(k)/(2.d0*pi)))**2.d0
      
      call GauLeg(0.0d0,mid2,p,w,10) 
      do j=1,10
      sum(k) =sum(k)+t(l)*(p(j)*(mid2-&
      &p(j)))**(1.d0/2.d0)/((eta)*2.d0*pi**2.d0)*w(j)
      end do
      end do
    
      else if(mid<0.d0.and.int(mid)/=n) then
      
      m=int(mid)-1
      do i=n,m
      mid2=eta-(2.0d0*pi*t(l)*(i+psi(k)/(2.0d0*pi)))**2
      call GauLeg(0.0d0,mid2,p,w,10)

      do j=1,10
      sum(k) =sum(k)+t(l)*(p(j)*(mid2-&
      &p(j)))**(1.d0/2.d0)/((eta)*2.0d0*pi**2)*w(j)
      end do

      end do

     
      else
      sum(k)=0.0d0
      end if
      
      open(10,file='1.txt')
      !open(20,file='2.txt')
      write(10,*)sum(k)
      !write(*,*) k,sum(k)
      !write(20,*)k
      !sumt1(l)=sum(k)
      end do
      
      do k=1,100
            thta=cmplx(0.D0,psi(k))
            sumt=sumt+exp(-thta)*sum(k)/100.0d0
      end do
      sumt1(l)=real(sumt)
      open(30,file='3.txt')
      write(30,*) sumt1(l)

      !stop
      end do
      
      !do h=1,399
       !     dsumt(h)=(sumt1(h+1)-sumt1(h))/1.d-3
        !    open(40,file='4.txt')
         !   write(40,*) dsumt(h)
      !end do
            
      end program
      
      subroutine GauLeg(x1,x2,x,w,n)

       implicit none  

       double precision, intent(in) :: x1, x2

       integer, intent(in) :: n

       double precision, dimension(n) :: x, w

       double precision, parameter :: eps=3.d-11

       double precision, parameter :: pi=3.141592653589793

       integer :: m, i, j, k

       double precision :: xm, xl, z, z1, p1, p2, p3, pp



       m=(n+1)/2 

       xm=0.5d0*(x1+x2)

       xl=0.5d0*(x2-x1)   

       z1=0.0d0   

         do i=1,m

           z=dcos(pi*(i-0.25d0)/(n+0.5d0)) 

           k=1

            do while(k .lt. 5)      

              p1=1.0d0

              p2=0.0d0

	            do j=1,n

                  p3=p2

                  p2=p1

                  p1=((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/dble(j)

                 end do

              pp=dble(n)*(z*p1-p2)/(z*z-1.0d0) 

              z1=z

              z=z1-p1/pp  

	           if(dabs(z-z1) .le. eps) k=10

             end do   

           x(i)=xm-xl*z

           x(n+1-i)=xm+xl*z

           w(i)=2.0d0*xl/((1.0d0-z*z)*pp*pp)

           w(n+1-i)=w(i)

         end do   

     end subroutine GauLeg


 
