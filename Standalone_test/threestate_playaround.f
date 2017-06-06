
      program threestate
 
c     Initialize everything.
c     -----------------------------------------------------------
c     To get the matrix multiplication function to work, need to 
c     redefine arrays as 43 x 1 matrices?

      real*8 N_overlap,  x_overlap, x_maxoverlap
      real*8 N_D1, N_D2, N_bound, N_bound_sum
      real*8 a1, a2, a3
      real*8 a_on_rate, a_off_rate, sum_of_rates, num_timesteps, t_p
      real*8 kd1, kd2, kd3, kd4, delta_D2
      real*8 kA_D2(41), kD2_A(41), bin_loc(41), bin_pops(41),
     &x_cb0(43,1), x_cb(43,1), moved_bins(41), x_cb_interp(41,1)
      real*8 M(43, 43)
      integer numsteps

      real*8:: N_a_active= 0.0, Ca = 0.0, k_mlcp = 3.29
      real*8:: k_force = 1.12*10E-4, T=288, k_b = 1.381E-23, hsl =
     & 1130.0
      real*8:: l_thin = 1012.0, l_thick = 815.0, l_bare = 80.0, hslc
      real*8:: time = 0.0, timestep = 0.001, t_act=0.1, hs_force = 0
      real*8:: k1 = 7.17, k3 = 13.1, k40 = 6.06, k41 =1.0,kD2_D1 = 100.0
      real*8:: k_cb = 0.002 , x_ps = 1.4, compliance = 0.5, cb_densit
     &y=6.9E16

      real*8:: a_on = 10E6
      real*8:: a_off =  10.0
      real*8:: k_thin_coop = 1.0
      integer p=1
      integer H=1

c     Will be specified at higher level in FE program
c      real*8:: hsl = 1130.0

c     Initialize initial conditions for states
      do 1 i = 1,41
           x_cb0(i,1) = 0
1     continue
      x_cb0(42,1) = 0.9
      x_cb0(43,1) = 0.1

      integer N_dstates = 2
      
c     Initialize Rate Constants kA_D2, kD2_A
c     --------------------------------------
      do 10 i = 1,41
            bin_loc(i) = ((i+1)*0.5-11)
c            kd1 = -1*k_cb*bin_loc(i)**2
c            kd2 = 2*k_b*T
c            kd3 = kd1/kd2
c            kd4 = exp(kd3)
c            kD2_A(i) = k3*kd4*(N_a_active-N_a_bound)

c       Attempting to do it all in one step. Seems to work
            kD2_A(i)=(k3)*exp((-1*k_cb*bin_loc(i)**2)/(2*k_b*T))
     &*(N_a_active-N_bound)

c           Attempting to set detachment as constant of same magnitude
c           as attachment
            kA_D2(i) = 7
c           kA_D2(i) = k40+k41*((bin_loc(i)-x_ps)**4)
10    continue
            kD1_D2 = (k1 + H*k_mlcp)*(1+k_force*hs_force)

c     Initialize DE Matrix
c     --------------------
      do 21 i = 1,43
        do 22 j = 1,43
           M(i,j) = 0.0
22      continue
21    continue
      do 20 i=1,41
            M(i,i) = -1*kA_D2(i)
            M(i,43) = kD2_A(i)
            M(43,i) = kA_D2(i)
20    continue
c     Don't know if there is a sum function. Will sum all D2->A rates
      do 30 i=1,41
            sum_of_rates = sum_of_rates + kD2_A(i)
30    continue
            
      M(42,42) = -1*kD1_D2
      M(43,43) = -1*(kD2_D1 + sum_of_rates)
      M(42,43) = kD2_D1
      M(43,42) = kD1_D2
      
      num_timesteps = int(.40/timestep)
      numsteps = num_timesteps
      

c     Open files to write data to for visualization. 12 is for actin
c     rates, 13 for myosin rate functions, and 14 for populations, force
c     ------------------------------------------------------------------
      open(11, file="calcium.txt")
      open(12, file="actin_rates.txt")
      open(13, file="kD1_D2.txt")
      open(14, file="kD2_D1.txt")
      open(16, file="kD2_A.txt")
      open(17, file="kA_D2.txt")
      open(18, file="populations_force.txt")
      open(19, file="mbound.txt")
      open(111, file="pop_distribution.txt")
      open(99, file="check.txt")
c     ---------------------
c     ---------------------
c     Begin Time Step Here
c     ---------------------
c     ---------------------

      do 1000 j = 1, numsteps      

c     reset this so it can be used again when updating DE matrix
      sum_of_rates = 0.0
      N_bound_sum = 0.0
      hs_force = 0.0

c     Calculate Calcium (Just a twitch for now)
c     -----------------------------------------
         
c      if (time.ge.t_act) then
c         t_p=t_act+0.01
c         if (time.ge.t_p) then
c           pCa=0.5*exp(-((time-t_p)*25)**2)
c         else
c           pCa=0
c         endif
c         Ca=(10**(-4.5))*sin(3.14*pCa)
c         if (time.eq.0.2) then
c            hslc = -50.0
c         else
c            hslc = 0.0
c         endif
c      else
c         Ca = 0.1
c      endif

       if (time.ge.t_act) then
           Ca = 10**(-4.5)
       else
           Ca = 0.0
       endif
       if (time.eq.0.2) then
           hslc = -50.0
       else
           hslc = 0.0
       endif   

c     Calculate N_overlap
c     -------------------
      x_overlap = l_thick + l_thin - hsl
      x_maxoverlap = l_thick - l_bare
      if (x_overlap.le.0.0) then
              N_overlap = 0.0
      elseif (x_overlap.ge.0.0 .AND. x_overlap.le. x_maxoverlap) then
              N_overlap = x_overlap/x_maxoverlap
      elseif (x_overlap.ge.x_maxoverlap) then
              N_overlap = 1.0
      endif
      
c     Calculate Number of Actin Active
c     --------------------------------
      a1 = (N_a_active/N_overlap)
      a2 = 1 + (k_thin_coop*a1)
      a3 = N_overlap - N_a_active
      a_on_rate = a_on*Ca*a3*a2

c       Trying to calculate a_on_rate in one sweep
c      a_on_rate = a_on*Ca*(N_overlap-N_a_active)*
c     &(1+(k_thin_coop*(N_a_active/N_overlap)))*0.1
      
      ao1 = N_a_active - N_bound
      ao2 = (N_overlap-N_a_active)/N_overlap
      ao3 = ao1*(1+ao2)
      a_off_rate = a_off*ao1*ao3
  
  
c      a_off_rate = a_off*(N_a_active - N_a_bound)*(1+k_thin_coop*(((N_ov
c     &erlap-N_a_active)/N_overlap)**p))
        
      N_a_active = N_a_active + (a_on_rate - a_off_rate)*timestep

c     Update DE Matrix (including rate constants)
c     -------------------------------------------
      do 40 i = 1,41
                kd1 = -1*k_cb*bin_loc(i)*bin_loc(i)*1e-18
                kd2 = (2*k_b*T)
                kd3 = kd1/kd2
                kd4 = exp(kd3)
                kD2_A(i) = k3*kd4*(N_a_active-N_bound)   

                sum_of_rates = sum_of_rates + kD2_A(i)
                kD1_D2 = (k1 + H*k_mlcp)*(1.0+k_force*hs_force)
                M(i,43) = kD2_A(i)
40    continue
                M(43,43) = -1*(kD2_D1 + sum_of_rates)
     

c     Runge Kutta for Updating number in each state
c     ---------------------------------------------
      call RK2(M, time, timestep, x_cb0, x_cb)
      N_D1 = x_cb(42,1)
      N_D2 = x_cb(43,1)
      
c     Assign everything from RK appropriately
c     ---------------------------------------
      
      do 41 i = 1,41
      N_bound_sum = N_bound_sum + x_cb(i,1)
41    continue 
      N_bound = N_bound_sum

c     Interpolate
c     -----------
      do 51 i = 1,41
         moved_bins(i) = bin_loc(i) - hslc*1e7*0.5
51    continue
         
      call pwl_interp_1d(41,bin_loc,x_cb,41,moved_bins,
     & x_cb_interp)

c     Move detached heads back into D2
c     --------------------------------
      do 63 i = 1,41
         delta_D2 = delta_D2 + (x_cb(i,1)-x_cb_interp(i,1))
63    continue
      x_cb(43,1) = x_cb(43,1) + delta_D2
      N_bound = N_bound - deltaD2

c     Calculate the half-sarc force
c     -----------------------------
c     Must sum force from each bin
      
      do 50 i = 1,41
        hs_force = hs_force+cb_density*k_cb*x_cb_interp(i,1)*(bin_loc(i)
     &+x_ps)

50    continue

      do 52 i = 1,41
      bin_pops(i) = x_cb_interp(i,1)
52    continue
      do 53 i = 1,41
      x_cb0(i,1) = x_cb_interp(i,1)
53    continue
      x_cb0(42,1) = x_cb(42,1)
      x_cb0(43,1) = x_cb(43,1)
c     Write out to text file for plotting, increment time. Formatting is
c     difficult because some rate functions are arrays whereas others
c     are a single value. Using separate text files
c     ----------------------------------------------------
      write(11, *) Ca, time
      write(12, *) a_on_rate, a_off_rate, N_a_active
      write(13, *) kD1_D2
      write(14, *) kD2_D1
      write(16, *) kD2_A
      write(17, *) kA_D2
      write(18, *) N_D1, N_D2, hs_force
      write(19, *) N_bound
      write(111, *) bin_pops
      time = time+timestep
      
1000  continue

c     Close files after simulation is finished
c     ----------------------------------------
      close(11)
      close(12)
      close(13)
      close(14)
      close(16)
      close(17)
      close(18)
      close(19)
      close(99)
      close(111)
      end program threestate
c--------------------------------------------------------------------------------
c--------------------------------------------------------------------------------

      subroutine RK(M, t, h, init_conds, yout)
         real*8 M(43,43), t, h, init_conds(43,1), yout(43,1)
         real*8 y_values(43,1), dyt(43,1), yt(43,1), dym(43,1),hh,h6,th
         real*8 dydt(43,1)
         y_values = init_conds
         hh = h*0.5
         h6 = h/6
         th = t+hh

c       Figure out how to code matrix multiplication
c       This will be hard coded keeping in mind the nature of M
         
c        M*y_values (Step 1)        
         do 15 i = 1,42
         dydt(i,1) = M(i,i)*y_values(i,1) + M(i,43)*y_values(43,1)
15       continue
         do 25 i = 1,43  
         dydt(43,1) = dydt(43,1) + M(43,i)*y_values(i,1)
25       continue
       
        
         

         do 35 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dydt(i,1)
35       continue
        
c        M*yt (Step 2)
         do 45 i = 1,42
         dyt(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
45       continue
         do 55 i = 1,43
         dyt(43,1) = dyt(43,1) + M(43,i)*yt(i,1)
55       continue

         do 65 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dyt(i,1)
65       continue
         
c        (Step 3)
         do 31 i = 1,42
         dym(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
31       continue
         do 32 i = 1,43
         dym(43,1) = dym(43,1) + M(43,i)*yt(i,1)
32       continue
         
         do 75 i = 1,43
         yt(i,1) = y_values(i,1) + h*dym(i,1)
75       continue

         do 33 i = 1,43
         dym(i,1) = dym(i,1)+dyt(i,1)
33       continue

c        (Step 4)
         do 36 i = 1,42
         dyt(i,1) = M(i,i)*yt(i,1) + M(i,43)*yt(43,1)
36       continue
         do 37 i = 1,43
         dyt(43,1) = dyt(43,1) + M(43,i)*yt(i,1)
37       continue

         do 85 i = 1,43
         yout(i,1) = y_values(i,1) + h6*(dydt(i,1)+dyt(i,1)+2*dym(i,1))
85       continue

 
         return 
         
      end

      subroutine pwl_interp_1d(nd,xd,yd,ni,xi,yi)
      integer nd,ni,i,k
      real*8 xd(nd),yd(nd,1),xi(ni),yi(ni,1)
      real*8 ff
      do 200 i=1,ni
        yi(i,1)=0.0
200    continue 
      if (nd.eq.1) then
        do 201 i=1,ni
          yi(i,1)=yd(1,1)
201      continue      
        return
      endif
      do 202 i=1,ni
        if (xi(i).lt.xd(1).or.xi(i).gt.xd(nd)) then
           yi(i,1)=0.0
        else
          do 203 k=2,nd
            if (xd(k-1).le.xi(i).and.xi(i).le.xd(k)) then
              ff=(xi(i)-xd(k-1))/(xd(k)-xd(k-1))
              yi(i,1)=(1.0-ff)*yd(k-1,1)+ff*yd(k,1)
            endif
203        continue
        endif
202    continue 
      return
      end 


c     Subroutine RK2 which will utlize matrix_mult
c     --------------------------------------------
      subroutine RK2(M, t, h, init_conds, yout)
         real*8 M(43,43), t, h, init_conds(43,1), yout(43,1)
         real*8 y_values(43,1), dyt(43,1), yt(43,1), dym(43,1), hh,h6,th
         real*8 dydt(43,1)
         
         y_values = init_conds
         hh = h*0.5
         h6 = h/6
         th = t+hh

c       Figure out how to code matrix multiplication
c       This will be hard coded keeping in mind the nature of M
         
c        M*y_values (Step 1)    

         call matrix_mult(M,43,43,y_values,43,1,dydt)


         do 35 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dydt(i,1)
35       continue
        

c        M*yt (Step 2)

         call matrix_mult(M,43,43,yt,43,1,dyt)

         do 65 i = 1,43
         yt(i,1) = y_values(i,1) + hh*dyt(i,1)
65       continue

c        (Step 3)

         call matrix_mult(M,43,43,yt,43,1,dym)
         
         do 75 i = 1,43
         yt(i,1) = y_values(i,1) + h*dym(i,1)
75       continue

         do 33 i = 1,43
         dym(i,1) = dym(i,1)+dyt(i,1)
33       continue

c        (Step 4)

         call matrix_mult(M,43,43,yt,43,1,dyt)

         do 85 i = 1,43
         yout(i,1) = y_values(i,1) + h6*(dydt(i,1)+dyt(i,1)+2*dym(i,1))
85       continue
         write(99, *) M
         return 
         
      end

      

c     Trying to implement my matrix algorithm subroutine
      subroutine matrix_mult(A,m,n,B,p,q,C)
      integer i,j,k,m,n,p,q,counter1,counter2
      real*8 A(m,n), B(p,q), C(m,q)

c      i = 1
c      j = 1
c      k = 1
      

c     Initialize resultant matrix to zeros, else it populates weird
      do 400 counter1 = 1,m
        do 401 counter2 = 1,q
          C(counter1,counter2) = 0.0
401     continue
400   continue
                
      if (n.eq.p) then
        continue
      else
        print *, "Error"
      endif

      do 300 i = 1,m

        do 301 j = 1,q

          do 302 k = 1,n

            C(i,j) = C(i,j) + A(i,q)*B(q,j)

302       continue

301     continue

300   continue
      
      return
      end
