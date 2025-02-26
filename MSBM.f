         subroutine vumat(
C Read only -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     3  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only -
     5  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C

C
C
C All arrays dimensioned by (*) are not used in this algorithm
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock) 


      character*80 cmname
C     onum: number of orientations to consider
C     dnum: number of initial lengths to consider
C     Fil_len: array of filament lengths
C     Fil_Freq: array of filament length frequency
C     lamnum: number of stretch ticks in one direction
      integer onum, modnum, dnum, lamnum, istage
      parameter (onum = 3, dnum = 12, modnum = 5000, lamnum = 5000)
      real*8 pi
      parameter (pi = 3.1415926)
      
      real*8 Fil_len(dnum), Fil_Freq(dnum)
C     tau_a: mean stress resultant
C     mu_g: shear modulus
C     area: deformation area, area_inv: area invariant
C     zrho: initial network density
      real*8 l1, l2 ! Principal stretches
      real*8 tau1, tau2 ! Principal stresses
	  real*8 tau_a, tau_s, tau_ai ! Stress resultants and membrane tension
      real*8 mu_g, area, theta
      real*8 U(4), B(4), B_trace, F
      real*8 UNFOLD_MODEL
      real*8 lammin, lammax, dlam
      real*8 dimernum, A0, r0, dx, maxtau
      integer index1, index2, i, n, m, k, l, ii, jj
      real*8 x(modnum),y(modnum), x2(modnum),y2(modnum)
      real*8 tau1_dat(lamnum,lamnum),tau2_dat(lamnum,lamnum)
      real*8 thick, Kbilayer, zrho
      real*8 Fmax, Fmin
      
      save istage
      save tau1_dat,tau2_dat
      save dlam, lammin, lammax
      
      
      Fmax = 2500.0 ! Maxiumum force
      Fmin = -2.654 ! Minimum force
      dimernum = 1778.7*2.0 ! Number of dimers in RVE
      A0 = 1.4*1.4 ! Area of RVE (um^2)
      zrho = dimernum/A0 ! Network density (1/um^2)
      r0 = 0.0460 ! Initial dimer length (um)

      Kbilayer = props(1) ! Area modulus of lipid bilayer
      thick = props(2) ! Thickness of NL layer (um)

      ! Pre-computing a discrete force-extension dataset for the lamin dimer
      if(istage .eq. 0) then
      do 10 i = 1, modnum
        y(i) = (Fmin + (Fmax)*(i - 1)/(modnum*1.0))
        x(i) = UNFOLD_MODEL(y(i))
   10 continue
  
      ! Creating a 2D lookup table of principal stresses based on principal stretches
      lammax = x(modnum)*1.0
      lammin = x(1)*1.0
      
      dlam = (lammax - lammin)/(1.0*lamnum)
      do 20 m = 1, lamnum
      do 30 n = 1, lamnum
      
      l1 = lammin + dlam*1.0*m
      l2 = lammin + dlam*1.0*n
      if (abs(l1 - l2) .lt. 1.0E-6) then
      l1 = l1*(1 + 1.0E-6)
      l2 = l2/(1 + 1.0E-6)
      end if

      tau1 = 0.0
      tau2 = 0.0
          
      do 40 k = 1, onum  ! distribution of orientation
        theta = k*pi/(1.0*onum)
        r_r0 = (l2*sin(theta))**2
        r_r0 = r_r0 + (l1*cos(theta))**2
        r_r0 = sqrt(r_r0)                                           

        do 50 l = 1, (modnum - 1)
          if (r_r0 .ge. x(l) .and. r_r0 .le. x(l+1)) then
          F = (r_r0 - x(l + 1))/(x(l) - x(l + 1))*y(l)
          F = F + (r_r0 - x(l))/(x(l + 1) - x(l))*y(l + 1)
          exit
          endif
   50   continue
   
        tau1 = tau1 + F*(1/r_r0)*cos(theta)**2
        tau2 = tau2 + F*(1/r_r0)*sin(theta)**2

   40 continue
      tau1 = tau1*l1*r0*zrho/(1.0*onum)/l2/thick
      tau2 = tau2*l2*r0*zrho/(1.0*onum)/l1/thick

      tau_a = (tau1 + tau2)/2.0
      tau_s = (tau1 - tau2)/2.0
      mu_g = 2.0*tau_s*(l1*l2)**2/(l1**2 - l2**2)

      tau1_dat(m,n) = tau1
      tau2_dat(m,n) = tau2

   30 continue
   20 continue
      istage = 1
      endif
      
      do 60 i = 1, nblock
        tau1 = 0.0
        tau2 = 0.0 
        l1 = stretchNew(i,1)
        l2 = stretchNew(i,2)

        if (abs(l1 - l2) .lt. 1.0E-5) then
        l1 = l1*(1.0 + 1.0E-5)
        l2 = l2/(1.0 + 1.0E-5)
        end if
        ii = ceiling((l1 - lammin)/dlam)
        jj = ceiling((l2 - lammin)/dlam)
        
	
        ! Bilinear interpolation
        temp1 = ((dlam*(ii + 1) + lammin) - l1)/dlam*tau1_dat(ii,jj)  
        temp1 = temp1 + (l1 - (dlam*ii + lammin))/dlam*tau1_dat(ii+1,jj)
        temp2 = ((dlam*(ii + 1) + lammin) - l1)/dlam*tau1_dat(ii,jj+1)
        temp2 = temp2 + (l1 - (dlam*ii + lammin))/dlam*tau1_dat(ii+1,jj+1)
          
        tau1 = temp1*((dlam*(jj + 1) + lammin) - l2)/dlam
        tau1 = tau1 + temp2*(l2 - (dlam*jj + lammin))/dlam
        temp1 = ((dlam*(ii + 1) + lammin) - l1)/dlam*tau2_dat(ii,jj)
        temp1 = temp1 + (l1 - (dlam*ii + lammin))/dlam*tau2_dat(ii+1,jj)
        temp2 = ((dlam*(ii + 1) + lammin) - l1)/dlam*tau2_dat(ii,jj+1)
        temp2 = temp2 + (l1 - (dlam*ii + lammin))/dlam*tau2_dat(ii+1,jj+1)
          
        tau2 = temp1*((dlam*(jj + 1) + lammin) - l2)/dlam
        tau2 = tau2 + temp2*(l2 - (dlam*jj + lammin))/dlam

        tau_a = (tau1 + tau2)/2.0   !!- c0/(l1*l2)
        tau_s = (tau1 - tau2)/2.0
        mu_g = 2.0*tau_s*(l1*l2)**2/(l1**2 - l2**2)

        area_inv = l1*l2 - 1.0
        tau_ai = Kbilayer/thick*area_inv
        stateNew(i,1) = tau1
        stateNew(i,2) = tau2
        stateNew(i,3) = mu_g
        stateNew(i,4) = l1
        stateNew(i,5) = l2
        stressNew(i,1)=tau_ai + tau1
        stressNew(i,2)=tau_ai + tau2
        
   60 continue
       end

      real*8 function UNFOLD_MODEL(F)

        include 'vaba_param.inc'
        real*8 F_half, DDx_star, Pu, Pc, kB, T
        real*8 r0, Lc, Lu, fc, tu, phi_u
        real*8 rc_r0, ru_r0, pi, mu0c, mu0u, F
        parameter (pi = 3.1415926)
        r0 = 0.0468 ! Dimer length (um)
        
        kB = 1.38065e-5 ! Boltzmann's constant (pN.um/K)
        T = 300.0 ! Temperature (K)
        
        F_half = 280.0/4.0 ! Force to unfold half on alpha-domains (pN)
        DDx_star = 0.0003 ! Activation length difference (um)
        
        Lc = 0.0495 ! Contour length of folded region (um)
        Pc = 0.1596 ! Persistence length of folded region (um)
        mu0c = 50000.0 ! Extension constant of coiled dimer (pN/um)
        
        Lu = 0.12925 ! Contour length of unfolded region (um)
        Pu = 0.0001154 ! Persistence length of unfolded region (um)
        mu0u = 90000.0 ! Extension constant of uncoiled dimer (pN/um)
        
        
        tu = F*Pu/kB/T
        fc = pi**2*kB*T*Pc/Lc**2

        rc_r0 = (1 + F/Lc/mu0c)*(1 - sqrt((kB*T)/(pi*Pc*(F + fc))))*Lc/r0

        ru_r0 = (4.0/3.0 - 4.0/(3.0*sqrt(tu + 1.0)) -
     1          10.0*exp((900.0/tu)**(0.25))/
     2          (sqrt(tu)*(exp((900.0/tu)**(0.25)) - 1.0)**2) +
     3          (tu)**1.62/(3.55 + 3.8*tu**2.2) + F/mu0u)*Lu/r0

        phi_u = exp((F - F_half)*DDx_star/(kB*T))
        phi_u = phi_u/(1.0+ phi_u)
        
        
        if (F .gt. 0.0 .and. F .lt. 5000.0) then
        
        UNFOLD_MODEL = rc_r0*(1 - phi_u) + ru_r0*phi_u
        
        elseif (F .ge. 5000.0) then
        
        UNFOLD_MODEL = ru_r0
        
        elseif (F .le. 0.0) then
        UNFOLD_MODEL = rc_r0

        endif
        
        return
      end
