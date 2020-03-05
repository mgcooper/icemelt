c mgc run2, added ro_snow to outputs, changed ice_avail to dy_p, 
c and compute column-average density to update spect_ext_coeffs

c Modified by mgc for calendar year 2016 Greenland run

c icemelt.f

c This FORTRAN code simulates the temperature within snow and ice
c   in response to the surface energy balance and penetration of
c   solar radiation into the snow and ice.
c
c The general model equations are described in the paper:
c   Below-surface ice melt on the coastal Antarctic ice sheet, by
c   Glen E. Liston and 4 others, Journal of Glaciology, 1999,
c   Vol. 45, No. 150, pages 273-285.
c
c The author of this code is:
c   Dr. Glen E. Liston
c   InterWorks Consulting
c   15048 NCR 25E
c   Loveland, Colorado 80538
c
c   Voice: (970) 491-8220
c   Internet: liston@iceberg.atmos.colostate.edu
c
c The model should run with any FORTRAN 77 compiler.

c mgc moved this down further, iyear is used to determine isleap
c maxiter is the number of time steps in the model run.
c      parameter (maxiter=8784) ! mgc leap year
c      parameter (maxiter=8760)
c     parameter (maxiter=365)
c     parameter (maxiter=1)

c mgc addition
      integer iyear,imonth,iday  ! model year, month, and day
c mgc addition
	  
c nz=JJ equals the number of grid cells in the z direction.  The
c   reason this is like this is because my heat equation solver
c   calls the z-dir(k) the y-dir(j).
      parameter(nz=500)
      parameter(JJ=nz)
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
	  real ro_snow_z(JJ+2)
	  
c mgc add NOTE: before actually draining water and reducing density,
c I will try computing N* as though density had decreased
c	  real :: drainthresh	  
c	  real Nstar(JJ+2)
c mgc this should implicitly create Nstar(JJ+2)	  	  
c	  real :: ro_star, ro_star_old

c PROVIDE SOME OF THE RUN DEPENDENT INFORMATION.

c Number of times to loop through the year, to ensure convergence
c  of deep ice temperatures.
      max_annual_loops = 1
c     max_annual_loops = 3 ! mgc changed to 3

c Julian day of the model start.  Usually I start the melt runs
c   in the middle of winter.
c     J_day_start = 182
      J_day_start = 1 ! mgc start Jan 1 for Greenland

c Height of wind and temperature observations.
      z_windobs = 2.0

c Surface elevation.
      topo = 1271.0 ! mgc mean elevation at Kan-M weather station

c Surface roughness length.
      z_0 = 0.001

c Snow-water-equivalent depth.  Any non-zero value makes the model
c   simulate snow-ice conditions.
c     swe_depth = 10.0
      swe_depth = 0.0 ! mgc 

c Model time step.
c     dt = 86400.
      dt = 3600.

c Latitude of center of domain, in decimal degrees.
c     xlat = -71.4
      xlat = 67.04 ! mgc mean latitude at Kan-M during 2016

c Topographic slopes.
      slope_azimuth = 0.0
      terrain_slope = 0.0

c dz in ice model.
      deltaz = 0.03

c Initial conditions for ice profile.  This should be equal to
c   something like the mean annual air temperature.
c	  temp_ice_init_C = -17.0
c     temp_ice_init_C = 0.0
c	  temp_ice_init_C = -7.0 ! mgc mean annual ice temp at Kan-M
      temp_ice_init_C = -3.0 ! mgc mean annual Tair at Kan-M is -2.65
	  
c Fractional cloud cover and transmissivity.  These values can be
c   used to adjust the incoming solar radiation.
c     cloud_frac = 0.50
c     transmiss = 0.70
      cloud_frac = 0.0 ! mgc set to 0
      transmiss = 0.70

c Identify whether this is run includes a non-zero conduction term
c   (icond_flag = 0 = no conduction, 1 = conduction).
c   Include conduction for Antarctic simulations.
      icond_flag = 1

c Open the meteorological data file.
c     open (31,file='test_met_daily.dat',form='formatted')
      open (31,file='test_met_1hrly.dat',form='formatted')

c mgc add to determine if iyear is leap
	  read (31,*) iyear_tmp

      if (MOD(iyear_tmp,4) .eq. 0) then
        maxiter = 8784
      else
        maxiter = 8760
      endif
	  
	  rewind 31
c mgc end
	  
c Output files.
      open (18,file='enbal.dat')
      open (19,file='ice1.dat')
      open (20,file='ice2.dat')

c Density, grain radius, and albedo.  See the extcoefs subroutine
c   about how to determine the snow_grain_radius index.
      ro_snow = 800.0 ! mgc leaving this as 800
	  
      n_snowgrain_radius = 32 ! mgc index 32 = 2.5 mm
      albedo = 0.65

c Get the general constants to be used.
        CALL CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water,
     &    xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow)

c Supply the initial configurations for the ice-free lake model.
        CALL ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ,
     &    Tf,water_frac,gamma,xk_snow,water_depth_old,temp_ice_init_C,
     &    deltaz,ro_snow_z)

c Calculate the solar radiation within the snow/ice matrix.
c mgc add ro_star (removed for now)	 
        CALL EXTCOEFS(nz,deltaz,albedo,ro_snow,up,down,
     &    n_snowgrain_radius,total_solar)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc add to highlight the main loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc this section loops through max_annual_loops and calls ENBALANCE
c and ICE_ENERGY for 1:maxiter and writes the output

        do kkk=1,max_annual_loops
c          print *,'Annual Loop Number =',kkk

          do iter=1,maxiter
c            print *,'WORKING ON DAY =',iter ! mgc suppressed all print statements

c Read the atmospheric forcing data for this time step.
cc          read (31,rec=iter) Tair,rh,windspd
c           read (31,*) Tair,rh,windspd
c mgc want yyyy mmm dd written out 
c            read (31,*) dum1,dum2,dum3,ihour,Tair,rh,windspd,dum5,dum6
            read (31,*) iyear,imonth,iday,ihour,Tair,rh,windspd,dum5,dum6
			Tair = Tair + Tf

c mgc temporary fix 
			if (rh.lt.60.0) rh = 60.0
			if (windspd.lt.0.1) windspd = 0.1
c mgc temporary fix

            CALL ENBALANCE(Tair,windspd,rh,
     &        Tsfc,Qsi,Qli,Qle,Qh,
     &        Qe,Qc,Qm,balance,Qf,
     &        swe_depth,topo,z_windobs,dt,gamma,
     &        T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0,
     &        J_day_start,iter,xlat,slope_az,terrain_slope,
     &        transmiss,ihour)

            CALL ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf)

c mgc NOTE this is roughly where Hoffman updates density based on
c drained water but it may be embedded in the k=1,JJ loop
	 
c Save the energy balance data.
c           write (18,*) Tair-Tf,Tsfc-Tf,rh,windspd,
c    &        Qsi,Qli,Qle,Qh,Qe,Qc,Qm,balance,albedo
c            write (18,88) iter,Tair-Tf,Tsfc-Tf,rh,windspd,
c     &        Qsi,Qli,Qle,Qh,Qe,Qc,Qm,balance,albedo

c mgc added write out year month day and hour
            write (18,88) iter,iyear,imonth,iday,ihour,Tair-Tf,
     &        Tsfc-Tf,rh,windspd,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,balance,albedo
c mgc end

c mgc use this to debug 
c            print *, iter,Tair-Tf,Tsfc-Tf,rh,windspd,
c     &        Qsi,Qli,Qle,Qh,Qe,Qc,Qm,balance,albedo
c mgc end            
			
c mgc add ro_snow			
c Save the ice data.
            write (19,*) Tsfc-Tf,Qsi,Qsip,water_depth,water_flux,ro_snow

c mgc commenting this out and also k below			
c            write (20,*) 'iteration = ',iter

            do k=1,JJ			
c              write (20,*) k,T_old(k)-Tf,gamma(k),xmelt(k),water_frac(k)
              write (20,*) T_old(k)-Tf,gamma(k),xmelt(k),water_frac(k),
     &          ro_snow_z(k)
c mgc add
c			  print *, iter,k,water_frac(k)
            enddo

          enddo
        enddo

  88  format (i6,i5,i4,i4,i6,f8.1,f10.4,f8.1,f10.4,7f9.3,f14.9,f8.2)

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ICE SECTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc this section computes energy to melt/freeze ice at each depth j
c and the total column water depth (and flux) by solving the ice temp
c equation (ICEHEAT) and then computing meltwater produced by extra 
c available energy or freezing by deficit (ICEMF) and then computes 
c a new temp profile (GETNEWT)

      SUBROUTINE ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &  water_depth_old,water_flux,xLf)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xmelt(JJ+2)
      real T_tmp(JJ+2)

c Save a copy of the temperature at the previous time step.  This
c   will be used to compute the ice temperature profile for the
c   case where liquid water is present, after a computation is 
c   made to account for the amount of ice melted or water frozen.
      do j=1,JJ
        T_tmp(j) = T_old(j)
      enddo

c Solve the ice temperature equation.
      CALL ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water)

c Correct for ice temperatures above freezing, and compute the
c   meltwater produced by the extra available energy.  Also deal
c   with the case of refreezing water.
      CALL ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice,
     &  ro_snow,water_frac,flag)

c If water is present, recompute the temperature profile.  Making
c   sure the water areas are at Tf.
      if (flag.eq.1.0) then
        CALL GETNEWT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &    dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &    ro_snow,Cp_snow,T_tmp,Tf,xk_water,water_frac,
     &    up,down,total_solar,xLs,Rv,ro_water)
      endif

c Compute a total-column water depth.
      water_depth = 0.0
      do j=1,JJ
        water_depth = water_depth + dy_p(j) * water_frac(j)
      enddo

c Compute the water flux by looking at the difference between the
c   current water depth and that at the previous time step.
      water_flux = water_depth - water_depth_old
      water_depth_old = water_depth
	  
c mgc update ro_snow, need to change 800 to ro_snow_init or something
	  ro_snow = 800 * (1 - water_depth/15.0) ! total column depth is 15.	  

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc calls main constants and computes thermal cond. (xk_snow)

      SUBROUTINE CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water,
     &  xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow)

      xLs = 2.500e6
      xLf = 3.34e5
      Tf = 273.16
      ro_water = 1000.0
      Cp_water = 4180.0
      xk_water = 0.552
      xk_ice = 2.10
      ro_ice = 917.0
      Cp_snow = 2106.0
      Rv = 461.0

c Compute the thermal conductivity from the snow density.
      if (ro_snow.lt.156.0) then
        xk_snow = 0.023 + 0.234 * (ro_snow/1000.0)
      else
        xk_snow = 0.138 - 1.01 * (ro_snow/1000.0) + 3.233 *
     &    (ro_snow/1000.0)**2
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc computes water_frac(j) NOTE: I think this is where I would 
c re-define ice density based on the drained water, or it may be
c better to do it at the top, after the timestep completes and 
c before the next iteration

      SUBROUTINE ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice,
     &  ro_snow,water_frac,flag)

      real dy_p(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real freeze(JJ+2)
      real ice_avail(JJ+2)
      real water_frac(JJ+2)

c Compute the maximum ice available to be melted.
c mgc note ro_ice = 917, ro_snow is fixed at 800
      do j=1,JJ
        ice_avail(j) = ro_snow / ro_ice * dy_p(j)
      enddo

      flag = 0.0
      extramelt = 0.0
      do j=1,JJ
        if (T_old(j).ge.Tf  .or.  extramelt.gt.0.0) then

c Compute the amount of water produced by the energy represented
c   by the temperature above freezing.
          xmelt(j) = Cp_snow * dy_p(j) * (T_old(j) - Tf) / xLf
          totalmelt = xmelt(j) + extramelt
		  
c Add the new water to the old water.
c mgc this might be a mistake, water_frac should = totalmelt/dy_p		  
c          water_frac(j) = water_frac(j) + totalmelt / ice_avail(j)
          water_frac(j) = water_frac(j) + totalmelt / dy_p(j)		  

c Assume that energy producing a water fraction above 1.0 goes 
c   into melting the ice below (as 'extramelt').
          extramelt = max(0.0,water_frac(j) - 1.0) * ice_avail(j)
          water_frac(j) = min(1.0,water_frac(j))

c mgc curious how much 'extramelt' is produced
c		  if (extramelt.gt.0.0) then
c		  	print *,iter,extramelt
c		endif
c Because ice is melting, the temperature must be Tf.
          T_old(j) = Tf
          flag = 1.0

        else

c Three cases are possible, 1) there is no water to freeze, so do
c   nothing, 2) all of the water freezes and the extra energy drops
c   the  energy below Tf, or 3) only some of the water freezes and
c   the temperature remains at Tf.

          if (water_frac(j).gt.0.0) then

c Compute the amount of water frozen by the energy represented
c   by the temperature below freezing.
          freeze(j) = Cp_snow * dy_p(j) * (Tf - T_old(j)) / xLf

            if (freeze(j).le.water_frac(j)*ice_avail(j)) then
c Case 3.
              water_frac(j) = water_frac(j) - freeze(j) / ice_avail(j)
              T_old(j) = Tf
              flag = 1.0

            else
c Case 2.
              freeze(j) = water_frac(j) * ice_avail(j)
              water_frac(j) = 0.0
              T_old(j) = Tf - freeze(j) * xLf / (Cp_snow * dy_p(j))
              flag = 1.0

            endif
          endif

        endif
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETNEWT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,T_tmp,Tf,xk_water,water_frac,
     &  up,down,total_solar,xLs,Rv,ro_water)

      real gamma(JJ+2)
      real g_b_ns(JJ+2)
      real f_n(JJ+2)
      real aN(JJ+2)
      real aP0(JJ+2)
      real aS(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)
      real b_vector(JJ+2)
      real T_old(JJ+2)
      real T_tmp(JJ+2)
      real Sc(JJ+2)
      real Sp(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xynet(JJ+2)

c Compute the solar radiation penetrating the lake surface.
      CALL SOLARPEN(Qsip,Qsi,albedo)

      xTsfc = Tsfc

c Compute gamma and the general equation coefficients.
      CAll GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water)
      CALL GAMMA1(g_b_ns,gamma,f_n,JJ)
      CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ,
     &  ro_snow,Cp_snow)

c Define the upper and lower boundary conditions.
      T_S = xTsfc
      bc_S = aS(1) * T_S
      bc_N = 0.0
      aN(JJ) = 0.0

c Provide the source terms.
c Build an array of values on the c.v. boundaries.
      xynet(1) = up(1) - down(1)
      do j=2,JJ
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + down(j+1))/2.0
      enddo
      xynet(JJ+1) = up(JJ+2) - down(JJ+2)

c Force the source terms to produce Tf at the positions with
c   water.  Let the solar radiation exist in other regions.
      do j=1,JJ
        if (T_old(j).eq.Tf) then
          Sc(j) = 10e30 * Tf
          Sp(j) = -10e30
        else

c This is dq/dz for q=Q0*exp(-extcoef*z).
c         Sc(j) = Qsip * extcoef * exp(- extcoef * y_crds(j+1))

c This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
c   are scaled by the ratio of the incoming solar to the solar
c   used to get up and down.
c         Sc(j) = - (xynet(j) - xynet(j+1)) / dy_p(j)
          Sc(j) = - Qsi / total_solar *
     &      (xynet(j) - xynet(j+1)) / dy_p(j)

          Sp(j) = 0.0
        endif
      end do

c Start the temperature computation over from the previous time
c   step.
      do j=1,JJ
        T_old(j) = T_tmp(j)
      end do

c Configure the information for the matrix solver.
      CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,
     &  dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

c Solve the system of equations.
      CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,JJ)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLARPEN(Qsip,Qsi,albedo)

      Qsip = (1.0 - albedo) * Qsi

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ,
     &  Tf,water_frac,gamma,xk_snow,water_depth_old,temp_ice_init_C,
     &  deltaz,ro_snow_z)

      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)
      real T_old(JJ+2)
      real gamma(JJ+2)
      real water_frac(JJ+2)
	  real ro_snow_z(JJ+2)

c Provide values of Control Volume size in the y direction, and
c   compute c.v. size information.
      CALL GETCV(deltaz,dy_p,JJ)
      CALL CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

c Supply the initial conditions.
      do j=1,JJ
        T_old(j) = temp_ice_init_C + Tf
        water_frac(j) = 0.0
        gamma(j) = xk_snow
		ro_snow_z(j) = ro_snow
      end do

      water_depth_old = 0.0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water)

      real gamma(JJ+2)
      real g_b_ns(JJ+2)
      real f_n(JJ+2)
      real aN(JJ+2)
      real aP0(JJ+2)
      real aS(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)
      real b_vector(JJ+2)
      real T_old(JJ+2)
      real Sc(JJ+2)
      real Sp(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      real xynet(JJ+2)

c Compute the solar radiation penetrating the lake surface.
      CALL SOLARPEN(Qsip,Qsi,albedo)

      xTsfc = Tsfc

c Compute gamma and the general equation coefficients.
      CAll GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water)
      CALL GAMMA1(g_b_ns,gamma,f_n,JJ)
      CALL GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ,
     &  ro_snow,Cp_snow)

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c Account for the boundary conditions.
c   South boundary condition:
c     For T_S = known, define 
c       bc_S = aS(1) * T_S;         where T_S = known
c     For dT_S/dn = 0, define
c       bc_S = 0.0
c       aS(1) = 0.0
c   North boundary condition:
c     For T_N = known, define 
c       bc_N = aN(JJ) * T_N;        where T_N = known
c     For dT_N/dn = 0, define
c       bc_N = 0.0
c       aN(JJ) = 0.0
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      T_S = xTsfc
      bc_S = aS(1) * T_S
      bc_N = 0.0
      aN(JJ) = 0.0

c Provide the source terms.
c Build an array of values on the c.v. boundaries.
      xynet(1) = up(1) - down(1)
      do j=2,JJ
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + down(j+1))/2.0
      enddo
      xynet(JJ+1) = up(JJ+2) - down(JJ+2)

      do j=1,JJ

c This is dq/dz for q=Q0*exp(-extcoef*z).
c       Sc(j) = Qsip * extcoef * exp(- extcoef * y_crds(j+1))

c This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
c   are scaled by the ratio of the incoming solar to the solar
c   used to get up and down.
c       Sc(j) = - (xynet(j) - xynet(j+1)) / dy_p(j)
        Sc(j) = - Qsi / total_solar * (xynet(j) - xynet(j+1)) / dy_p(j)

!c     print *,j,Sc(j)

        Sp(j) = 0.0
      end do

c Configure the information for the matrix solver.
      CALL PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,
     &  dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

c Solve the system of equations.
      CALL TRISOLVE(T_old,A_sub,A_main,A_super,b_vector,JJ)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETCV(deltaz,dy_p,JJ)

      real dy_p(JJ)

c Provide values of Control Volume size in the y direction.
      do j=1,JJ
        dy_p(j) = deltaz
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water)

      real gamma(JJ)
      real water_frac(JJ)
      real T_old(JJ)

c Over water.
c     A = 6.1121 * 100.0
c     B = 17.502
c     C = 240.97
c Since this is within the ice, I am using ice coefficients,
c   eventhough at higher temperatures there is a mix of water and
c   ice.
c Over ice.
      A = 6.1115 * 100.0
      B = 22.452
      C = 272.55

      do j=1,JJ
        esi = A * exp((B * (T_old(j) - Tf))/(C + (T_old(j) - Tf)))
        desi_dT = esi * B * C/((C + (T_old(j) - Tf))**2)

        De = 9.0e-5 * (T_old(j) / Tf)**14

c Increase the water vapor diffusion coeff by a factor of 10 to
c   be consistent with Anderson (1976).
c       xk_vapor = xLs * De / (Rv * T_old(j)) * desi_dT
        xk_vapor = 10.0 * xLs * De / (Rv * T_old(j)) * desi_dT

c This is a change from the original paper.  Here I am assuming
c   that the snow/ice density determines the air and ice/water
c   fractions.
        air_frac = 1.0 - ro_snow/ro_water
        gamma(j) = (1.0-air_frac) *
     &    ((1.0-water_frac(j)) * xk_snow + water_frac(j) * xk_water) +
     &    air_frac * xk_vapor

c       gamma(j) = (1.0 - water_frac(j)) * xk_snow +
c    &    water_frac(j) * xk_water + xk_vapor
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PREPSOLVE(A_sub,A_super,A_main,b_vector,T_old,
     &  dy_p,bc_S,bc_N,Sc,Sp,aN,aS,aP0,JJ)

      real aP(JJ+2)
      real aN(JJ+2)
      real aS(JJ+2)
      real Sp(JJ+2)
      real Sc(JJ+2)
      real aP0(JJ+2)
      real dy_p(JJ+2)
      real T_old(JJ+2)
      real b_vector(JJ+2)
      real A_sub(JJ+2)
      real A_super(JJ+2)
      real A_main(JJ+2)

c Compute matrix diagonal and b coeffs.
      do j=1,JJ
        aP(j) = aN(j) + aS(j) + aP0(j) - Sp(j) * dy_p(j)
        b_vector(j) = Sc(j) * dy_p(j) + aP0(j) * T_old(j)
      end do

c Modify b to account for dirichlet boundary conditions.
      b_vector(1) = b_vector(1) + bc_S
      b_vector(JJ) = b_vector(JJ) + bc_N

c Prepare to call the tridiagonal solver.
      do j=1,JJ-1
        A_sub(j) = - aS(j+1)
        A_super(j) = - aN(j)
      end do

      do j=1,JJ
        A_main(j) = aP(j)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

      real dy_pbc(JJ+2)
      real dely_p(JJ+2)
      real f_n(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)

c PRESSURE CONTROL VOLUME SIZE AND POSITION INFORMATION

c Include exterior boundary pressure grid points.
      dy_pbc(1) = 0.0
      do j = 2,JJ+1
        dy_pbc(j) = dy_p(j-1)
      end do
      dy_pbc(JJ+2) = 0.0

c Compute the distance between pressure grid points.
      do j = 1,JJ+1
        dely_p(j) = .5 * (dy_pbc(j) + dy_pbc(j+1))
      end do

c Compute the distance between the pressure grid points and the control
c   volume wall.  (The following is true because the grid points do
c   pressure are defined to be in the center of the control volume.)
c   And then compute f_e and f_n.  These two steps are combined below.
      do j = 1,JJ+1
        f_n(j) = .5 * dy_pbc(j+1) / dely_p(j)
      end do

c Compute the x and y coordinates of the pressure c.v. grid points,
c   including boundaries.
      temp = 0.0
      do j = 1,JJ+2
        y_crds(j) = temp + .5 * dy_pbc(j)
        temp = temp + dy_pbc(j)
      end do

c Compute the x and y coordinates of the pressure c.v. walls.
      y_wall(1) = 0.0
      do j = 2,JJ+1
        y_wall(j) = y_wall(j-1) + dy_p(j-1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GAMMA1(g_b_ns,gamma,f_n,JJ)

      real g_b_ns(JJ+1)
      real gamma(JJ)
      real g_ns(JJ+2)
      real f_n(JJ+2)

c This provides gamma information on c.v. walls.

c Include gamma just outside of n, s boundaries.
      g_ns(1) = gamma(1)
      do j = 2,JJ+1
        g_ns(j) = gamma(j-1)
      end do
      g_ns(JJ+2) = gamma(JJ)

c Compute gamma (diffusion coefficient) at the n, s control
c   volume boundaries using equation 4.9, p. 45.
      do j=1,JJ+1
        g_b_ns(j) = 1.0/((1.0 - f_n(j))/g_ns(j) + f_n(j)/g_ns(j+1))
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GE_COEF(aN,aS,aP0,dy_p,dely_p,g_b_ns,dt,JJ,
     &  ro_snow,Cp_snow)

      real aN(JJ+2)
      real aS(JJ+2)
      real aP0(JJ+2)
      real dely_p(JJ+2)
      real g_b_ns(JJ+2)
      real dy_p(JJ+2)

c CALCULATE THE COEFFICIENTS aP, for the general phi equation.
      do j = 2,JJ+1
        aN(j-1) = g_b_ns(j)   / dely_p(j)
        aS(j-1) = g_b_ns(j-1) / dely_p(j-1)
      end do

      do j=1,JJ
        aP0(j) = ro_snow * Cp_snow * dy_p(j) / dt
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc solve the system of equations to get the up/down flux at each depth
c mgc NOTE: first argument, x, is 'rad' in GETUPDOWN and SOLVETWOSTREAM
c SOLVETWOSTREAM call: TRISOLVE(rad,A_sub,A_main,A_super,b_vector,nz)

      SUBROUTINE TRISOLVE(x,asub,amain,asuper,b,JJ)

      real asub(JJ+2)
      real asuper(JJ+2)
      real amain(JJ+2)
      real b(JJ+2)
      real x(JJ+2)
      real z(JJ+2)
      real lmain(JJ+2)
      real lsub(JJ+2)
      real usuper(JJ+2)

      lmain(1) = amain(1)
      usuper(1) = asuper(1)/lmain(1)

      do j=2,JJ-1
        lsub(j-1) = asub(j-1)
        lmain(j) = amain(j) - lsub(j-1) * usuper(j-1)
        usuper(j) = asuper(j) / lmain(j)
      end do

      lsub(JJ-1) = asub(JJ-1)
      lmain(JJ) = amain(JJ) - lsub(JJ-1) * usuper(JJ-1)
      z(1) = b(1) / lmain(1)

      do j=2,JJ
        z(j) = 1.0 / lmain(j) * (b(j) - lsub(j-1) * z(j-1))
      end do

      x(JJ) = z(JJ)

      do j=JJ-1,1,-1
        x(j) = z(j) - usuper(j) * x(j+1)
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ENERGY-BALANCE SECTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ENBALANCE(Tair,windspd,rh,
     &    Tsfc,Qsi,Qli,Qle,Qh,
     &    Qe,Qc,Qm,balance,Qf,
     &    swe_depth,topo,z_windobs,dt,gamma,
     &    T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0,
     &    J_day_start,iter,xlat,slope_az,terrain_slope,
     &    transmiss,ihour)

c These are here to handle the case of non-zero conduction.
      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)

c Define the constants used in the computations.
        CALL CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity,xLs,
     &    xkappa,xLf,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice,
     &    ihrs_day)

c Atmospheric vapor pressure from relative humidity data.
        CALL VAPPRESS(ea,rh,Tair,Tf)

c Compute the average station pressure.
        CALL PRESSURE(Pa,one_atmos,scale_ht,topo)

c Compute the incoming longwave radiation.
        CALL LONGIN(Qli,ea,Tair,Stef_Boltz)

c Compute the incoming solar radiation.
        CALL SOLARIN(J_day_start,iter,dt,Qsi,xlat,cloud_frac,
     &    slope_az,terrain_slope,ihrs_day,transmiss,ihour)

c Compute the turbulent exchange coefficients.
        CALL EXCOEFS(De_h,z_0,z_windobs,windspd,xkappa)

c Compute the flux contribution due to conduction.
        CALL CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ)

c Solve the energy balance for the surface temperature.
        CALL SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,
     &    De_h,Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,
     &    Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc)

c Make sure the snow surface temperature is <= 0 C.
        CALL MELTTEMP(Tsfc,Tf,swe_depth)

c Compute the stability function.
        CALL STABLEFN(stability,Tair,Tsfc,windspd,z_windobs,
     &    gravity,xkappa,z_0)

c Compute the water vapor pressure at the surface.
        CALL VAPOR(es0,Tsfc,Tf)

c Compute the latent heat flux.
        CALL LATENT(Qe,De_h,stability,ea,es0,ro_air,
     &    xLs,Pa)

c Compute the sensible heat flux.
        CALL SENSIBLE(Qh,De_h,stability,Tair,Tsfc,
     &    ro_air,Cp)

c Compute the longwave flux emitted by the surface.
        CALL LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

c Compute the energy flux available for melting or freezing.
        CALL MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,
     &    Qc,Qm,Qf,Tsfc,Tf,Tair,windspd,z_windobs,
     &    gravity,De_h,ea,ro_air,xLs,Pa,Cp,emiss_sfc,
     &    Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p,
     &    JJ,icond_flag)

c Decrease the swe depth by the swe melt depth.
c   Turn this off for blue-ice simulations.

c mgc turned this back on and it worked but need to set swe_depth
C       CALL SNOW_UPDATE(swe_depth,Qm,dt,ro_ice,xLf)

c Perform an energy balance check.
        CALL ENBAL(balance,albedo,Qsi,Qli,Qle,Qh,
     &    Qe,Qc,Qm)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LONGIN(Qli,ea,Ta,Stef_Boltz)

c Compute Qli.
      emiss_cloud = 1.08 * (1.0 - exp(-(0.01 * ea)**(Ta/2016.)))
      Qli = emiss_cloud * Stef_Boltz * Ta**4

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLARIN(J_day_start,iter,dt,Qsi,xlat,cloud_frac,
     &  slope_az,terrain_slope,ihrs_day,transmiss,ihour)

c Compute the incoming solar radiation.  Here I am going to assume
c   that we are using daily time steps.  If you have some other time
c   step, see MicroMet for ideas about what to do.
      J_day = iter + J_day_start - 1
      Qsi_sum = 0.0
      if (dt.eq.86400.0) then
        do ihour=1,ihrs_day
          xhour = real(ihour)
          call solar_rad(Qsi_tmp,J_day,xlat,cloud_frac,
     &      xhour,slope_az,terrain_slope,transmiss)
            Qsi_sum = Qsi_sum + Qsi_tmp
        enddo
        Qsi = Qsi_sum / real(ihrs_day)
      else
        J_day = iter/24 + J_day_start - 1
        xhour = real(ihour)
!        print *,J_day,xhour
        call solar_rad(Qsi,J_day,xlat,cloud_frac,
     &    xhour,slope_az,terrain_slope,transmiss)
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLAR_RAD(Qsi,J_day,xlat,cloud_frac,
     &  xhour,slope_az,terrain_slope,transmiss)

      implicit none

      integer J_day

      real solar_const,days_yr,Trop_Can,solstice,pi,deg2rad,
     &  Qsi_direct,Qsi_diffuse,cos_i,cos_Z,Qsi,xlat,sin_z,xhour,
     &  cloud_frac,slope_az,terrain_slope,sol_dec,hr_angl,
     &  trans_direct,trans_diffuse,Qsi_trans_dir,Qsi_trans_dif,
     &  sun_azimuth,slope_az_S0,transmiss

c Required constants.
      solar_const = 1370.
      days_yr = 365.25
      Trop_Can = 0.41
      solstice = 173.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

c COMPUTE THE BASIC SOLAR RADIATION PARAMETERS.

c Compute the solar declination angle (radians).
      sol_dec = Trop_Can *
     &  cos(2.*pi * (real(J_day) - solstice)/days_yr)
      
c Compute the sun's hour angle (radians).
      hr_angl = (xhour * 15.0 - 180.0) * deg2rad

c Compute cos_Z.  Note that the sin of the solar elevation angle,
c   sin_alfa, is equal to the cosine of the solar zenith angle,
c   cos_Z.
      cos_Z = sin(sol_dec) * sin(xlat * deg2rad) + 
     &  cos(sol_dec) * cos(xlat * deg2rad) * cos(hr_angl)
      cos_Z = max(0.0,cos_Z)

c Account for clouds, water vapor, pollution, etc.
      trans_direct = transmiss * (0.6 + 0.2 * cos_Z) * (1.0-cloud_frac)
      trans_diffuse = transmiss * (0.3 + 0.1 * cos_Z) * cloud_frac

c Compute the solar radiation transmitted through the atmosphere.
      Qsi_trans_dir = solar_const * trans_direct
      Qsi_trans_dif = solar_const * trans_diffuse

c COMPUTE THE CORRECTIONS TO ALLOW FOR TOPOGRAPHIC SLOPE AND ASPECT.

c The sine of the solar zenith angle.
      sin_Z = sqrt(1.0 - cos_Z*cos_Z)

c Azimuth of the sun, with south having zero azimuth.
      sun_azimuth = 
     &  asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))

c Make the corrections so that the angles below the local horizon
c   are still measured from the normal to the slope.
      if (hr_angl.lt.0.0) then
        if (hr_angl.lt.sun_azimuth) sun_azimuth = - pi - sun_azimuth
      elseif (hr_angl.gt.0.0) then
        if (hr_angl.gt.sun_azimuth) sun_azimuth = pi - sun_azimuth
      endif

c Build, from the variable with north having zero azimuth, a 
c   slope_azimuth value with south having zero azimuth.
      if (slope_az.ge.180.0) then
        slope_az_S0 = slope_az - 180.0
      else
        slope_az_S0 = slope_az + 180.0
      endif

c Compute the angle between the normal to the slope and the angle
c   at which the direct solar radiation impinges on the sloping
c   terrain (radians).
      cos_i = cos(terrain_slope * deg2rad) * cos_Z + 
     &  sin(terrain_slope * deg2rad) * sin_Z * 
     &  cos(sun_azimuth - slope_az_S0 * deg2rad)

c Adjust the topographic correction due to local slope so that
c   the correction is zero if the sun is below the local horizon 
c   (i.e., the slope is in the shade) or if the sun is below the
c   global horizon.
      if (cos_i.lt.0.0) cos_i = 0.0
      if (cos_Z.le.0.0) cos_i = 0.0

c Adjust the solar radiation for slope, etc.
      Qsi_direct = cos_i * Qsi_trans_dir
      Qsi_diffuse = cos_Z * Qsi_trans_dif

c Combine the direct and diffuse solar components.
      Qsi = Qsi_direct + Qsi_diffuse

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE PRESSURE(Pa,one_atmos,scale_ht,topo)

c Compute the average station pressure.
      Pa = one_atmos * exp(- topo / scale_ht)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SNOW_UPDATE(swe_depth,Qm,dt,ro_ice,xLf)

c Calculate the swe melt depth for this time step.
      swe_melt = Qm * dt / (ro_ice * xLf)

c Decrease the swe depth by the swe melt depth.
      swe_depth = swe_depth - swe_melt
      swe_depth = max(0.0,swe_depth)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE VAPPRESS(ea,rh,Tair,Tf)

c Also see the VAPOR subroutine.

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

c Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
c Over ice.
c       A = 6.1115 * 100.0
c       B = 22.452
c       C = 272.55

c Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc,
     &  Tf,Tair,windspd,z_windobs,gravity,De_h,ea,ro_air,xLs,Pa,Cp,
     &  emiss_sfc,Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p,
     &  JJ,icond_flag)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)

c If Qm is > 0, then this is the energy available for melting.
c   If Qm is < 0, then this is the energy available for freezing
c   liquid water in the snowpack.
      if (swe_depth.gt.0.0 .and. Tsfc.eq.Tf) then
        Qm = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc
      else
        Qm = 0.0
      endif

      if (Tsfc.lt.Tf) then
        xTsfc = Tf
        CALL STABLEFN(xstability,Tair,xTsfc,windspd,z_windobs,
     &    gravity,xkappa,z_0)
        CALL VAPOR(xes0,xTsfc,Tf)
        CALL LATENT(xQe,De_h,xstability,ea,xes0,ro_air,xLs,Pa)
        CALL SENSIBLE(xQh,De_h,xstability,Tair,xTsfc,ro_air,Cp)
        CALL LONGOUT(xQle,xTsfc,emiss_sfc,Stef_Boltz)
        CALL CONDUCT(icond_flag,xQc,gamma,T_old,dy_p,JJ)
        Qf = (1.0-albedo) * Qsi + Qli + xQle + xQh + xQe + xQc
      else
        Qf = 0.0
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MELTTEMP(Tsfc,Tf,swe_depth)

      if (swe_depth.gt.0.0 .and. Tsfc.gt.Tf) Tsfc = Tf

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)

      if (icond_flag.eq.0) then
        Qc = 0.0
      else
        Qc = - (gamma(1) + gamma(2))/2.0 * (T_old(1) - T_old(2)) /
     &    (dy_p(1) + dy_p(2))/2.0
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity,xLs,
     &  xkappa,xLf,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice,
     &  ihrs_day)

      emiss_sfc = 0.98
      Stef_Boltz = 5.6696e-8
      ro_air = 1.275
      Cp = 1004.
      gravity = 9.81
      xLs = 2.500e6
      xkappa = 0.4
      xLf = 3.34e5
      Tf = 273.16
      ro_water = 1000.0
      one_atmos = 101300.0
      scale_ht = 8500.0
      Cp_water = 4180.0
      ro_ice = 917.0
      ihrs_day = 24

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ENBAL(balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm)

      balance = (1.0-albedo) * Qsi + Qli + Qle + Qh + Qe + Qc - Qm 

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE STABLEFN(stability,Tair,Tsfc,windspd,z_windobs,
     &  gravity,xkappa,z_0)

      C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0)))**2 *
     &  sqrt(z_windobs/z_0)
      C2 = gravity * z_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      if (Tsfc.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(Tsfc - Tair)
        stability = 1.0 + B1 * (Tsfc - Tair) / B3
      elseif (Tsfc.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - Tsfc))**2)
      else
c Neutrally stable case.
        stability = 1.0
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

      Qle = - emiss_sfc * Stef_Boltz * Tsfc**4

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SENSIBLE(Qh,De_h,stability,Tair,Tsfc,ro_air,Cp)

      Qh = ro_air * Cp * De_h * stability * (Tair - Tsfc)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,
     &  Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity,
     &  xLs,xkappa,z_0,Tf,Qc)

      AAA = ro_air * Cp * De_h
      CCC = 0.622 / Pa
      DDD = emiss_sfc * Stef_Boltz
      EEE = (1.0-albedo) * Qsi + Qli + Qc
      FFF = ro_air * xLs * De_h

c Compute the constants used in the stability coefficient
c   computations.
      C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0)))**2 *
     &  sqrt(z_windobs/z_0)
      C2 = gravity * z_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      CALL SOLVE(Tsfc,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVE(xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf)

      tol = 1.0e-4
      maxiter = 20
      old = Tair

c Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
c Over ice.
c       A = 6.1115 * 100.0
c       B = 22.452
c       C = 272.55

      do i=1,maxiter

c This section accounts for an increase in turbulent fluxes
c   under unstable conditions.
        other1 = AAA * (Tair - old)
        es0 = A * exp((B * (old - Tf))/(C + (old - Tf)))
        other2 = FFF*CCC*(ea-es0)

        dother1 = - AAA
        dother2 = - FFF*CCC*es0*B*C/((C + (old - Tf))**2)

      if (old.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(old - Tair)
        stability = 1.0 + B1 * (old - Tair) / B3
        dstability = B1/B3 - (B1*B2*(old-Tair))/
     &    (2.0*B3*B3*sqrt(old-Tair))
        fprime1 = - 4.0*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      elseif (old.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - old))**2)
        dstability = 2.0 * B8 / ((1.0 + B8 * (Tair - old))**3)
        fprime1 = - 4.0*DDD*old**3
        fprime2 = stability * dother1 + other1 * dstability
        fprime3 = stability * dother2 + other2 * dstability
        fprime4 = - 0.0

      else
c Neutrally stable case.
        stability = 1.0
        fprime1 = - 4.0*DDD*old**3
        fprime2 = dother1
        fprime3 = dother2
        fprime4 = - 0.0
      endif

        funct = EEE - DDD*old**4 + AAA*(Tair-old)*stability +
     &    FFF*CCC*(ea-es0)*stability +
     &    0.0
        fprime = fprime1 + fprime2 + fprime3 + fprime4

        xnew = old - funct/fprime

        if (abs(xnew - old).lt.tol) return
        old = xnew

      end do

c mgc commented next two lines out
c      write (*,102)
c  102 format('max iteration exceeded when solving for Tsfc')
      xnew = Tair

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LATENT(Qe,De_h,stability,ea,es0,ro_air,xLs,Pa)

      Qe = ro_air * xLs * De_h * stability * (0.622/Pa * (ea - es0))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EXCOEFS(De_h,z_0,z_windobs,windspd,xkappa)

      De_h = (xkappa**2) * windspd / ((log(z_windobs/z_0))**2)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE VAPOR(es0,Tsfc,Tf)

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

c Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
c Over ice.
c       A = 6.1115 * 100.0
c       B = 22.452
c       C = 272.55

c Compute the water vapor pressure at the surface.
      es0 = A * exp((B * (Tsfc - Tf))/(C + (Tsfc - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c EXTINCTION COEFFICIENT SECTION 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc define the grain radius, the 118x47 bands x grainsize grid, and 
c call GETSCATTERCOEFS, GETDWAVELEN, GETSOLAR, GETZ, SPECTEXTCOEF, 
c BULKEXTCOEF, GET_A_AND_R, and SOLVETWOSTREAM
c mgc NOTE: the reason the 118x47 grid is read in is because 
c SPECTEXTCOEF uses (k,n_snowgrain_radius) to index into each of the
c mie parameters in mie.dat, which are each 118x47

      SUBROUTINE EXTCOEFS(nz,deltaz,albedo,ro_snow,up,down,
     &  n_snowgrain_radius,total_solar)

c This model is described by Equations 7-14 in the paper:
c   Below-surface ice melt on the coastal Antarctic ice sheet, by
c   Glen E. Liston and 4 others, Journal of Glaciology, 1999,
c   Vol. 45, No. 150, pages 273-285.
c
c The author of this code is:
c   Dr. Glen E. Liston
c   InterWorks Consulting
c   15048 NCR 25E
c   Loveland, Colorado 80538
c
c   Voice: (970) 491-8220
c   Internet: liston@iceberg.atmos.colostate.edu
c
c The model should run with any FORTRAN 77 compiler.
c
c Note that the outputs correspond to the center of each level,
c   and the top surface of the top grid cell, and the bottom
c   surface of the bottom grid cell (thus the 502 values, for
c   the current setup).

c There should be no limit to how coarse (like 10 cm) or how fine
c   (like 1 mm) you want to define your levels to be.  There is
c   an upper limit of 500 levels hard coded in the program that
c   could be changed if you want/need to.
c
c The simulation domain is defined by using a combination of deltaz
c   (the thickness of each level) and nz (the number of levels).
c
c In the GETSOLAR subroutine you will find where I use Jerry
c   Harrington's solar radiation spectrum.  If you read in observed
c   radiation values and corresponding wavelengths, this code
c   will do the interpolation to the 118 bands used by the model.
c
c I suspect that the code could be modified to have different
c   densities and grain sizes at different levels, but I have not
c   looked into this.  Right now it assumes that these are the same
c   throughout the vertical domain.
c
c I talked with someone the other day who suggested that the albedo
c   used in the model could (maybe should!) be made to be wavelength
c   dependent.
c
c Note that: down = solar down, up = solar up, up/down = albedo,
c   and up-down = net solar (or - solar_absorbed).

      implicit none

      integer nvalues,nclasses,nz,n_snowgrain_radius

c The number of wavelength bands that are used.
      parameter(nvalues=118)

c The number of grain radii that can be used.
      parameter(nclasses=47)

      real deltaz,ro_snow,r_snow,Qsi,albedo,total_solar


      real radii(nclasses)
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real wavelength(nvalues,nclasses)

      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)

      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)

      real z_without_bc(nz)
      real z_with_bc(nz+2)

      real down(nz+2)
      real up(nz+2)
      real rad(nz)

      real a(nz+2)
      real r(nz+2)

c These are the possible radii that can be used for the model
c   simulations.  They are in mm, and must be converted to meters
c   before they can be used in the model.  To pick a radius, you
c   just pick the array position corressponding to the radius
c   you are interested in.  The model will extract the value and
c   convert to meters.
      data radii/0.005, 0.007, 0.010, 0.015, 0.020,
     &           0.030, 0.040, 0.050, 0.065, 0.080,
     &           0.100, 0.120, 0.140, 0.170, 0.200,
     &           0.240, 0.290, 0.350, 0.420, 0.500,
     &           0.570, 0.660, 0.760, 0.870, 1.000,
     &           1.100, 1.250, 1.400, 1.600, 1.800,
     &           2.000, 2.500, 3.000, 3.500, 4.000,
     &           4.500, 5.000, 5.500, 6.000, 6.500,
     &           7.000, 7.500, 8.000, 8.500, 9.000,
     &           9.500,10.000/

      r_snow = radii(n_snowgrain_radius) / 1000.0

c mgc commented out next two lines
c      write (6,102) r_snow * 1000.0
c  102 format ('you have picked a grain radius (mm) = ',f9.3)

c Read in the wavelength-dependent scattering coefficient arrays.
      CALL GETSCATTERCOEFS(nvalues,nclasses,g,qext,ss_coalb,
     &  wavelength)

c Generate a delta_wavelength array.
      CALL GETDWAVELEN(wavelength,nvalues,dwavelen,nclasses)

c Produce a downward solar spectrum.
      CALL GETSOLAR(nvalues,nclasses,wavelength,solar,
     &  Qsi,dwavelen,total_solar)

c Make a z-depth arrays, in meters.
      CALL GETZ(z_with_bc,z_without_bc,deltaz,nz)

c Compute the spectral extinction coefficients as a function of
c   wavelength.
      CALL SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow,
     &  r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses)

c Compute the downward bulk extinction coefficient.
      CALL BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen,
     &  z_without_bc,spect_extcoef_snow,bulk_ext_snow,
     &  bulk_ext_snowbc)

c Compute the a and r coefficients from knowledge of the
c   surface albedo and the extinction coefficient.
      CALL GET_A_AND_R(nz,bulk_ext_snowbc,a,r,albedo)

c Solve the system of equations.
      CALL SOLVETWOSTREAM(nz,a,r,deltaz,bulk_ext_snowbc,Qsi,rad,
     &  up,down)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc read in mie.dat and store g, qext, omega, and lambda

      SUBROUTINE GETSCATTERCOEFS(nvalues,nclasses,g,qext,ss_coalb,
     &  wavelength)

      implicit none

      integer nvalues,nclasses,i,j

      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real wavelength(nvalues,nclasses)

c Glen's unix binary arrays.

c Read in the wavelength-dependent scattering coefficient arrays.
c     open (46,file='mie.gdat',
c    &  form='unformatted',access='direct',recl=4*nvalues*nclasses)
c     read (46,rec=1) ((g(i,j),i=1,nvalues),j=1,nclasses)
c     read (46,rec=2) ((qext(i,j),i=1,nvalues),j=1,nclasses)
c     read (46,rec=3) ((ss_coalb(i,j),i=1,nvalues),j=1,nclasses)
c     read (46,rec=4) ((wavelength(i,j),i=1,nvalues),j=1,nclasses)

c Text version.
      open (46,file='mie.dat',form='formatted')
      do j=1,nclasses
        read (46,301) (g(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (qext(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (ss_coalb(i,j),i=1,nvalues)
      enddo
      do j=1,nclasses
        read (46,301) (wavelength(i,j),i=1,nvalues)
      enddo
  301 format(118f10.6)

      close (46)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc get the bulk absorption and reflection coefficients at each 
c layer in the z-direction 

      SUBROUTINE GET_A_AND_R(nz,bulk_ext_snowbc,a,r,albedo)

      implicit none

      integer k,nz

      real albedo
      real bulk_ext_snowbc(nz+2)
      real a(nz+2)
      real r(nz+2)

c Compute the a and r coefficients from knowledge of the
c   albedo and the bulk extinction coefficient.
      do k=1,nz+2
        a(k) = (1.0 - albedo) /
     &    (1.0 + albedo) * bulk_ext_snowbc(k)
        r(k) = 2.0 * albedo * bulk_ext_snowbc(k) /
     &    (1.0 - albedo**2)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc define a z-depth array, the centers of each grid cell

      SUBROUTINE GETZ(z_with_bc,z_without_bc,deltaz,nz)

      implicit none

      integer k,nz

      real deltaz
      real z_without_bc(nz)
      real z_with_bc(nz+2)

c Make a z-depth array, in meters.  These are z-levels, the centers
c   of each grid cell.  z_with_bc includes the top and bottom edges of
c   the top and bottom grid cells.
      z_with_bc(1) = 0.0
      do k=2,nz+1
        z_with_bc(k) = deltaz * real(k-1) - deltaz / 2.0
        z_without_bc(k-1) = z_with_bc(k)
      enddo
      z_with_bc(nz+2) = deltaz * real(nz)

c     do k=1,nz+2
c       print *, k,z_with_bc(k)
c     enddo

c     do k=1,nz
c       print *, k,z_without_bc(k)
c     enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc compute the bulk (spectrally integrated) extinction coefficient
c mgc NOTE: check the multiplication by dwavelen(k), confirm I do it
c this way also

      SUBROUTINE BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen,
     &  z_without_bc,spect_extcoef_snow,bulk_ext_snow,
     &  bulk_ext_snowbc)

      implicit none

      integer k,kk,nvalues,nz

      real sum1,sum2,deltaz
      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)
      real z_without_bc(nz)
      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)

c Compute the downward bulk extinction coefficient.
      do kk=1,nz
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nvalues
          sum1 = sum1 + solar(k) *
     &      exp(- spect_extcoef_snow(k) * (z_without_bc(kk)+deltaz)) *
     &      dwavelen(k)
          sum2 = sum2 + solar(k) *
     &      exp(- spect_extcoef_snow(k) * z_without_bc(kk)) *
     &      dwavelen(k)
        enddo
        bulk_ext_snow(kk) = - (1.0 / deltaz) * log(sum1/sum2)
      enddo

c Cast in a form that can be used in the two-stream
c   computation (add the boundaries).  Here I have assumed that
c   it is okay to call the boundaries equal to the value at the
c   center of that grid cell (the value prevails thoughout the
c   cell).
      bulk_ext_snowbc(1) = bulk_ext_snow(1) 
      do k=2,nz+1
        bulk_ext_snowbc(k) = bulk_ext_snow(k-1) 
      enddo
      bulk_ext_snowbc(nz+2) = bulk_ext_snow(nz) 

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc compute the spectral extinction coefficient
c mgc NOTE: this is where i would update N
c mgc add ro_star (removed for now)
      SUBROUTINE SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow,
     &  r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses)

      implicit none

      integer k,n_snowgrain_radius,nvalues,nclasses
c mgc add ro_star
      real ro_snow,ro_ice,sigma_e,r_snow
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real spect_extcoef_snow(nvalues)

      ro_ice = 917.0

      do k=1,nvalues
        sigma_e = 3.0/4.0 * qext(k,n_snowgrain_radius)/r_snow *
     &    ro_snow/ro_ice
        spect_extcoef_snow(k) = sigma_e *
     &    sqrt(ss_coalb(k,n_snowgrain_radius) -
     &    ss_coalb(k,n_snowgrain_radius) * g(k,n_snowgrain_radius) +
     &    ss_coalb(k,n_snowgrain_radius)**2 * g(k,n_snowgrain_radius))
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc get dlambda i.e. bandwidths

      SUBROUTINE GETDWAVELEN(wavelength,nvalues,dwavelen,nclasses)

      implicit none

      integer k,nvalues,nclasses

      real wavelength(nvalues,nclasses)
      real dwavelen(nvalues)

      dwavelen(1) = 2.0 * (wavelength(2,1) - wavelength(1,1))
      do k=2,nvalues-1
        dwavelen(k) = (wavelength(k+1,1) - wavelength(k-1,1)) / 2.0
      enddo
      dwavelen(nvalues) = 2.0 *
     &  (wavelength(nvalues,1) - wavelength(nvalues-1,1))

      do k=1,nvalues
c       print *,k,dwavelen(k),wavelength(k,1)
c       print *,wavelength(k,1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc read in solar.dat, interpolate it to the 118 bands of mie.dat, 
c and integrate it to compute total_solar

      SUBROUTINE GETSOLAR(nvalues,nclasses,wavelength,solar,
     &  Qsi,dwavelen,total_solar)

      implicit none

      integer isolarvals,k,icount,i,nvalues,nclasses

      parameter(isolarvals=250)

      real x,x1,x2,y1,y2,Qsi,total_solar
      real wavelength(nvalues,nclasses)
      real solar(nvalues)
      real dwavelen(nvalues)
      real wavel_tmp(isolarvals)
      real solar_tmp(isolarvals)

c This was just a dummy distribution used for testing.
c     data wavel_tmp/0.0,0.1,  0.4,  0.6,  1.0,  1.5, 2.0, 2.5,3.0/
c     data solar_tmp/0.0,0.0,760.0,600.0,300.0,100.0,30.0,10.0,0.0/

c Glen's unix binary arrays.

c Read Jerry Harrington's solar spectrum data.
c     open (76,file='solar.gdat',
c    &  form='unformatted',access='direct',recl=4*isolarvals)
c     read (76,rec=1) (wavel_tmp(i),i=1,isolarvals)
c     read (76,rec=2) (solar_tmp(i),i=1,isolarvals)

c Text version.
      open (76,file='solar.dat',form='formatted')
      do i=1,isolarvals
        read (76,*) wavel_tmp(i),solar_tmp(i)
      enddo

c Generate a dummy downward solar spectrum using the above _tmp
c   data strings, interpolating to the wavelengths of interest.
      do k=1,nvalues
        x = wavelength(k,1)
          do i=1,isolarvals-1
            if (x.gt.wavel_tmp(i)) then
              icount = i
            endif
          enddo
        x1 = wavel_tmp(icount)
        x2 = wavel_tmp(icount+1)
        y1 = solar_tmp(icount)
        y2 = solar_tmp(icount+1)

        solar(k) = y1 + (x - x1) * (y2 - y1)/(x2 - x1)
      enddo

c Integrate the solar radiation.
      Qsi = 0.0
      do k=1,nvalues
        Qsi = Qsi + solar(k) * dwavelen(k)
      enddo
      total_solar = Qsi
c     print *, 'total solar radiation = ',total_solar

      close (76)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc solve the system of equations to get the up/down flux at each depth

      SUBROUTINE GETUPDOWN(a,r,deltaz,nz,rad,up,down,Qsi)

      implicit none

      integer k,nz

      real alfa,deltaz,Qsi
      real down(nz+2)
      real up(nz+2)
      real down_tmp(nz+2)
      real up_tmp(nz+2)
      real rad(nz)
      real a(nz+2)
      real r(nz+2)

c Add the boundary conditions to rad.
      alfa = 1.0 / (a(1) + r(1))
      up(1) = alfa / (deltaz + alfa) * rad(1) +
     &  alfa * deltaz * r(1) * Qsi / (deltaz + alfa)
      do k=2,nz+1
        up(k) = rad(k-1)
      enddo
      up(nz+2) = 0.0

c Reconstruct y.
      down(1) = Qsi
      do k=2,nz+1
        down(k) = (a(k) + r(k)) / r(k) * up(k) - (up(k+1) - up(k-1)) /
     &    (2.0 * deltaz * r(k))
      enddo
      down(nz+2) = (a(nz+2) + r(nz+2)) / r(nz+2) * up(nz+2) - 
     &  (up(nz+2) - up(nz+1)) / (deltaz * r(nz+2))

c Smooth any small bumps in the up and down curve.  This will assist
c   in any delta up, delta down computations.
c Do the interior.
      do k=2,nz+1
        down_tmp(k) = (down(k) + 0.5 * (down(k-1) + down(k+1))) / 2.0
        up_tmp(k) = (up(k) + 0.5 * (up(k-1) + up(k+1))) / 2.0
      enddo
c Do the ends.
      down_tmp(1) = (down(2) + 0.5 * down(1)) / 1.5
      down_tmp(nz+2) = (down(nz+1) + 0.5 * down(nz+2)) / 1.5
      up_tmp(1) = (up(2) + 0.5 * up(1)) / 1.5
      up_tmp(nz+2) = (up(nz+1) + 0.5 * up(nz+2)) / 1.5
c Rewrite the arrays.
      do k=1,nz+2
        down(k) = down_tmp(k)
        up(k) = up_tmp(k)
      enddo
c Now adjust to get back the Qsi at the surface.
      do k=1,nz+2
        down(k) = down(k) * Qsi / down_tmp(1)
        up(k) = up(k) * Qsi / down_tmp(1)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc solve the two-stream model following Schlatter, 1972 method, and
c call TRISOLVE and GETUPDOWN to complete the solution

      SUBROUTINE SOLVETWOSTREAM(nz,a,r,deltaz,xmu,Qsi,rad,
     &  up,down)

      implicit none

      integer k,nz

      real alfa,xk1,gamma,deltaz,Qsi
      real b_vector(nz)
      real A_sub(nz-1)
      real A_super(nz-1)
      real A_main(nz)
      real a(nz+2)
      real r(nz+2)
      real xmu(nz+2)
      real down(nz+2)
      real up(nz+2)
      real rad(nz)

c Compute matrix diagonal and b coeffs.
      do k=1,nz
        A_main(k) = (2.0 + deltaz**2 * xmu(k+1)**2) -
     &    (deltaz / (2.0 * r(k+1)) * (a(k+1) * (r(k+2) - r(k)) -
     &    r(k+1) * (a(k+2) - a(k))))
        b_vector(k) = 0.0
      end do

c Account for the boundary conditions on the main diagonal.
      alfa = 1.0 / (a(1) + r(1))
      xk1 = 1.0 + (r(3) - r(1)) / (4.0 * r(2))
      gamma = xk1 * alfa / (deltaz + alfa)
      A_main(1) = A_main(1) - gamma

c Modify b to account for the boundary conditions.
      b_vector(1) = xk1 * alfa * r(1) * Qsi * deltaz / (deltaz + alfa)
      b_vector(nz) = 0.0

c Prepare to call the tridiagonal solver.
      do k=1,nz-1
        A_sub(k) = - (1.0 + (r(k+3) - r(k+1)) / (4.0 * r(k+2)))
        A_super(k) = - (1.0 - (r(k+2) - r(k)) / (4.0 * r(k+1)))
      end do

c Solve the system of equations.
      CALL TRISOLVE(rad,A_sub,A_main,A_super,b_vector,nz)

c Add the boundary conditions to up and reconstruct down.
      CALL GETUPDOWN(a,r,deltaz,nz,rad,up,down,Qsi)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

