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

ccccccc needed to read year on cmd line - DOS
ccccccc	USE DFLIB
c needed to call system for shell commands - DOS
c	USE DFPORT  !moved to inc file

c this inc file allows code to be compiled on dos or linux with no changes
c (the OS-specific stuff is in the include.  there are a few shell commands
c  written for each OS that will just do nothing on the other OS besides cause complaints)
	include 'os_info.inc'  

c maxiter is the number of time steps in the model run.
c      parameter (maxiter=365)
c     parameter (maxiter=1)
c replace this with leap year check below

c nz=JJ equals the number of grid cells in the z direction.  The
c   reason this is like this is because my heat equation solver
c   calls the z-dir(k) the y-dir(j).
	integer nz,JJ,nx,ny
      parameter(nz=70) !170, 37, 71, 70
      parameter(JJ=nz)
	parameter(nx=200)
	parameter(ny=140)
c	parameter(iwordlength=4)  !moved to inc file
	real deltaz(nz)
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
	real water_frac_old(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
      integer :: maxiter
	character*(2) c_yearstart
	character*(2) c_yearend
	character*(80) mm_met_file
	character*(80) stn_met_file
	character*(31) albedo_file
	character*(26) Pa_file
	character*(4) c_year
	character*(3) glaccode
	character*(80) runnametext
	character*(80) runname
	integer ablation_output
	character*(4) yeararg
	integer i_yearstart
	integer hr
	character*(80) glacier_cells_file
	character*(80) TmeanAnnual_file
	character*(80) topo_file
	character*(3) c_i
	character*(3) c_j
	character*(2) c_snowgrain_radius
	character*(9) c_z_0
!	character*(5) c_deltaz1
	character*(6) c_deltaz1
	character*(80) c_md_string
	character*(80) cjunk
	integer runcell(nx)
	logical runcond
	integer runmin, runmax
	integer cellcount
	integer topo_array(nx,ny)
	real Tannualmean(nx,ny)
	integer immstart(2,28), immoffset
	real xinternal_heating_corr
	integer :: strlen
	integer stnx(8),stny(8),glacnum
	real xdur
	real xmmdata(6,175000),xstndata(6,175000),xpadata(2,175000)
	real day_melt_abl_out(3,7000)
	integer i2,j2,iarraypos
	real xdataout(30,175000),subout(100)
	double precision totalheat,totalheat2
	character*(80) nlfname
	real endofsummerdensity(JJ)
        integer :: ierr
        real :: z_0, dz1, drainthresh, tempadd, windmult, albedo_offset
        integer :: n_snowgrain_radius

        real albedo_evo_poly_a,albedo_evo_poly_b,albedo_evo_poly_c
        real albedo_evo_poly_d,snow_albedo_thresh
c constants for calculating the best fit polynomial of albedo.
        parameter(albedo_evo_poly_a = 9.749256e-10)
        parameter(albedo_evo_poly_b = -7.432486e-07)
        parameter(albedo_evo_poly_c = 2.099152e-04)
        parameter(albedo_evo_poly_d = -0.025216)
        parameter(albedo_evo_poly_e = 1.614333)

	data stnx/53,63,65,133,143,127,164,0/
	data stny/36,43,45,93,66,88,114,0/

 
        namelist /params/ glacnum, z_0, dz1, n_snowgrain_radius,
     &                      runmin, runmax, runnametext,
     &                      tempadd, windmult, albedo_offset,
     &                      maxiter


c PROVIDE SOME OF THE RUN DEPENDENT INFORMATION.
c	nzz = nz + 1
c	nzzz = nz + 2
c	undef = -9999.0
	xdur=secnds(0.0)

c===============================================================
c Surface roughness length. 
c      z_0=0.00020
c Thickness of 'surface' layer
c	deltaz(1)=0.030
c Density, grain radius, and albedo.  See the extcoefs subroutine about how
c to determine the snow_grain_radius index.  jon=31; best fit to pyra=18
c BW1993: 100um = 11; Glen1999: snow=.35mm=18 blueice=5.0mm=37
c      n_snowgrain_radius = 16
c      n_snowgrain_radius = 10
c Glacier to run
c	glacnum = 5
c Description
c	runnametext='scalar_mmlwsw'
c===============================================================

!!c Load params from a file?
!!	if (IARGC() .eq. 1) then
!!	CALL GETARG(1, paramfname)
!!	else
!!	paramfname='icemelt.par'
!!	endif
        if (IARGC() .eq. 1) then
          CALL GETARG(1, nlfname)
        else
          nlfname='namelist.input'
        endif

!----------------------------------------------        
! Specify default values for variables read from namelist                  
!----------------------------------------------        
        glacnum = -1
        z_0 = 1.0  ! mm
        dz1 = 1.0  ! cm
        n_snowgrain_radius = 10  ! index
        runmin = -1
        runmax = -1
        runnametext = "DEFAULT_NAME"

c mgc drainage threshold
		
        drainthresh = 0.10  ! the water frac at which water is removed from subsurface
!	drainthresh = 0.99999  ! the water frac at which water is removed from subsurface
c adjustments to turn on or off for adjustments to general met data
        tempadd      = 0.0   ! FLOOR=1.5   WALL=0.5
        windmult     = 1.00  ! FLOOR=0.33  WALL=0.67
c Albedo Correction.  This value is added to measured albedo.
        albedo_offset= 0.0 ! FLOOR=-0.2  WALL=-0.05

! maxiter        
!c maxiter is the number of time steps in the model run.
!c in the hourly model it is the number of days, hours are handled later
!c consider LEAP years
!	if (MOD(i_yearend,4) .eq. 0) then
!        maxiter = 366
!	else
!        maxiter = 365
!	endif
c run the model from 1995/7/1 to 2006/6/30 = 4018
c	maxiter=4018
c run the model from 1995/7/1 to 2008/1/22 = 4589
c	maxiter=4589
c run the model from 1995/7/1 to 2008/6/30 = 4749
c	maxiter=4749
c run the model from 1995/7/1 to 2009/1/15 = 4948
c	maxiter=4948
c run the model from 1995/7/1 to 2013/1/31 = 6425
c       maxiter=6425
c run the model from 1995/7/1 to 2013/2/01 = 6426
       maxiter=6426



!----------------------------------------------        
! Read namelist
!----------------------------------------------        
         open (11, file=nlfname, status='old',form='formatted')
         read (11,params,iostat=ierr)
         if (ierr > 0) then
             write(0,*) 'Error reading namelist!'
         endif
         close(11)

! Set dz1 to the array
         deltaz(1) = dz1         


!!c runmin and runmax are necessary, but do nothing in a stn run
!!c they refer to the basins input grid and indicate what range
!!c of basin indices (inclusive) to run, which allows a run for
!!c just a particular basin or glacier.
!!	open (11,file=paramfname,form='formatted')
!!	read(11,*) glacnum, z_0, deltaz(1), n_snowgrain_radius
!!	read(11,*) runmin, runmax
!!	read(11,*) runnametext
!!	close (11)

c write hourly melt file?
	iwritehourlymelt=1

c factors to adjust met variables along cliffs
	cliffwindmult = 0.68   ! CAA westside 0.68
	clifftempadd = 0.5		!0.5


c directory to store this run
	select case (glacnum)
		case (1)
		glaccode = 'TAR'
		case (2)
		glaccode = 'TR2'
		case (3)
		glaccode = 'BFS'  ! Blood Falls Met cliff
		case (4)
		glaccode = 'CAA'
		case (5)
		glaccode = 'HOD'
		case (6)
		glaccode = 'LHC'  ! Canada cliff near LH camp
		case (7)
		glaccode = 'COH'
		case (-1)  ! cliff model run
		glaccode = 'XCL'
		case (0)   ! spatial run
		glaccode = 'XXX'
		case default
		print *,'improper glacnum!'
		stop
	end select

	if (glacnum.eq.-1) then
		iscliff=1
	else
		iscliff=0
	endif

	if (glacnum.ge.1) then
		isstn=1  ! run at a single point with stn met data
		makedetailedoutput = 1
	else
		isstn=0  ! run spatially
		makedetailedoutput = 0
	endif


c nz=170 goes to 15 m.
	if (nz.eq.170) then

		do i=2,50
		deltaz(i) = 0.01
		enddo
		do i=51,100
		deltaz(i) = 0.02
		enddo
		do i=101,130
		deltaz(i) = 0.05
		enddo
		do i=131,150
		deltaz(i) = 0.10
		enddo
		do i=151,169
		deltaz(i) = 0.50
		enddo
		do i=170,170
		deltaz(i)=0.50-deltaz(1)+0.01
		enddo

	elseif (nz.eq.71) then

		deltaz(2:20) = 0.01
		do j=21,71
			deltaz(j)=deltaz(j-1)*1.101
		enddo
		deltaz(71)=deltaz(71)-deltaz(1)+0.01

	elseif (nz.eq.70) then
c my new setup that gives 1 cm cells to 30 cm
c then gives a total of 70 cells to exactly 15 m
c the lowest cell is 1.80 m thick
c cell thickness is 0.10 m at 1.0m depth.
		deltaz(2:30) = 0.01
		do j=31,70
			deltaz(j)=deltaz(j-1)*1.13863
		enddo
		deltaz(70)=deltaz(70)-deltaz(1)+0.01

	elseif (nz.eq.37) then

	deltaz(2:15) = 0.01
	do j=16,37
		deltaz(j)=deltaz(j-1)*1.3047
	enddo
		deltaz(37)=deltaz(37)-deltaz(1)+0.01

	else
	print *,'dz not defined properly!'
	stop
	endif


ccc Range of years used in the model.  This is used for
ccc the output data files and info displayed to console.
ccc Manually change the value here to run a new year.
ccc Ideally this will be input or use a param file.
	yeararg='1995'
c	CALL GETARG(1, yeararg)
	iyr1=ichar(yeararg(1:1))
	iyr2=ichar(yeararg(2:2))
	iyr3=ichar(yeararg(3:3))
	iyr4=ichar(yeararg(4:4))
	i_yearstart=(iyr1-48)*1000+(iyr2-48)*100+(iyr3-48)*10+iyr4-48

	i_yearend=i_yearstart+1
c create a char version of yearstart for file naming purposes
	write(c_year,'(i4.4)') i_yearstart
	if (i_yearstart.ge.2000) then
		write(c_yearstart, '(i2.2)') i_yearstart-2000
	else
		write(c_yearstart, '(i2.2)') i_yearstart-1900
	endif
	if (i_yearend.ge.2000) then
		write(c_yearend, '(i2.2)') i_yearend-2000
	else
		write(c_yearend, '(i2.2)') i_yearend-1900
	endif
c	print *,'Working on year',i_yearstart

c Record number in MM met file for July 1, Hour 0 of each year
c (1994 starts at 1994/12/1/1
	data immstart / 
     &	1991, -9999,
     &	1992, -9999,
     &	1993, -9999,
     &	1994, 1,
     &	1995, 5088,
     &	1996, 13872,
     &	1997, 22632,
     &	1998, 31392,
     &	1999, 40152,
     &	2000, 48936,
     &	2001, 57696,
     &	2002, 66456,
     &	2003, 75216,
     &	2004, 84000,
     &	2005, 92760,
     &	2006, 101520,
     &	2007, 110280,
     &  2008, 119064,
     &  2009, 127824,
     &  2010, 136584,
     &  2011, 145344,
     &  2012, 154128,
     &  2013, 162888,
     &  2014, 171648,
     &  2015, 180408,
     &  2016, 189192, 
     &  2017, 197952,
     &  2018, 206712/

	immoffset=immstart(2,i_yearstart-1990)
c	immoffset=1

c Number of times to loop through the year, to ensure convergence
c  of deep ice temperatures.
      max_annual_loops = 1
c Julian day of the model start.  Usually I start the melt runs
c   in the middle of winter.
      J_day_start = 182
c     J_day_start = 1
c Height of wind and temperature observations.
      z_windobs = 3.0

c The factor to multiply by net solar rad in the energy balance to account for the fact
c that some of the energy is absorbed at depth.  Set how many layers you want to 'black out'
c and the DARKENLAYERS subroutine will calculate the equivalent energy (qsfactor).
	ndarklayers=1
c Internal Heating Correction.  Multiplies internal heat from solar rad by a factor
	xinternal_heating_corr=1.00

c Density, grain radius, and albedo.  See the extcoefs subroutine
c   about how to determine the snow_grain_radius index.
c      ro_ice = 870.0
	ro_ice = 870.0
	ro_snow = ro_ice

c initialize albedo to 0.65 for the EXTCOEFS calculation
	albedo = 0.562

c Snow-water-equivalent depth.  Any non-zero value makes the model
c   simulate snow-ice conditions.
      swe_depth = 10.0

c Model time step.  day=86400, hr=3600 sec.
      dt = 3600.0

c Latitude of center of domain, in decimal degrees.
c TAR=-77.74
      xlat = -77.74

c Fractional cloud cover and transmissivity.  These values can be
c   used to adjust the incoming solar radiation.
c clear_sky is a var jon added.
      cloud_frac = 0.05
      transmiss = 1.08
	clear_sky = 0.87

c Identify whether this is run includes a non-zero conduction term
c   (icond_flag = 0 = no conduction, 1 = conduction).
c   Include conduction for Antarctic simulations.
      icond_flag = 1

c Build the run name==============================
	write(c_snowgrain_radius,'(I2.2)') n_snowgrain_radius
	write(c_z_0,'(f9.7)') z_0
c	write(c_deltaz1,'(f5.3)') deltaz(1)
	write(c_deltaz1,'(f6.4)') deltaz(1)
	i=strlen(runnametext)
	runname = glaccode//'_'//runnametext(1:i)//'_'//c_z_0//'_'//
     &	c_deltaz1//'_'//c_snowgrain_radius 
	i=strlen(runname)
	if (DOS .eq. 1) then
c DOS shell command
	c_md_string='md output\' // runname // ''
	print *, c_md_string
	junk= system(c_md_string)
	c_md_string='cd output\' // runname(1:i)// '*.*'
	junk= system(c_md_string)
	c_md_string='copy icemelt_8hr_Qc_spatial_xy.f .\output\' 
     &	// runname(1:i)// '\'
	junk= system(c_md_string)
	else
c linux shell command
	c_md_string='mkdir output/' // runname
	junk= system(c_md_string)
	c_md_string='rm -f ./output/' // runname(1:i) // '/*'
	print *,c_md_string
	junk= system(c_md_string)
	endif

c Get the general constants to be used.
        CALL CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water,
     &    xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow,ro_pure_ice)

c Supply the initial configurations for the ice-free lake model.
        CALL ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ,
     &    Tf,water_frac,gamma,xk_snow,water_depth_old,
     &    temp_ice_init_array,deltaz)

c Calculate the solar radiation within the snow/ice matrix.
c Run extcoefs, or use last run's extcoefs results?
	i_run_extcoefs=1
	if (i_run_extcoefs.eq.1) then
        CALL EXTCOEFS(nz,deltaz,albedo,ro_snow,up,down,
     &    n_snowgrain_radius,total_solar,ro_pure_ice,y_crds,runname)
	else
		open (101,file='./output/downupext.out' )
		read (101,*) total_solar
		do nnn=1,JJ+2
			read (101,*) down(nnn),up(nnn)
c		down(pp)=down(pp)*1.000
c		up(pp)=up(pp)*1.000
		enddo
		close (101)
	endif

c Calculate how much energy blacked out in the upper n layers
	CALL DARKENLAYERS(nz,dy_p,up,down,ndarklayers,qsfactor)	

	write (6,203) z_0,qsfactor,ndarklayers
203	format ('roughness (m), qsfactor, ndarklayers = ',2f10.5,i5)

c Save qsfactor into a file for output processing
	open (81,file='./output/'//runname(1:strlen(runname))
     &	//'/qsfactor')
	write(81,'(f10.7)') qsfactor
	close(81)



c ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c START SPATIAL LOOP HERE cccccccccccccccccccccccccccccccccc
c cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Open grids needed for surface or cliff domains
	if (isstn.eq.1) then
		glacier_cells_file='./input/tv_landcovermetstake.txt'	
	else
		if (iscliff.eq.1) then
c			glacier_cells_file='./input/cliff_vegtype.txt'
			glacier_cells_file='./input/tv_basins_cliff.txt'
		else
			glacier_cells_file='./input/tv_basins_surface.txt'
		endif
	endif
	open (50,file=glacier_cells_file,form='formatted')
	
	topo_file='./input/tv_dem250.txt'
	open (51,file=topo_file,form='formatted')
	TmeanAnnual_file='./input/T_avg_all.txt' !contains both cliff & surf cells
	open (52,file=TmeanAnnual_file,form='formatted')

	iheader = 6
c read through headers
	do k=1,iheader
		read (50,*) cjunk
		read (51,*) cjunk
		read (52,*) cjunk
	enddo

c read topo array
	do j=ny,1,-1
		read (51,*) (topo_array(i,j),i=1,nx)
		read (52,*) (Tannualmean(i,j),i=1,nx)
	enddo
	close (51)
	close (52)
	TannualmeanRef = Tannualmean(stnx(1),stny(1)) !TAR stn is ref

	cellcount=0

c   START SPATIAL LOOP
	
	do jjj=ny,1,-1
		read(50,*) runcell
		do iii=1,nx

		
c check if we should run this cell================================
	if (isstn.eq.0) then
		runcond = ( (runcell(iii).ge.runmin) .and. 
     &		(runcell(iii).le.runmax) )
c	  runcond = ( (runcell(iii).eq.24).or.(runcell(iii).eq.25).or.
c     &  (runcell(iii).eq.27).or.(runcell(iii).eq.28) ) 
	else
	  runcond = ((iii.eq.stnx(glacnum)).and.(jjj.eq.stny(glacnum)))
	endif

	if (runcond) then
		cellcount=cellcount+1
        print *,'WORKING ON CELL = ',iii, ' , ' , jjj,',num',cellcount

        ! We want to use a lower albedo for Howard and Commonwealth glaciers
        ! Set that here, but only for 'clean' ice
        if (runcell(iii)>=50) then
           if (albedo_offset == 0.0) then
                albedo_offset = -0.05
           endif
        endif

c reset everything for each cell to use clean

c MM should be accounting for slope and aspect for radiation
c Use slope=0 when doing a MM run.
c but we still need elevation to calculate Pa
		slope_az = 0.0
		terrain_slope = 0.0
		topo = topo_array(iii,jjj)

c Output options can be put here.  Maybe later in a param file
c this needs to start off, so we only write on the last iteration	
	ablation_output=0

c Open the atmospheric forcing data data files for the year.
c Make sure start date is July 1st in data file!
c output name indicates location in ascii grid dimensions
	write(c_i,'(i3.3)') iii
	write(c_j,'(i3.3)') jjj


c Input files.==============================================

c Data out of MicroMet is binary and has 6 variables
c Note the wordlength (4) is dependent on compiler settings!
	if (iscliff.eq.1) then
	mm_met_file='./input_cliff/' //   c_i // c_j // '.bin'
	else
	mm_met_file='./input/' //   c_i // c_j // '.bin'
	endif
c	open (31,file=mm_met_file,access='direct',form='unformatted',
c     &  recl=iwordlength*6)
	open (31,file=mm_met_file,access='direct',form='unformatted',
     &  recl=iwordlength*6*(immoffset-1+maxiter*24))
c Read entirety of binary input files into memory to speed program
	read (31,rec=1) ((xmmdata(i2,j2),i2=1,6)
     &	,j2=1,immoffset-1+maxiter*24)

c For runs of Station data, we also need this (in same structure as mm file)
c	stn_met_file='./input/TAR_stn.bin'
	if ((isstn.eq.1).and.(glacnum.ne.2)) then
	stn_met_file='./input/'//glaccode//'_stn.bin'
c	open (32,file=stn_met_file,access='direct',form='unformatted',
c     &	recl=iwordlength*6)
		open (32,file=stn_met_file,access='direct',form='unformatted',
     &		 recl=iwordlength*6*(immoffset-1+maxiter*24))
		read (32,rec=1) ((xstndata(i2,j2),i2=1,6)
     &		,j2=1,immoffset-1+maxiter*24)     
	endif

	SELECT CASE (glacnum)
		case (0) !running all cells
			if (runmax .ge. 30) then
c Use CAA albedo in Hoare and Fryxell basins
c Note this will not work right for runs that span the basins
				albedo_file = './input/9513_alb.CAA'
			else
c Use TAR albedo in Bonney basin
				albedo_file = './input/9513_alb.TAR'
			endif
		case (-1)
			albedo_file = './input/9513_alb.clf'
		case (2)  !at TAR2 use TAR albedo
			albedo_file = './input/9513_alb.TAR'
		case (3)
			albedo_file = './input/9509_alb.BFS'
		case (6)
			albedo_file = './input/9509_alb.BFS'
		case default
			albedo_file = './input/9509_alb.' // glaccode
	end select
	open (33,file=albedo_file,form='formatted')
c read off albedo file until we reach the start date
	do i=1,immoffset-1
		if (mod(i,24).eq.1) read (33,*) junk1,junk2,junk3,xjunk4
	enddo

c Read entire Pa file
	Pa_file = './input/hoe_pa.bin'
c	open (36,file=Pa_file,access='direct',form='unformatted',
c     &  recl=iwordlength*2)
	open (36,file=Pa_file,access='direct',form='unformatted',
     &  recl=iwordlength*2*(immoffset-1+maxiter*24))
	read (36,rec=1) ((xpadata(i2,j2),i2=1,2)
     &	,j2=1,immoffset-1+maxiter*24)  
c Elevation of Lake Hoare station
	elev_ref=77.1

c FOR NOW: I will just use a random icetempinit.txt file from a TAR run
c it shouldn't vary that much once we get to summer I hope.
		open (39,file='./input/icetempinit2008good.txt')
c Supply the initial conditions.
		do j=1,JJ
		 read (39,'(f10.4)') T_old(j)
c Shift ice temp column based on mean annual air temp
	T_old(j)=T_old(j)+Tannualmean(iii,jjj)-TannualmeanRef
		 water_frac(j) = 0.0
		 gamma(j) = xk_snow
		end do	
		close (39)	
c	print *,'Ttop,Tbottom: ',T_old(1),T_old(170)


c Output files.==============================================

c DETAILED OUTPUT - for met stations, and maybe stakes

c	if ((runcell(iii).eq.25).or.(runcell(iii).eq.28)) then
	
	if (makedetailedoutput .eq. 1) then

c	open (18,file='./output/'//runname(1:strlen(runname))// 
c     &	'/' // c_i // c_j //'.enbal')
	open (26,file='./output/' // runname(1:strlen(runname))
     &	//'/'//c_i// c_j // '.subsurf'  
     &	, access='direct',form=
     &	'unformatted',recl=iwordlength*1*100 )
c	open (26,file='./output/' // runname(1:strlen(runname))
c     &	//'/'//c_i// c_j // '.subsurfT'  
c     &	, access='direct',form=
c     &	'unformatted',recl=iwordlength*1*100 )
c	open (27,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i//c_j//'.subsurfMelt' 
c     &	, access='direct',form=
c     &	'unformatted',recl=iwordlength*1*100)
	open (28,file='./output/'//runname(1:strlen(runname))
     &	//'/'//c_i//c_j//'.out' 
     &	, access='direct',form=
     &	'unformatted',recl=iwordlength*20)
c	open (27,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i//c_j//'.drain' 
c     &	, form='formatted')

c to write all data in 1 record, use recl=iwordlength*maxiter*24*20
	endif

c GENERAL OUTPUT - for all cells
c optional file to write out melt, ablation, submelt hourly
	if (iwritehourlymelt.eq.1) then
	open (21,file='./output/'//runname(1:strlen(runname))
     &	//'/'//c_i // c_j //'.ablation.hourly'
     &	, access='direct',form=
     &	'unformatted',recl=iwordlength*3)
	endif
	open (20,file='./output/'//runname(1:strlen(runname))
     &	//'/'//c_i // c_j //'.ablation.daily'
     &	, access='direct',form=
     &	'unformatted',recl=iwordlength*3*maxiter)

c a file to write out the end of summer density profile for each year.
c since it's short, just make it ascii
	open (66, file='./output/'//runname(1:strlen(runname))
     &	//'/'//c_i // c_j //'.densityprofile')

c     	open (79,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i // c_j //'.notes' )
c     	open (81,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i // c_j //'.totalheat1' )
c     	open (82,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i // c_j //'.totalheat2' )
c     	open (83,file='./output/'//runname(1:strlen(runname))
c     &	//'/'//c_i // c_j //'.totalheat3' )

     
      water_depth_old = 0.0

c Jon added this var.  I should probably move it somewhere cleaner
	ro_snow_on_top = ro_water
c Jon added this too.  He has it read a value from the previous year stored
c in a file.  That seems overly complex for now, so for now manually set it.	
	snow_cover_depth_old = 0.0
	gravity = 9.81

c Need to add initial conditions for snow_cover_depth_old
c Jon has convoluted system to read a file just to get the last day's value

		
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                      START ANNUAL LOOP - ITERATION

c Note the annual iteration will be phased out.  So don't add anything
c important to this little section
        do kkk=1,max_annual_loops
          print *,'Annual Loop Number =',kkk
c Jon has initial conditions of each iteration come from the last day
c of the previous iteration for:
c T_old, gamma, water_frac, snow_cover_depth_old, water_depth_old
c Maybe I can ignore that.  Or maybe those values will be retained anyway
c I intend to phase out the annual loop eventually, but for now
c I will duplicate his code
	if (kkk .gt. 1) then
c I am assuming these values are retained from last day of prev. iteration
c		do j=1,nz
c			T_old
c			gamma
c			water_frac
c		enddo
c		ro_snow_on_top = ro_water    c  this is now in SNOWEVENT
		snow_cover_depth_old = 0.0
c		water_depth_old  assume this retains its value too!

c Rewind files on each annual iteration - 31 only if ascii
c		rewind(31)
c		rewind(32)
		rewind(33)
c		rewind(34)
c		rewind(35)
		rewind(36)
		rewind(37)
		rewind(38)
c		rewind(40)
c		rewind(41)

	endif


	Tsfc=xmmdata(1,immoffset)+Tf !initialize tsfc for the brent solver

c set/reset the density to ice before starting the run
	do i=1,JJ
		endofsummerdensity(i)=ro_snow
	enddo


	if (kkk.eq.max_annual_loops) then
		ablation_output = 1
	endif
	
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc			START DAILY TIMESTEP LOOP
c Note: The grid cell loop will be nested into here!
      do iter=1,maxiter
		if ((mod(iter,1000).eq.0).or.(iter.le.3)) then
         print *,'WORKING ON DAY =',iter
		endif

c	print *,'WORKING ON DAY =',iter

		daymelt = 0.0	
		dayablation = 0.0
		daysubdrain = 0.0

c optional write-out of ice temps
c 4749 is 7/1/2008 assuming start on 7/1/1995
c		if (iter.eq.4749) then
c			open (87,file='./output/'//runname(1:strlen(runname))// 
c     &	'/' // c_i // c_j //'.finalicetemp')
c			do icnt=1,JJ
c				write(87,'(f10.2)') T_old(icnt)
c			enddo
c		close (87)
c		endif
c		if (iter.eq.4530) then
c			open (87,file='./output/'//runname(1:strlen(runname))// 
c     &	'/' // c_i // c_j //'.2007_11_24icetemp')
c			do icnt=1,JJ
c				write(87,'(f10.2)') T_old(icnt)
c			enddo
c		close (87)
c		endif


c Albedo is constant for each day------------
		read (33,*) junk1,junk2,junk3,albedo
c		albedo=0.57
		albedo=albedo + albedo_offset



c ====================SNOW SECTION===========================
c 1. Check for snow presence! ------------------------

c This equation is used to determine the best fit 4th order
c polynomial of albedo during the summer months (Days 65-280) for an 
c ice-exposed surface.  (from 13 summers of TAR, with snow manually removed)
c The + 0.10 term is the increase in albedo that will trigger snow.
c calculate the Day of model Year, assuming model starts on July 1.
c This might be off by a day now and then, but that's fine
c	itemp=365*INT(real(iter)/365.25)
        itemp=365*(iter/365) !integer division will round down (FLOOR)
        doy=iter-itemp
        if ((doy .ge. 75) .and. (doy .le. 265)) then
          snow_albedo_thresh=albedo_evo_poly_a * real(doy)**4 + 
     &     albedo_evo_poly_b * real(doy)**3 +
     &     albedo_evo_poly_c * real(doy)**2 + 
     &     albedo_evo_poly_d * real(doy) + albedo_evo_poly_e + 0.10
         else
           snow_albedo_thresh=0.65 + 0.10
         endif

         if ((albedo .lt. snow_albedo_thresh).or.(iscliff.eq.1)) then
            snow_cover_depth_old = 0.0
         else
            snow_cover_depth_old = 0.05 !constant 5 cm - any positive value is fine here
         endif
c		snowdepth = 0.0 !just make it 0 instead!


c ====================END SNOW SECTION===========================


c--------------Start hourly timestep-------------------
	do hr=0,23
c	         print *,'WORKING ON HOUR =',hr
c Read the atmospheric forcing data for this time step.
c	read (31,*) Tair,rh,windspd,Qsi,Qli
c Read MM data
c		read (31,rec=(immoffset + (iter-1)*24 + hr) )
c     &	Tair,rh,windspd,winddir,Qsi,Qli
	iarraypos=immoffset+(iter-1)*24+hr
	Tair=xmmdata(1,iarraypos)
	rh=xmmdata(2,iarraypos)
	windspd=xmmdata(3,iarraypos)
	winddir=xmmdata(4,iarraypos)
	Qsi=xmmdata(5,iarraypos)
	Qli=xmmdata(6,iarraypos)

c Read Station data and use it if good
c If Stn data is bad, then use MM data
c (this could be done more efficiently, but this is easier)
	if ((isstn.eq.1).and.(glacnum.ne.2)) then
	if (xstndata(1,iarraypos).gt.-9999.0) Tair=xstndata(1,iarraypos)
	if (xstndata(2,iarraypos).gt.-9999.0) rh=xstndata(2,iarraypos)
	if (xstndata(3,iarraypos).gt.-9999.0)
     &	windspd=xstndata(3,iarraypos)
c	if (xstndata(5,iarraypos).gt.-9999.0) Qsi=xstndata(5,iarraypos)
c	if (xstndata(6,iarraypos).gt.-9999.0) Qli=xstndata(6,iarraypos)
	endif

		if (Qsi .lt. 1.0) then
			Qsi=0.0
		end if

		Tair = Tair + Tf
c Windspeed needs to be above a threshold value to avoid computational problems. 0.1
c MM already does this, but station data does not
		if (windspd .lt. 1.0) then
			windspd = 1.0
		endif
c Calc Pressure using Lk Hoare Pa measurements
c There are 36 hours in the whole 14 years with Pa missing at LH
c For those times, just use the previous time step value.
c		read (36,rec=(immoffset + (iter-1)*24 + hr) )
c     &		Pa_ref,T_ref
		Pa_ref=xpadata(1,iarraypos)		
		T_ref=xpadata(2,iarraypos)		
		if ((Pa_ref .gt. 0.0) .and. (T_ref .gt. -9000)) then
		Pa=(Pa_ref*100.0) * exp( (topo-elev_ref) / 
     &		(-287.04*0.5*(T_ref+Tf + Tair)/gravity) )
		endif

c Sensitivity adjustments here=======================
		Tair=Tair+0.0
		windspd=windspd*1.0
		Qli=Qli+0.0

c Cliff Met adjustments
	if (((iscliff.eq.1).or.(glacnum.eq.3)).or.(glacnum.eq.6)) then
		windspd = windspd * cliffwindmult
		if (Qsi.gt.50.0) then
			Tair = Tair + clifftempadd
		endif
	else
c Adjustments anywhere else
		windspd = windspd * windmult
		if (Qsi.gt.50.0) then
			Tair = Tair + tempadd
		endif

	endif

c Jon has CRYCON stuff here! I am leaving out for now.

            CALL ENBALANCE(Tair,windspd,rh,
     &        Tsfc,Qsi,Qli,Qle,Qh,
     &        Qe,Qc,Qm,balance,Qf,
     &        swe_depth,topo,z_windobs,dt,gamma,
     &        T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0,
     &        J_day_start,xlat,slope_az,terrain_slope,
     &        transmiss,clear_sky,
     &	snow_cover_depth_old,surface_melt,ro_snow_on_top,
     &	ablation,iter,xLs,xLf,i_yearstart,
     &	ablation_output,stability,hr,Qsi_fraction,
     &y_crds,f_n,dely_p,Qsip,extcoef,xk_snow,ro_snow,
     &Cp_snow,xk_water,water_frac,up,down,total_solar,
     &Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,ro_pure_ice,kkk,
     &	xinternal_heating_corr,qsfactor,ndarklayers,Pa)

c mgc ifinalcall is used after ICEHEAT to decide whether to 
c write out totalheat1. I am not sure the difference b/w totalheat0
c and totalheat1, it might be vertical resolved vs column total
     
	ifinalcall = 1 !Iterating is complete
            CALL ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &        ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf,ro_pure_ice,
     &	xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qc)

c mgc totalheat
	 
	totalheat2=totalheat
	totalheat=dble(0.0)
	do mmm=1,JJ
	totalheat=totalheat+dble((T_old(mmm)-270.0)*
     &	dble(Cp_snow)*dble(dy_p(mmm))*dble(ro_snow))
	totalheat=totalheat+dble(water_frac(mmm))*dble(dy_p(mmm))*
     &	dble(ro_snow)*dble(xLf)
	enddo
	
c	print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
c	print '(f7.2, f12.1, f15.1, 3f12.1)',dble(iter)+hr/24.0,
c     &	totalheat-totalheat2 + dble(Qc*dt) -
c     &	dble(Qsi*(1-albedo)*(1-qsfactor)*dt),
c     &     totalheat,
c     &	totalheat-totalheat2,Qc*dt,
c     &	Qsi*(1-albedo)*(1-qsfactor)*dt
	
	
c mgc water_frac
		
	do i=1,JJ !update for next time step
	water_frac_old(i) = water_frac(i)
	enddo

c mgc drainage 
	
c Calculate drainage amount
	subdrain=0.0 !total amount of water we've removed for the column
	do i=1,JJ
	if (water_frac(i).gt.drainthresh) then
		subdrain=subdrain+(water_frac(i)-drainthresh)*dy_p(i)*100.0  
     &		*ro_snow/ro_water !in cm weq
c adjust density profile by the amount that drained from each depth this hour
		endofsummerdensity(i)=endofsummerdensity(i) - 
     &		(water_frac(i)-drainthresh) * ro_snow
	endif
	enddo

c mgc surface melt and ablation
	
c calcs for daily total
	daymelt = daymelt + surface_melt
	dayablation = dayablation + ablation
	daysubdrain = daysubdrain + subdrain

c Write output here
c Only write output on final iteration, if doing multiple.**********
	if (kkk.eq.max_annual_loops) then

	iarraypos=(iter-1)*24+hr+1

	xdataout(1,iarraypos)=real(iter)+real(hr)/24
	xdataout(2,iarraypos)=Tair-Tf
	xdataout(3,iarraypos)=Tsfc-Tf
	xdataout(4,iarraypos)=Qsi
	xdataout(5,iarraypos)=Qli
	xdataout(6,iarraypos)=Qle
	xdataout(7,iarraypos)=Qh
	xdataout(8,iarraypos)=Qe
	xdataout(9,iarraypos)=Qc
	xdataout(10,iarraypos)=Qm
	xdataout(11,iarraypos)=balance
	xdataout(12,iarraypos)=albedo
	xdataout(13,iarraypos)=stability
	xdataout(14,iarraypos)=surface_melt
	xdataout(15,iarraypos)=ablation
	xdataout(16,iarraypos)=snow_cover_depth_old
	xdataout(17,iarraypos)=water_depth
	xdataout(18,iarraypos)=water_flux
c	xdataout(19,iarraypos)=rh
	xdataout(19,iarraypos)=subdrain
	xdataout(20,iarraypos)=windspd

c Write General Output - hourly? no - save it for daily totals
c	write(21,rec=iarraypos) (xdataout(i2,iarraypos),i2=14,15)
c NEED TO ADD SUBSURFACE MELT IN SOME WAY


c Write Detailed Ouput---------------
c	if (((runcell(iii).eq.25).or.(runcell(iii).eq.28))
c     &  .and.(isstn.eq.1)) then
	if (makedetailedoutput .eq. 1) then
c combine t_old and water_frac into 1 array for storage
	do k=1,100
		subout(k)=t_old(k)-Tf
		if (t_old(k).eq.Tf) subout(k) = water_frac(k)
	enddo

		write(28,rec=iarraypos) 
     &		(xdataout(i2,iarraypos),i2=1,20)
		write (26,rec=(iter-1)*24 + hr +1) (subout(k),k=1,100)




	endif !detailed output-------------

c optional print out melt, ablation and drained submelt for every hour
	if (iwritehourlymelt.eq.1) then
	write(21,rec=iarraypos) surface_melt, ablation, submelt
	endif  !hourly melt file
	

	endif !final annual loop check -> output  ******************


c mgc water_frac 
	
c REMOVE EXCESS WATER NOW that we have calculated this hour's output!
	do i=1,JJ
	if (water_frac(i).gt.drainthresh) then
		water_frac(i)=drainthresh
	endif
	enddo
    
	enddo !^----------End hourly timestep-------------^

c mgc below this he writes out daymelt, dayablation, and daysubdrain

c Write daily totals+++++++++++++++++++++++++++++++++++++++
c save to array, then write it all at once at end of run
	if (kkk.eq.max_annual_loops) then
	day_melt_abl_out(1,iter) = daymelt
	day_melt_abl_out(2,iter) = dayablation
	day_melt_abl_out(3,iter) = daysubdrain
	endif
c                  +++++++++++++++++++++++++++++++++++++++

c mgc below this he checks if it's the last day and if so, 
c writes an end-of-summer density profile for the upper 30 cm

c end of summer density profile - write out and reset at end of each summer
c do this each year on jan 31 (approx end of summer date).  
c If model starts on july 1, then jan 31 is 215 days later
	if (mod(real(iter+365-215),365.25).lt.1.0)  then
c write out just the upper 30 layers
		write(66,'(30f8.1)') (endofsummerdensity(i),i=1,30)

c no reset profile to density
		do i=1,JJ
			endofsummerdensity(i)=ro_snow
		enddo
	endif



      enddo !^---------- End Daily Loop ------------------^
      enddo !^----------- End Annual Loop -------------------^

c Write out final data as binary
c don't bother with all the stuff
c	write(28,rec=1) ((xdataout(i2,j2),i2=1,20),j2=1,maxiter*24)
c just write the 3 ablation quantities in cm
	write(20,rec=1) ((day_melt_abl_out(i2,j2),i2=1,3),j2=1,maxiter)
c write a final end of summer density in case the last season ends before jan 31
		write(66,'(30f8.1)') (endofsummerdensity(i),i=1,30)

c close files and end spatial loops
	close (18)
	close (19)
	close (20)
	close (21)
c	close (31)
	close (33)
	close (36)
	close (37)
	close (38)
	close (26)
	close (27)
	close (28)
	close (66)
	endif !check to run cell
	enddo !spatial loop
	enddo !spatial loop

	close (50)
	
c  88  format (i6,f8.1,f10.4,f8.1,f10.4,7f9.3,f14.9,f8.2)
c  89  format (f9.3,2f10.4,6f9.3,f11.5,f14.9,f6.2,f7.3)

	print *,runnametext
	print *,'z0,dz1,ngrainrad: ', z_0,deltaz(1),n_snowgrain_radius

	xdur=secnds(xdur)
	PRINT *,'ELAPSED TIME,CELLS,TIME/CELL:',
     &	xdur,cellcount,xdur/cellcount
 
      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c ICE SECTION
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ICE_ENERGY(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &  water_depth_old,water_flux,xLf,ro_pure_ice,
     &	xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qcout)

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
	integer ndarklayers
	real Qc1,Qc2,Qc3,Qcout
      real Sc(JJ+2),T1(JJ+2)

	double precision totalheat0,totalheat1,totalheat2,totalheat3
	integer icondflag

	icondflag = 1
c Save a copy of the temperature at the previous time step.  This
c   will be used to compute the ice temperature profile for the
c   case where liquid water is present, after a computation is 
c   made to account for the amount of ice melted or water frozen.
      do j=1,JJ
        T_tmp(j) = T_old(j)
      enddo

c mgc totalheat0 = sensible + latent heat at each layer 1:JJ

	totalheat0=0.0
	do mmm=1,JJ
	totalheat0=totalheat0+
     &	dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
	totalheat0=totalheat0+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
	enddo

c Solve the ice temperature equation.
      CALL ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water,ro_pure_ice,xinternal_heating_corr,
     &  ndarklayers)

c mgc added call to CONDUCT
     
	call CONDUCT(icondflag,Qc1,gamma,T_old,dy_p,JJ,Tsfc)

c mgc totalheat1, I think he is checking if there's a difference
c before/after the calls to CONDUCT
	if(ifinalcall.eq.1) then
	totalheat1=0.0
	do mmm=1,JJ
	T1(mmm)=T_old(mmm)
	totalheat1=totalheat1+
     &	dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
	totalheat1=totalheat1+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
	enddo	
c	print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
c	print '(i4, f12.1, f15.1, 3f12.1)',1,
c     &	totalheat1-totalheat0 + dble(Qc1*dt -
c     &	Qsi*(1-albedo)*(1-qsfactor)*dt),
c     &     totalheat1,
c     &	totalheat1-totalheat0,Qc1*dt,
c     &	Qsi*(1-albedo)*(1-qsfactor)*dt
c	write (81,*)  totalheat1-totalheat0,dble(Qc1)*dble(dt),
c     &	dble(Qsi)*dble(1-albedo)*dble(1-qsfactor)*dble(dt)
	endif

c Correct for ice temperatures above freezing, and compute the
c   meltwater produced by the extra available energy.  Also deal
c   with the case of refreezing water.
      CALL ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice,
     &  ro_snow,water_frac,flag,ro_pure_ice)

	call CONDUCT(icondflag,Qc2,gamma,T_old,dy_p,JJ,Tsfc)

c mgc totalheat2
	if(ifinalcall.eq.1) then
	totalheat2=0.0
	do mmm=1,JJ
	totalheat2=totalheat2+
     &	dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
	totalheat2=totalheat2+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
	enddo	
c	print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
c	print '(i4, f12.1, f15.1, 3f12.1)',2,
c     &	totalheat2-totalheat0 + dble(Qc2*dt -
c     &	Qsi*(1-albedo)*(1-qsfactor)*dt),
c     &     totalheat2,
c     &	totalheat2-totalheat0,Qc2*dt,
c     &	Qsi*(1-albedo)*(1-qsfactor)*dt
c	write (82,*)  totalheat2-totalheat0,dble(Qc2)*dble(dt),
c     &	dble(Qsi)*dble(1-albedo)*dble(1-qsfactor)*dble(dt)
	endif

c If water is present, recompute the temperature profile.  Making
c   sure the water areas are at Tf.
    
c Disable Step 3 for now.  (Eventually should delete and clean up.)
      if (flag.eq.9999.0) then  ! start step 3 if-construct
        CALL GETNEWT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &    dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &    ro_snow,Cp_snow,T_tmp,Tf,xk_water,water_frac,
     &    up,down,total_solar,xLs,Rv,ro_water,ro_pure_ice,
     &    xinternal_heating_corr,ndarklayers,Sc)

	call CONDUCT(icondflag,Qc3,gamma,T_old,dy_p,JJ,Tsfc)

c mgc totalheat3
	if(ifinalcall.eq.1) then
	totalheat3=0.0
	sourceterm=0.0
	do mmm=1,JJ
	totalheat3=totalheat3+
     &	dble((T_old(mmm)-270.0)*Cp_snow*dy_p(mmm)*ro_snow)
	totalheat3=totalheat3+dble(water_frac(mmm)*dy_p(mmm)*ro_snow*xLf)
	if (Sc(mmm).lt.1e6) sourceterm=sourceterm+Sc(mmm)*dy_p(mmm)
	enddo	
c	print *,'day,hr',iter.hr,net heat,internal,internal change,Qc,Qsip
c	print '(i4, f12.1, f15.1, 3f12.1)',3,
c     &	totalheat3-totalheat0 + dble(Qc3*dt -
c     &	SOURCETERM*dt),
c     &     totalheat3,
c     &	totalheat3-totalheat0,Qc3*dt,
c     &	SOURCETERM*dt
c	write (83,*)  totalheat3-totalheat0,dble(Qc3)*dble(dt),
c     &	dble(SOURCETERM)*dble(dt)

	endif

	else
c write dummy for q3
c	if (ifinalcall.eq.1) write (83,*)  0.0,0.0,0.0
      endif  !end step 3 if-construct

c Compute a total-column water depth.
      water_depth = 0.0
      do j=1,JJ
        water_depth = water_depth + dy_p(j) * water_frac(j)
      enddo
      
c mgc checks if water_depth is >0
      
	if (water_depth .gt. 0.0) then
	hi =1
	endif

	if (water_depth .gt. 1.0) then
	write(*,*) 'water depth is greater than 1.0! ', water_depth
    !     print *, water_frac
    !     stop
	endif

c Compute the water flux by looking at the difference between the
c   current water depth and that at the previous time step.
      water_flux = water_depth - water_depth_old
      water_depth_old = water_depth

c Choose which conducton value to return
c	if (flag.eq.1.0) then
c		Qcout=Qc3
c	else
		Qcout=Qc2
c	endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONSTS_ICE(xLs,xLf,Tf,ro_water,Cp_water,xk_water,
     &  xk_ice,ro_ice,Cp_snow,Rv,ro_snow,xk_snow,ro_pure_ice)

c Jon has these 2 for some reason I dont quite understand
c	ro_ice = 900.0
	ro_snow = ro_ice
	
	ro_pure_ice = 917.0
	xLv = 2.500e6
c      xLs = 2.500e6
      xLf = 3.34e5
	xLs = xLv + xLf
      Tf = 273.16
      ro_water = 1000.0
      Cp_water = 4180.0
      xk_water = 0.552
c      xk_ice = 2.10
	xk_ice = 1.8
c      ro_ice = 917.0
      Cp_snow = 2106.0
      Rv = 461.0

c Compute the thermal conductivity from the snow density.
      if (ro_snow.lt.156.0) then
        xk_snow = 0.023 + 0.234 * (ro_snow/1000.0)
      elseif ((ro_snow .lt. 600) .and. (ro_snow .lt. 156)) then
        xk_snow = 0.138 - 1.01 * (ro_snow/1000.0) + 3.233 *
     &    (ro_snow/1000.0)**2
	else 
		xk_snow = xk_ice
      endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ICEMF(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice,
     &  ro_snow,water_frac,flag,ro_pure_ice)

      real dy_p(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real freeze(JJ+2)
      real ice_avail(JJ+2)
      real water_frac(JJ+2)

c I don't think flag is used anymore...
c Compute the maximum ice available to be melted.
      do j=1,JJ
c        ice_avail(j) = ro_snow / ro_pure_ice * dy_p(j)
          ice_avail(j) = dy_p(j)
      enddo

      flag = 0.0
      extramelt = 0.0
      do j=1,JJ
        if (T_old(j).ge.Tf  .or.  extramelt.gt.0.0) then !==============

		extrameltprevious=extramelt ! for reference

c Compute the amount of water produced by the energy represented
c   by the temperature above freezing.
          xmelt(j) = Cp_snow * dy_p(j) * (T_old(j) - Tf) / xLf
          totalmelt = xmelt(j) + extramelt
	
c in some cases (if T is below freezing but extramelt is available from the layer above)
c totalmelt will be a negative number because xmelt is negative.  This is bad.  
c In that case use the 'melt' energy to raise the temperature of this cell instead. 

		if (totalmelt .gt. 0.0) then
c Add the new water to the old water.
			water_frac(j) = water_frac(j) + totalmelt / ice_avail(j)

c mgc everthing below here is essentially the same as liston, but hoffman
c added some checks as noted above, but no density or drainage stuff
c Assume that energy producing a water fraction above 1.0 goes 
c   into melting the ice below (as 'extramelt').
			extramelt = max(0.0,water_frac(j) - 1.0) * ice_avail(j)
c	extramelt = 0.0 ! assume water drains
			water_frac(j) = min(1.0,water_frac(j))
c Because ice is melting, the temperature must be Tf.
          		T_old(j) = Tf
          		flag = 1.0
		else ! This is the rare strange case - the change should be miniscule
c  this happens if there was extramelt, but the temp of this cell is far enough <0
c to use up all the extramelt and still have T<0
			T_old(j) = T_old(j) + extramelt * xLf / (Cp_snow*dy_p(j))
			extramelt=0.0 !all used up!
		endif

	if (water_frac(j).lt.0.0) then
	write (*,*) 'WARNING: Water fraction after adding water
     &	     at level # is neg: ', 	j,water_frac(j)
	endif

        else  !=====================

c Three cases are possible, 1) there is no water to freeze, so do
c   nothing, 2) all of the water freezes and the extra energy drops
c   the  energy below Tf, or 3) only some of the water freezes and
c   the temperature remains at Tf.
	if (water_frac(j).lt.0.0) then
	write (*,*) 'WARNING: Water fraction at level # is neg: ', 
     &	j,water_frac(j)
	endif

          if (water_frac(j).gt.0.0) then

c Compute the amount of water frozen by the energy represented
c   by the temperature below freezing.
          freeze(j) = Cp_snow * dy_p(j) * (Tf - T_old(j)) / xLf

            if (freeze(j).le.water_frac(j)*ice_avail(j)) then
c Case 3.
              water_frac(j) = water_frac(j) - freeze(j) / ice_avail(j)
              T_old(j) = Tf
              flag = 1.0

	if (water_frac(j).lt.0.0) then
	write (*,*) 'WARNING: Water frac after Case3 at level # is neg:', 
     &	j,water_frac(j)
	endif
            else
c Case 2.
c              freeze(j) = water_frac(j) * ice_avail(j)

              freeze(j) = freeze(j) - water_frac(j) * ice_avail(j)
              water_frac(j) = 0.0
              T_old(j) = Tf - freeze(j) * xLf / (Cp_snow * dy_p(j))
              flag = 1.0

	if (water_frac(j).lt.0.0) then
	write (*,*) 
     &	'WARNING: Water frac after Case2 at level # is neg: ', 
     &	j,water_frac(j)
	endif
            endif
          endif

        endif ! if t_old>0 or extramelt>0    =========================

	if (water_frac(j).lt.0.0) then
	write (*,*) 'WARNING: Water fraction at end at level # is neg: ', 
     &	j,water_frac(j)
	endif

      enddo !z-layer loop

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc not sure what role this plays, its a bit shorter than ICEMF but
c is not used afaik

      SUBROUTINE ICEMF2(T_old,JJ,dy_p,xmelt,Cp_snow,xLf,Tf,ro_ice,
     &  ro_snow,water_frac,flag,ro_pure_ice)

      real dy_p(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real freeze(JJ+2)
      real ice_avail(JJ+2)
      real water_frac(JJ+2)


c Compute the maximum ice available to be melted.
      do j=1,JJ
c        ice_avail(j) = ro_snow / ro_ice * dy_p(j)
        ice_avail(j) = dy_p(j) !assuming density never changes
      enddo

      flag = 0.0
      extramelt = 0.0
      do j=JJ,1,-1
        if (T_old(j).ge.Tf  .or.  extramelt.gt.0.0) then

c Compute the amount of water produced by the energy represented
c   by the temperature above freezing.
          xmelt(j) = Cp_snow * dy_p(j) * (T_old(j) - Tf) / xLf
          totalmelt = xmelt(j) + extramelt

c Add the new water to the old water.
          water_frac(j) = water_frac(j) + totalmelt / ice_avail(j)

c Assume that energy producing a water fraction above 1.0 goes 
c   into melting the ice below (as 'extramelt').
          extramelt = max(0.0,water_frac(j) - 1.0) * ice_avail(j)
          water_frac(j) = min(1.0,water_frac(j))

c Because ice is melting, the temperature must be Tf.
          T_old(j) = Tf
          flag = 1.0

        else

c Three cases are possible, 1) there is no water to freeze, so do
c   nothing, 2) all of the water freezes and the extra energy drops
c   the  energy below Tf, or 3) only some of the water freezes and
c   the temperature remains at Tf.
	if (water_frac(j).lt.0.0) then
	write (*,*) 'Water fraction at level # is neg: ', j,water_frac(j)
	endif

          if (water_frac(j).gt.0.0) then

c Compute the amount of water frozen by the energy represented
c   by the temperature below freezing.
          freeze(j) = Cp_snow * dy_p(j) * (Tf - T_old(j)) / xLf

            if (freeze(j).le.water_frac(j)*ice_avail(j)) then
c Case 3.
              water_frac(j) = water_frac(j) - freeze(j) / ice_avail(j)
              T_old(j) = Tf
              flag = 1.0

	if (water_frac(j).lt.0.0) then
	write (*,*) 'Froze the pooch!Waterfrac# is neg: '
     &	, j,water_frac(j)
	water_frac(j) = 0.0
	endif

            else
c Case 2.
c              freeze(j) = water_frac(j) * ice_avail(j)

              freeze(j) = freeze(j) - water_frac(j) * ice_avail(j)
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
     &  up,down,total_solar,xLs,Rv,ro_water,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,Sc)

c mgc everything folded is identical to liston
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


	if (Qsi.gt.100) then
	hi=1
	endif


c Compute the solar radiation penetrating the lake surface.
      CALL SOLARPEN(Qsip,Qsi,albedo)

      xTsfc = Tsfc

c Compute gamma and the general equation coefficients.
      CALL GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)
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

c mgc internal heating correction

          Sc(j) = - Qsip / (total_solar * (1-0.562)) *
     &      (xynet(j) - xynet(j+1)) / dy_p(j)
	Sc(j)=xinternal_heating_corr * Sc(j)

c          Sc(j) = - Qsi / total_solar *
c     &      (xynet(j) - xynet(j+1)) / dy_p(j)

c mgc skin heating

c Turn off heating in the top layers that contribute to SEB directly
	if (j .le. ndarklayers) then
		Sc(j) = 0.0
	end if

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
c mgc identical to liston
      SUBROUTINE SOLARPEN(Qsip,Qsi,albedo)
	implicit none
	real Qsip,Qsi,albedo

      Qsip = (1.0 - albedo) * Qsi

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc identical to liston
      SUBROUTINE ICEINIT(Tsfc,T_old,dely_p,f_n,y_crds,y_wall,dy_p,JJ,
     &  Tf,water_frac,gamma,xk_snow,water_depth_old,temp_ice_init_C,
     &  deltaz)

      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real y_wall(JJ+2)
      real T_old(JJ+2)
      real gamma(JJ+2)
      real water_frac(JJ+2)
      real deltaz(JJ)

c Provide values of Control Volume size in the y direction, and
c   compute c.v. size information.
      CALL GETCV(deltaz,dy_p,JJ)
      CALL CV_INFO(dely_p,f_n,y_crds,y_wall,dy_p,JJ)

c Supply the initial conditions.
      do j=1,JJ
        T_old(j) = temp_ice_init_C + Tf
        water_frac(j) = 0.0
        gamma(j) = xk_snow
      end do

      water_depth_old = 0.0

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE ICEHEAT(gamma,T_old,Tsfc,JJ,dy_p,y_crds,
     &  dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &  ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &  xLs,Rv,Tf,ro_water,ro_pure_ice,xinternal_heating_corr,
     &  ndarklayers)

c mgc folded is identical to liston
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

	if (Qsi.gt.100) then
	hi=1
	endif
      xTsfc = Tsfc

c Compute gamma and the general equation coefficients.
      CALL GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)
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
	
c	bc_S = 0.0
c	aS(1) = 0.0
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


c mgc see notes below about different version of dQdz equation

c This is dq/dz for q=Q0*exp(-extcoef*z).
c       Sc(j) = Qsip * extcoef * exp(- extcoef * y_crds(j+1))

c This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
c   are scaled by the ratio of the incoming solar to the solar
c   used to get up and down.
c After further review, mh change this to first version (see J 1119)
c        Sc(j) = - (xynet(j) - xynet(j+1)) / dy_p(j)
ccc        Sc(j) = - Qsi / total_solar * (xynet(j) - xynet(j+1)) / dy_p(j)
	Sc(j) = - Qsip / (total_solar *(1-0.562))
     &	 * (xynet(j) - xynet(j+1)) / dy_p(j)
	Sc(j)=xinternal_heating_corr * Sc(j)

c Eliminate internal heating in the upper few layers and let that heat be part of the 
c surface energy balance instead.  
c The Qsi is adjusted by this amount in the TSFC subroutine.
	if (j .le. ndarklayers) then
		Sc(j) = 0.0
	end if

c output internal heating
c	write (67,*) j,Sc(j)*dy_p(j),Qsi,albedo
           
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
      real deltaz(JJ)

c Provide values of Control Volume size in the y direction.
      do j=1,JJ
        dy_p(j) = deltaz(j)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETGAMMA1(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice)

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
c        air_frac = 1.0 - ro_snow/ro_water
        air_frac = 1.0 - ro_snow/ro_pure_ice
        gamma(j) = (1.0-air_frac) *
     &    ((1.0-water_frac(j)) * xk_snow + water_frac(j) * xk_water) +
     &    air_frac * xk_vapor

c		gamma(j) = gamma(j)*4.0

c       gamma(j) = (1.0 - water_frac(j)) * xk_snow +
c    &    water_frac(j) * xk_water + xk_vapor
      end do

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc hoffman added new methods to compute thermal K 

      SUBROUTINE GETGAMMA(gamma,JJ,xk_snow,xk_water,water_frac,
     &  T_old,xLs,Rv,Tf,ro_snow,ro_water,ro_pure_ice,albedo)

      real gamma(JJ)
      real water_frac(JJ)
      real T_old(JJ)

c MH: The original GETGAMMA subroutine generates conductivity values
c that seem very low, resulting in a subsurface temp. increases being very high.
c Since that subroutine is based on Sturm's research on snow, and we only have ice
c I have decided to use empirical equations based on ice.  These come from Paterson, p.205.
c After further review, the low GETGAMMA values may be due to a low value for ice (1.8),
c and not the routine itself.  In any case, this is an alternate method.

      do j=1,JJ
c Determine K for pure ice at the given temperature.
		gammapureice=9.828*exp(-5.7e-3*T_old(j))

c Determine K for ice at given density at that temperature (Schwerdtfeger, 1963 - upper limit)
		gammaice1 = 2*gammapureice*ro_snow / (3*ro_pure_ice - ro_snow)

c  van Dusen 1929 formula (lower limit)
	gammaice2=(2.1e-2) + (4.2e-4) * ro_snow + (2.2e-9) * ro_snow**3

c Choose one or the other or the average
		gammaice=(gammaice1+gammaice2)/2.0
		
c Assume no air fraction, but account for water fraction
	
      gamma(j) = (1.0-water_frac(j))*gammaice + water_frac(j)*xk_water

      end do

	if (albedo.gt.0.70) then  !simulate thermal covering of snow
!0.73 is value needed to achieve diffusivity corresponding to cond=0.3 but with ice density
c		gamma(1)=0.6
c		gamma(2)=0.6
	endif

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

      SUBROUTINE TRISOLVE(x,asub,amain,asuper,b,JJ)

	implicit none

      integer JJ,j
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

c mgc appears essentially same as liston, but with a bunch of additional
c variables declared and called by ENBALANCE, SFCTEMP, and MFENERGY,
c many of which relate to the subsurface, so could be clues into the coupling

      SUBROUTINE ENBALANCE(Tair,windspd,rh,
     &    Tsfc,Qsi,Qli,Qle,Qh,
     &    Qe,Qc,Qm,balance,Qf,
     &    swe_depth,topo,z_windobs,dt,gamma,
     &    T_old,dy_p,JJ,icond_flag,cloud_frac,albedo,z_0,
     &    J_day_start,xlat,slope_az,terrain_slope,
     &    transmiss,clear_sky,
     &	snow_cover_depth_old,surface_melt,ro_snow_on_top,
     &	ablation,model_day,xLs,xLf,i_yearstart,
     &	ablation_output,stability,hr,Qsi_fraction,
     &y_crds,f_n,dely_p,Qsip,extcoef,xk_snow,ro_snow,
     &Cp_snow,xk_water,water_frac,up,down,total_solar,
     &Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,ro_pure_ice,kkk,
     &	xinternal_heating_corr,qsfactor,ndarklayers,Pa)

	integer ablation_output,model_day,hr,ndarklayers

c These are here to handle the case of non-zero conduction.
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

c Define the constants used in the computations.
        CALL CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity,
     &    xkappa,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice,
     &    ihrs_day)

c Atmospheric vapor pressure from relative humidity data.
        CALL VAPPRESS(ea,rh,Tair,Tf)

c Jon calcs ro_air here.  Constants could be moved to CONSTS, but this is easier
		real_Ma = 32.56e-3
		R = 8.314
		epsilon = 0.62201
		ro_air = Pa * real_Ma / (R * Tair) * (1.0 + (epsilon - 1.0) 
     &		* (ea/Pa))

c Compute the incoming longwave radiation.
        CALL LONGIN(Qli,ea,Tair,Stef_Boltz)

c Compute the incoming solar radiation.
c        CALL SOLARIN(J_day_start,model_day,dt,Qsi,xlat,cloud_frac,
c     &    slope_az,terrain_slope,ihrs_day,transmiss,
c     &	clear_sky,Pa,one_atmos,i_yearstart,Qsi_fraction,hr)

c Compute the turbulent exchange coefficients. AND z_0
c        CALL EXCOEFS(D_h,D_e,z_0,z_windobs,windspd,xkappa,model_day)
      CALL EXCOEFS_SCALAR(D_h,D_e,z_0,z_windobs,windspd,xkappa,
     &	model_day)

c Compute the flux contribution due to conduction.
c Move this into SFCTEMP
c        CALL CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ)

c Solve the energy balance for the surface temperature.
        CALL SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,
     &    D_h,D_e,Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,
     &    Stef_Boltz,gravity,xLs,xkappa,z_0,Tf,Qc,model_day,
     &	icond_flag,gamma,T_old,dy_p,JJ,hr,
     &y_crds,dt,f_n,dely_p,Qsip,extcoef,xk_snow,
     &ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,xLf,ro_pure_ice,kkk,
     &	xinternal_heating_corr,qsfactor,ndarklayers)

c Make sure the snow surface temperature is <= 0 C.
        CALL MELTTEMP(Tsfc,Tf,swe_depth)

c Compute the stability function.
        CALL STABLEFN(stability,Tair,Tsfc,windspd,z_windobs,
     &    gravity,xkappa,z_0)
	
c Compute the water vapor pressure at the surface.
        CALL VAPOR(es0,Tsfc,Tf)

c Compute the latent heat flux.
        CALL LATENT(Qe,D_e,stability,ea,es0,ro_air,
     &    xLs,Pa)

c Compute the sensible heat flux.
        CALL SENSIBLE(Qh,D_h,stability,Tair,Tsfc,
     &    ro_air,Cp)

c Compute the longwave flux emitted by the surface.
        CALL LONGOUT(Qle,Tsfc,emiss_sfc,Stef_Boltz)

c Get some additional info about snow events (Jon)
c	CALL SNOWEVENT(model_day,
c     &  snow_cover_depth_old,i_yearstart,ro_snow_on_top,ro_water,
c     &  Qsi,Qsi_fraction,albedo,rh,windspd)

c Check for snow on surface
c	if (albedo.ge.0.65) then
c	snow_cover_depth_old = 1.0
c	else
c	snow_cover_depth_old = 0.0
c	endif
	
c Compute the energy flux available for melting or freezing.
        CALL MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,
     &    Qc,Qm,Qf,Tsfc,Tf,Tair,windspd,z_windobs,
     &    gravity,De_h,ea,ro_air,xLs,Pa,Cp,emiss_sfc,
     &    Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p,
     &    JJ,icond_flag,
     &    xLf,ro_water,surface_melt,ro_snow_on_top,ablation,
     &    snow_cover_depth_old,model_day,ablation_output,hr,qsfactor,dt)

c Decrease the swe depth by the swe melt depth.
c   Turn this off for blue-ice simulations.
c       CALL SNOW_UPDATE(swe_depth,Qm,dt,ro_ice,xLf)

c Perform an energy balance check.
        CALL ENBAL(balance,albedo,Qsi,Qli,Qle,Qh,
     &    Qe,Qc,Qm,qsfactor)
	if (balance.gt.2.0) then
	print *,'big balance! ', balance
	endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LONGIN(Qli,ea,Ta,Stef_Boltz)

c Read from file
c	read (37,*) Qli

c Check for the few bad days that exist.  If no reading available (-999)
c then use Glen's formulation to estimate Qli

	if (Qli .lt. 0.0) then

	print *,'Qli data missing'
	stop
c Compute Qli.
      emiss_cloud = 1.08 * (1.0 - exp(-(0.01 * ea)**(Ta/2016.)))
      Qli = emiss_cloud * Stef_Boltz * Ta**4

	endif


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc hoffman added a Qsi_fraction that corrects the computed Qsi based
c on station data

      SUBROUTINE SOLARIN(J_day_start,iter,dt,Qsi,xlat,cloud_frac,
     &  slope_az,terrain_slope,ihrs_day,transmiss,
     &  clear_sky,Pa,one_atmos,i_yearstart,Qsi_fraction,hr)

c Compute the incoming solar radiation.  Here I am going to assume
c   that we are using daily time steps.  If you have some other time
c   step, see MicroMet for ideas about what to do.
      J_day = iter + J_day_start - 1

c      Qsi_sum = 0.0
      if (dt.eq.3600.0) then
c        do ihour=1,ihrs_day
          xhour = real(hr)
          call SOLAR_RAD(Qsi_tmp,J_day,xlat,cloud_frac,
     &      xhour,slope_az,terrain_slope,transmiss,
     &	  clear_sky,Pa,one_atmos,i_yearstart)
c            Qsi_sum = Qsi_sum + Qsi_tmp
c        enddo
c Qsi here is Jon's Qsi_calculated
c        Qsi = Qsi_sum / real(ihrs_day)
		Qsi=Qsi_tmp
      else
        print *,'Need to fix the solar routines for this dt'
        stop
      endif

c Read the Qsi_fraction at the nearest lake station from file
c	read (32,*) Qsi_fraction
c Now transform Qsi from calculated to actual
	Qsi = Qsi * Qsi_fraction

c Note if we want Qsi_calc_daily output later, we will need to save the value
	
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc some slope/azimuth differences but mostly the same

      SUBROUTINE SOLAR_RAD(Qsi,J_day,xlat,cloud_frac,
     &  xhour,slope_az,terrain_slope,transmiss,
     &  clear_sky,Pa,one_atmos,i_yearstart)

      implicit none

      integer J_day,i_yearstart

      real solar_const,days_yr,Trop_Cap,solstice,pi,deg2rad,
     &  cos_i,cos_Z,Qsi,xlat,sin_z,xhour,
     &  cloud_frac,slope_az,terrain_slope,sol_dec,hr_angl,
     &  trans_direct,trans_diffuse,Qsi_trans_dir,Qsi_trans_dif,
     &  sun_azimuth,slope_az_S0,transmiss
c	real Qsi_diffuse,Qsi_direct

c additional vars that Jon uses
	real GGG,eccentricity,cos_AZ,sun_azimuth_deg,tmp
	real clear_sky,Pa,one_atmos

c Required constants.
      solar_const = 1370.
      days_yr = 365.25
      Trop_Cap = 0.41
c      solstice = 173.
      pi = 2.0 * acos(0.0)
      deg2rad = pi / 180.0

	if (MOD(i_yearstart,4) .eq. 0) then
	  solstice = 174
	else
	  solstice = 173
	endif

c Jon has these as well
      if (J_day.gt.365) then
		J_day = J_day - 365
      endif
	GGG = 2.0 * pi * (J_day - 1)/days_yr
	eccentricity = 1.00011 + 0.034221 * cos(GGG) + 0.00128 * sin(GGG)
     &	+ 0.000719 * cos(2*GGG) + 0.000077 * sin(2*GGG)


c COMPUTE THE BASIC SOLAR RADIATION PARAMETERS.

c Compute the solar declination angle (radians).
      sol_dec = Trop_Cap *
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
c      sun_azimuth = 
c     &  asin(max(-1.0,min(1.0,cos(sol_dec)*sin(hr_angl)/sin_Z)))
c Solar azimuth using equations from (Ebnet 2005) where north has zero
c azimuth.
	cos_AZ = (sin(sol_dec) * cos(xlat * deg2rad) - cos(sol_dec) * 
     &	sin(xlat * deg2rad) * cos(hr_angl))/sin_Z
	if ((xhour.ge.1).and.(xhour .le. 11)) then
	  sun_azimuth = acos(cos_AZ)
	elseif (xhour .eq. 12) then
	  sun_azimuth = 0
	elseif ((xhour .ge. 13).and.(xhour .le. 23)) then
	  sun_azimuth = 2 * pi - acos(cos_AZ)     
	elseif (xhour .eq. 24) then
	  sun_azimuth = pi
	endif
	sun_azimuth_deg = sun_azimuth * 180/pi

cc I believe the previous statements are meant to replace this bit
cc Make the corrections so that the angles below the local horizon
cc   are still measured from the normal to the slope.
c      if (hr_angl.lt.0.0) then
c        if (hr_angl.lt.sun_azimuth) sun_azimuth = - pi - sun_azimuth
c      elseif (hr_angl.gt.0.0) then
c        if (hr_angl.gt.sun_azimuth) sun_azimuth = pi - sun_azimuth
c      endif
c
cc Build, from the variable with north having zero azimuth, a 
cc   slope_azimuth value with south having zero azimuth.
c      if (slope_az.ge.180.0) then
c        slope_az_S0 = slope_az - 180.0
c      else
c        slope_az_S0 = slope_az + 180.0
c      endif
	slope_az_s0 = slope_az

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
c Jon added these:
      if (cos_Z.le.0.0) then
      	tmp = 0.0
      else 
      	tmp = solar_const * (eccentricity**2) * (clear_sky**(Pa / 
     &    (one_atmos * cos_Z))) * cos_i
      endif
      
c Jon skips the whole direct/diffuse thing!
cc Adjust the solar radiation for slope, etc.
c      Qsi_direct = cos_i * Qsi_trans_dir
c      Qsi_diffuse = cos_Z * Qsi_trans_dif
c
cc Combine the direct and diffuse solar components.
c      Qsi = Qsi_direct + Qsi_diffuse
c
	Qsi = tmp
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc hoffman reads in station pressure - wonder if this makes any difference 

      SUBROUTINE PRESSURE(Pa,one_atmos,scale_ht,topo,ro_air,Tair,
     &	Tf,gravity)

c Compute the average station pressure.
c Glen's Pa calc:
c      Pa = one_atmos * exp(- topo / scale_ht)

	read (36,*) Pa_ref,T_ref,elev_ref

	Pa=(Pa_ref*100.0) * exp( (topo-elev_ref) / 
     &	(-287.04*0.5*(T_ref+Tf + Tair)/gravity) )


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc no difference
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

cc Over water.
c        A = 6.1121 * 100.0
c        B = 17.502
c        C = 240.97
c Over ice.
       A = 6.1115 * 100.0
       B = 22.452
       C = 272.55

c Atmospheric vapor pressure from relative humidity data.
      ea = rh / 100.0 * A * exp((B * (Tair - Tf))/(C + (Tair - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MFENERGY(albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc,
     &  Tf,Tair,windspd,z_windobs,gravity,De_h,ea,ro_air,xLs,Pa,Cp,
     &  emiss_sfc,Stef_Boltz,swe_depth,xkappa,z_0,gamma,T_old,dy_p,
     &  JJ,icond_flag,
     &  xLf,ro_water,surface_melt,ro_snow_on_top,ablation,
     &	snow_cover_depth_old,model_day,ablation_output,hr,qsfactor,dt)

	implicit none

	real albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,Qf,Tsfc
	real Tf,Tair,windspd,z_windobs,gravity,De_h,ea,ro_air,xLs,Pa,Cp
	real emiss_sfc,Stef_Boltz,swe_depth,xkappa,z_0
	integer JJ,model_day,icond_flag,hr
      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)
	real xLf,ro_water,total_surface_melt,ro_snow_on_top
	real snow_cover_depth_old
	integer ablation_output
	
!	real albedo_evo_poly_a,albedo_evo_poly_b,albedo_evo_poly_c
!	real albedo_evo_poly_d,calculated_albedo
	real surface_melt,ablation
	real qsfactor,dt

c	real albedo_tom,calculated_albedo_tom,snow_cover_depth,ablation_tmp


c constants for calculating the best fit polynomial of albedo.
c See the start of the TIMESTEP loop for usage
!	albedo_evo_poly_a = 0.8011
!	albedo_evo_poly_b = -5.353e-3
!	albedo_evo_poly_c = 2.762e-5
!	albedo_evo_poly_d = -3.775e-8
c This equation is used to determine the best fit third order
c polynomial of albedo during the summer months (October 1st through March
c 31st) for an ice-exposed surface.  The regression line is used in the
c model for albedo of the terminal cliffs and when no snow exists on the
c surface.  
!		calculated_albedo=albedo_evo_poly_a + albedo_evo_poly_b * 
!     &		model_day + albedo_evo_poly_c * model_day**2 + 
!     &		albedo_evo_poly_d * model_day**3

c If Qm is > 0, then this is the energy available for melting.
c   If Qm is < 0, then this is the energy available for freezing
c   liquid water in the snowpack.
      if (swe_depth.gt.0.0 .and. Tsfc.eq.Tf) then
c MH: added qsfactor here to account for energy absorbed below surface
        Qm = (1.0-albedo) * Qsi * qsfactor + Qli + Qle + Qh + Qe + Qc
      else
        Qm = 0.0
      endif

c Jon seems to have rewritten the rest of this subroutine
c so i am duplicating it below
c      if (Tsfc.lt.Tf) then
c        xTsfc = Tf
c        CALL STABLEFN(xstability,Tair,xTsfc,windspd,z_windobs,
c     &    gravity,xkappa,z_0)
c        CALL VAPOR(xes0,xTsfc,Tf)
c        CALL LATENT(xQe,De_h,xstability,ea,xes0,ro_air,xLs,Pa)
c        CALL SENSIBLE(xQh,De_h,xstability,Tair,xTsfc,ro_air,Cp)
c        CALL LONGOUT(xQle,xTsfc,emiss_sfc,Stef_Boltz)
c        CALL CONDUCT(icond_flag,xQc,gamma,T_old,dy_p,JJ)
c        Qf = (1.0-albedo) * Qsi + Qli + xQle + xQh + xQe + xQc
c      else
c        Qf = 0.0
c      endif

c Jon has extra bit here:
c Calculate the ablation and melt of the ice surface.
c Surface melt...

c mgc ABLATION	  
	  
	if (snow_cover_depth_old.gt.0) then
		surface_melt = 0.0
	else
		surface_melt = Qm / (xLf * ro_water) * dt * 100.0
	endif
	total_surface_melt = total_surface_melt + surface_melt

c Check Values
	if (surface_melt.lt.0.0) then
		if (surface_melt .gt. -1e-4) then
			surface_melt = 0
		else
			print *,'WARNING: melt<0!', surface_melt
		endif
	endif

	if (surface_melt.gt.1.0) then
		print *,'WARNING: bigmelt!', surface_melt
	endif

c Ablation...(no melt can occur if snow is on top)
	if (snow_cover_depth_old.gt.0.0) then
c		ablation = -Qe/(xLs * ro_snow_on_top) * dt * 100.0
c Don't count ice sublimation if snow is on top		
		ablation = 0.0
	else
		ablation = surface_melt - Qe/(xLs * ro_water) 
     &		* dt * 100.0
	endif
	if (ablation.lt. 0.0) then
c		ablation = 0.0
	endif

c	if (snow_cover_depth_old.gt.0.0) then
c		ablation_tmp = 0.0
c	else
c		ablation_tmp = ablation
c	endif

c Now let us compute the new total snow depth.  
c	snow_cover_depth = 0.0
c	snow_cover_depth = snow_cover_depth_old - ablation
c	if (snow_cover_depth.lt.0.0) then
c		snow_cover_depth = 0.0
c	endif
c	if ((model_day.ge.115) .and. (model_day.le.215)) then
c Check tomorrow's albedo
c		read (33,*) albedo_tom
c		BACKSPACE (33)
c	calculated_albedo_tom=albedo_evo_poly_a + albedo_evo_poly_b * 
c    &		(model_day+1) + albedo_evo_poly_c * (model_day+1)**2 + 
c   &		albedo_evo_poly_d * (model_day+1)**3
c
c		if ((albedo_tom - calculated_albedo_tom .ge. 0.1) .and.
c     &		(snow_cover_depth_old .ne. 0.0)) then
c			snow_cover_depth_old = max(1.0,snow_cover_depth)
c		else
c			snow_cover_depth_old = snow_cover_depth
c		endif
c	else
c		snow_cover_depth_old = snow_cover_depth
c	endif


      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE MELTTEMP(Tsfc,Tf,swe_depth)

      if (swe_depth.gt.0.0 .and. Tsfc.gt.Tf) then
	Tsfc = Tf
	endif

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc see notes about using sfc and level 1 vs level 1 and 2 for ghf
      SUBROUTINE CONDUCT(icond_flag,Qc,gamma,T_old,dy_p,JJ,Tsfc)

      real gamma(JJ+2)
      real T_old(JJ+2)
      real dy_p(JJ+2)

      if (icond_flag.eq.0) then
        Qc = 0.0
      else
c This was MH attempt to make conduction based on Levels Sfc & 1,
c rather than 1 &2.  This greatly smoothed the Sfc Temp curve, but
c smoothed it too much relative to IRT data.
c		if (Tsfc.le.273) then
			Qc = - gamma(1) * (Tsfc-T_old(1)) / (dy_p(1)*0.5)

c	Qc = - (gamma(1)*dy_p(1)+gamma(2)*dy_p(2)*0.5)/
c     &	(dy_p(1)+dy_p(2)*0.5)*(Tsfc-T_old(2)) / (dy_p(1)+dy_p(2)*0.5)

c		else
c replace 1&2 with 3&4
c	        Qc = - (gamma(1) + gamma(2))/2.0 * (T_old(1) - T_old(2)) /
c     &		    ((dy_p(1) + dy_p(2))/2.0)
c		endif
      endif

c NOTE: The ice heat flux is calculated using the 1st and 2nd levels within
c the profile of the ice instead of the ice at the surface and the 1st
c level.  This allows for surface ice at the melting point to still have a
c positive ice heat flux where as surface ice at the melting point should
c have zero gain or loss of energy to the ice below.

      return
      end
      

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE CONSTS(emiss_sfc,Stef_Boltz,ro_air,Cp,gravity,
     &  xkappa,Tf,ro_water,one_atmos,scale_ht,Cp_water,ro_ice,
     &  ihrs_day)

      emiss_sfc = 0.98
      Stef_Boltz = 5.6696e-8
      ro_air = 1.275
      Cp = 1004.
      gravity = 9.81
c      xLs = 2.500e6
      xkappa = 0.4
c     xLf = 3.34e5
      Tf = 273.16
      ro_water = 1000.0
      one_atmos = 101300.0
      scale_ht = 8500.0
      Cp_water = 4180.0
      ro_ice = 917.0
      ihrs_day = 24
c Note Ma, epsilon, R added to PRESSURE - easier than putting them here
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc Chi parameter

      SUBROUTINE ENBAL(balance,albedo,Qsi,Qli,Qle,Qh,Qe,Qc,Qm,qsfactor)

      balance = qsfactor * (1.0-albedo) * Qsi + 
     &	Qli + Qle + Qh + Qe + Qc - Qm 

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

      SUBROUTINE SENSIBLE(Qh,D_h,stability,Tair,Tsfc,ro_air,Cp)

      Qh = ro_air * Cp * D_h * stability * (Tair - Tsfc)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc hoffman added qsfactor and used Brent's method instead of SOLVE

      SUBROUTINE SFCTEMP(Tsfc,Tair,Qsi,Qli,ea,albedo,D_h,D_e,
     &  Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity,
     &  xLs,xkappa,z_0,Tf,Qc,model_day,
     &  icond_flag,gamma,T_old,dy_p,JJ,hr,
     &y_crds,dt,f_n,dely_p,Qsip,extcoef,xk_snow,
     &ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,xLf,ro_pure_ice,kkk,
     &	xinternal_heating_corr,qsfactor,ndarklayers)
	
	integer JJ

      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

	integer hr, ndarklayers

c Solve the energy balance for the surface temperature.

c Optional output of details of Tsfc solution 
c	if (model_day.eq.2) then
c		CALL EBPLOT(Tsfc,Tair,Qsi,Qli,ea,albedo,De_h,
c     &  Pa,z_windobs,windspd,ro_air,Cp,emiss_sfc,Stef_Boltz,gravity,
c     &  xLs,xkappa,z_0,Tf,Qc,model_day,icond_flag,gamma,T_old,dy_p,JJ)
c	endif

      AAA = ro_air * Cp * D_h
      CCC = 0.622 / Pa
      DDD = emiss_sfc * Stef_Boltz
c      EEE = (1.0-albedo) * Qsi + Qli + Qc
cc	EEE = (1.0-albedo) * Qsi + Qli
	EEE = qsfactor * (1.0-albedo) * Qsi + Qli
      FFF = ro_air * xLs * D_e

c Compute the constants used in the stability coefficient
c   computations.
      C1 = 5.3 * 9.4 * (xkappa/(log(z_windobs/z_0)))**2 *
     &  sqrt(z_windobs/z_0)
      C2 = gravity * z_windobs / (Tair * windspd**2)
      B1 = 9.4 * C2
      B2 = C1 * sqrt(C2)

      CALL SOLVE_BRENT
     &	(Tsfc,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf,model_day,
     &	icond_flag,gamma,T_old,dy_p,JJ,hr,
     &y_crds,dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &xLs,Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,xLf,ro_pure_ice,kkk,
     &	xinternal_heating_corr,ndarklayers,Qc)
	 
c mgc Tsfc skin heating	 

c !!!remove the heat from the top layer that got put into the SEB
c	T_old(1)=T_old(1) - Qc/ (dy_p(1) * 2106.0)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVE_BRENT(xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,
     &	Tf,model_day,icond_flag,gamma,T_old,dy_p,JJ,hr,
     &y_crds,dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &xLs,Rv,ro_water,xmelt,ro_ice,water_depth,	
     &water_depth_old,water_flux,xLf,ro_pure_ice,kkk,
     &	xinternal_heating_corr,ndarklayers,Qc)
	implicit none

	integer JJ
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)
           
	real seb_val,swindow,EPS,tol,hightol,x1,x2
	real old,es0
	real a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1
	real xnew,Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf
	real dt,Qsip,Qsi,albedo,extcoef,xk_snow,ro_snow,Cp_snow
	real xk_water,total_solar,xLs,Rv,ro_water,ro_ice,water_depth
	real water_depth_old,water_flux,xLf,ro_pure_ice,Qc
	real xinternal_heating_corr
	integer maxiter,hr,model_day,icond_flag,kkk,ndarklayers,i

	EPS = 3.0e-8
	tol = 0.05
      hightol = 0.0001 !use higher tol when we converge close to or above Tf
	maxiter = 100
c	swindow = 30.0
c	a = Tair - swindow
c	b = Tair + swindow
	swindow = 15.0  ! 10 should be big enough when using tsfc
	a = xnew - swindow  ! last time's tsfc is a better guess than tair
	b = xnew + swindow

111	continue
c Calculate function at end points
	fa = seb_val(a,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea,
     &  icond_flag,gamma,T_old,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &      ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,Qc)
	fb = seb_val(b,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea,
     &  icond_flag,gamma,T_old,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &      ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,Qc)
	c = b
	fc = fb

	if (fa/abs(fa) .eq. fb/abs(fb)) then 
		swindow=swindow*2.0 ! double window each time
		a = xnew - swindow  ! last time's tsfc is a better guess than tair
		b = xnew + swindow
		print *,'Root not bracketed: 
     &		new window, day = ',swindow, model_day
			if (swindow.gt.60.0) then
				print *,'Window bigger than 60!'
				stop
			endif
		goto 111
	endif

c	if (fa/abs(fa) .eq. fb/abs(fb)) then
c	print *,'Root is not bracketed at model_day,iter=',model_day,i
c	swindow=swindow*2
c	stop
c	endif

      do i=1,maxiter ! start loop

	if ((fb.gt.0.0.and.fc.gt.0.0).or.(fb.lt.0.0.and.fc.lt.0.0)) then
c rename a,b,c and adjust bounding interval d
		c=a
		fc=fa
		d=b-a
		e=d
	endif
	if(abs(fc).lt.abs(fb)) then
		a=b
		b=c
		c=a
		fa=fb
		fb=fc
		fc=fa
	endif
	tol1=2.0*EPS*abs(b)+0.5*tol !convergence check

	xm=0.5*(c-b)
c increase tol if the roots bracket Tf or if both roots are close to or above Tf (melting)
	if(abs(xm).le.tol1 .and. 
     &	( (b-Tf)*(c-Tf).lt.0.0 .or. 
     &     (b.gt.Tf-1.0 .and. c.gt.Tf-1.0))  ) then
			tol=hightol !increase tol if we have converged close to Tf
			tol1=2.0*EPS*abs(b)+0.5*tol !update convergence check
	endif
	if(abs(xm).le.tol1 .or. fb.eq.0.0) then
		xnew =b		! FOUND THE ROOT
c			write(79,*) i
		return
	endif

	if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
		s=fb/fa  !attempt inv quad interp
		if(a.eq.c) then
			p=2.0*xm*s
			q=1.0-s
		else
			q=fa/fc
			r=fb/fc
			p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0))
			q=(q-1.0)*(r-1.0)*(s-1.0)
		endif

		if(p.gt.0.0) q=-q  !check whether in bounds
		p=abs(p)
		if(2.0*p .lt. min(3.0*xm*q-abs(tol1*q),abs(e*q))) then
			e=d  !accept interpolation
			d=p/q
		else
			d=xm  !interp failed - use bisection instead
			e=d
		endif
	else  ! bounds decreasing too slowly - use bisection
		d=xm
		e=d
	endif

	a=b !move last best guess to a
	fa=fb
	if(abs(d) .gt. tol1) then  !evaluate new trial root
		b=b+d
	else
		b=b+sign(tol1,xm)
	endif
	fb = seb_val(b,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea,
     &  icond_flag,gamma,T_old,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &      ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,Qc)
	enddo

	print *, 'zbrent exceeding maximum iteration!!!!!!!!!!!!!!!!!!'
	xnew=b
c		write (79,*) i
	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc hoffman's function to calculate the seb, used in BRENT_SOLVE

	FUNCTION seb_val(T_sfc,Tair,B1,B2,AAA,CCC,DDD,EEE,FFF,ea,
     &  icond_flag,gamma,T_old,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &      ro_snow,Cp_snow,xk_water,water_frac,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt,ro_ice,water_depth,
     &        water_depth_old,water_flux,xLf,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,Qc)

	implicit none

	integer JJ
      real gamma(JJ+2)
      real f_n(JJ+2)
      real dely_p(JJ+2)
      real dy_p(JJ+2)
      real y_crds(JJ+2)
      real T_old(JJ+2)
      real xmelt(JJ+2)
      real water_frac(JJ+2)
      real up(JJ+2)
      real down(JJ+2)

	real gamma_temp(JJ+2)
	real T_old_temp(JJ+2)
	real xmelt_temp(JJ+2)
	real water_frac_temp(JJ+2)
	           
	real T_sfc,es0
	real Tair,ea,AAA,CCC,DDD,EEE,FFF,B1,B2,Tf
	real dt,Qsip,Qsi,albedo,extcoef,xk_snow,ro_snow,Cp_snow
	real xk_water,total_solar,xLs,Rv,ro_water,ro_ice,water_depth
	real water_depth_old,water_flux,xLf,ro_pure_ice
	real xinternal_heating_corr
	integer icond_flag,kkk,ndarklayers,k,ifinalcall
	real seb_val
	real A,B,C,B3,B8,stability,water_flux_temp,water_depth_old_temp
	real T_sfcguess,Qc,qsfactor

           
c Over water.
	if (T_sfc.ge.Tf) then
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
	else
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55
	endif

c This section accounts for an increase in turbulent fluxes
c   under unstable conditions.
        es0 = A * exp((B * (T_sfc - Tf))/(C + (T_sfc - Tf)))
      if (T_sfc.gt.Tair) then
c Unstable case.
        B3 = 1.0 + B2 * sqrt(T_sfc - Tair)
        stability = 1.0 + B1 * (T_sfc - Tair) / B3
      elseif (T_sfc.lt.Tair) then
c Stable case.
        B8 = B1 / 2.0
        stability = 1.0 / ((1.0 + B8 * (Tair - T_sfc))**2)
      else
c Neutrally stable case.
        stability = 1.0
      endif

c !!! Calculate the ice temp profile with today's Tsfc guess, in order 
c to get a (hopefully) more stable value for Qc

c make copies of all variables that get written to within ICE_ENERGY
	do k=1,JJ+2
		gamma_temp(k)=gamma(k)                             
		T_old_temp(k)=T_old(k)
		xmelt_temp(k)=xmelt(k)
		water_frac_temp(k)=water_frac(k)
	enddo
	water_flux_temp=water_flux
	water_depth_old_temp=water_depth_old

c Be sure that ICE_ENERGY and CONDUCT aren't using a Tsfc>0
	if (T_sfc.gt.Tf) then
	T_sfcguess=Tf
	else
	T_sfcguess=T_sfc
	endif

	ifinalcall = 0 ! Let ICE_ENERGY know that this isn't the final call
	qsfactor = -99999.0 ! This can be a dummy value - make it crash if used

      CALL ICE_ENERGY(gamma_temp,T_old_temp,T_sfcguess,JJ,dy_p,y_crds,
     &        dt,f_n,dely_p,Qsip,Qsi,albedo,extcoef,xk_snow,
     &     ro_snow,Cp_snow,xk_water,water_frac_temp,up,down,total_solar,
     &        xLs,Rv,Tf,ro_water,xmelt_temp,ro_ice,water_depth,
     &        water_depth_old_temp,water_flux_temp,xLf,ro_pure_ice,
     &  xinternal_heating_corr,ndarklayers,ifinalcall,qsfactor,Qc)

cc      CALL CONDUCT(icond_flag,Qc,gamma_temp,T_old_temp,dy_p,
cc     &	JJ,T_sfcguess)


	  seb_val = EEE - DDD*T_sfc**4 + AAA*(Tair-T_sfc)*stability +
     &    FFF*CCC*(ea-es0)*stability + Qc + 0.0

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE LATENT(Qe,D_e,stability,ea,es0,ro_air,xLs,Pa)

      Qe = ro_air * xLs * D_e * stability * (0.622/Pa * (ea - es0))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	SUBROUTINE EXCOEFS(D_h,D_e,z_0,z_windobs,windspd,xkappa,theday)
	implicit none

	integer theday
	real D_h,D_e,z_0,z_windobs,windspd,xkappa

c	z_0 = 0.0
c	model_day = tempvar
c Taylor Glacier:
c	if ((theday .le. 104) .or. (theday .gt. 244)) then
c		z_zero = 0.0001
c	elseif ((theday .gt. 104) .and. (theday .le. 174)) then
c		z_zero = 0.00001414 * (theday - 104) + 0.0001
c	elseif ((theday .gt. 174) .and. (theday .le. 244)) then
c		z_zero = -0.00001414 * (theday -174) + 0.001
c	endif
		
	D_h = (xkappa**2) * windspd / ((log(z_windobs/z_0))**2)
	D_e = D_h
	return
	end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c mgc not sure why but hoffman added this

	SUBROUTINE EXCOEFS_SCALAR(D_h,D_e,z_0,z_windobs,windspd,xkappa,
     &     	theday)
	implicit none

	integer theday
	real D_h,D_e,z_0,z_windobs,windspd,xkappa
	real visc_air,alpha_h,alpha_q,Re,u_fric,z_t,z_q
	real CD,CH,CE
c	z_0 = 0.0
c	model_day = tempvar
c Taylor Glacier:
c	if ((theday .le. 104) .or. (theday .gt. 244)) then
c		z_zero = 0.0001
c	elseif ((theday .gt. 104) .and. (theday .le. 174)) then
c		z_zero = 0.00001414 * (theday - 104) + 0.0001
c	elseif ((theday .gt. 174) .and. (theday .le. 244)) then
c		z_zero = -0.00001414 * (theday -174) + 0.001
c	endif
		

	visc_air=1.461e-5
	alpha_h=1.0
	alpha_q=1.0
c Note: u_fric should account for stability!!!!!
	u_fric=xkappa*z_0 / (log(z_windobs/z_0))
	Re=u_fric*z_0/visc_air
	if (Re.lt.0.135) then
		z_t=z_0*exp(1.250)
		z_q=z_0*exp(1.610)
	elseif (Re.lt.2.5) then
		z_t=z_0*exp(0.149-0.55*log(Re))
		z_q=z_0*exp(0.351-0.628*log(Re))
	else
		z_t=z_0*exp(0.317-0.565*log(Re)-0.183*(log(Re))**2 )
		z_q=z_0*exp(0.396-0.512*log(Re)-0.180*(log(Re))**2 )
	endif
	CD=xkappa**2/(log(z_windobs/z_0))**2
	CH=alpha_h*xkappa*sqrt(CD)/(xkappa*CD**(-0.5) - log(z_t/z_0))
	CE=alpha_q*xkappa*sqrt(CD)/(xkappa*CD**(-0.5) - log(z_q/z_0))
c	Qh=ro_air*Cp*CH*wspd*stability*(Tair-T_sfc)
c	Qe=ro_air*xLs*CE*wspd*stability*CCC*(ea-es0)
	D_h = CH * windspd
	D_e = CE * windspd

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c mgc hoffman added statements to turn on water if Tsfc>Tf and ice if not

      SUBROUTINE VAPOR(es0,Tsfc,Tf)

c Coeffs for saturation vapor pressure over water (Buck 1981).
c   Note: temperatures for Buck's equations are in deg C, and
c   vapor pressures are in mb.  Do the adjustments so that the
c   calculations are done with temperatures in K, and vapor
c   pressures in Pa.

	if (Tsfc.ge.Tf) then
c Over water.
        A = 6.1121 * 100.0
        B = 17.502
        C = 240.97
	else
c Over ice.
        A = 6.1115 * 100.0
        B = 22.452
        C = 272.55
	endif

c Compute the water vapor pressure at the surface.
      es0 = A * exp((B * (Tsfc - Tf))/(C + (Tsfc - Tf)))

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c EXTINCTION COEFFICIENT SECTION 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE EXTCOEFS(nzext,deltazext,albedo,ro_snow,upext,downext,
     &  n_snowgrain_radius,total_solar,ro_pure_ice,y_crdsext,runname)

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

      integer nz,nvalues,nclasses,nzext,n_snowgrain_radius

c The number of z levels to use internally in the EXTCOEFS section
	parameter (nz=15000)
c The number of wavelength bands that are used.
      parameter (nvalues=118)

c The number of grain radii that can be used.
      parameter (nclasses=47)

      real ro_snow,r_snow,Qsi,albedo,total_solar

	real :: deltaz=.001
	real deltazext(nzext)
	real y_crdsext(nzext+2)

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

      real downext(nzext+2)
      real upext(nzext+2)
	real down(nz+2)
	real up(nz+2)
      real rad(nz)

      real a(nz+2)
      real r(nz+2)

	real Sc(nz)
	real xydiff(nz)
	real xynet(nz+2)

	real ro_pure_ice

	integer kk,kkext,j

	character*(80) runname
	integer runnamelen
	integer :: strlen

c These are the possible radii that can be used for the model
c   simulations.  They are in mm, and must be converted to meters
c   before they can be used in the model.  To pick a radius, you
c   just pick the array position corresponding to the radius
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
      write (6,102) r_snow * 1000.0
  102 format ('you have picked a grain radius (mm) = ',f9.3)

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
     &  r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses,ro_pure_ice)

c	open (76,file='./output/specextcoef.OUT')
c	do kk=1,nvalues
c			write (76,*) wavelength(kk,n_snowgrain_radius),
c     &	spect_extcoef_snow(kk),qext(kk,n_snowgrain_radius),
c     &ss_coalb(kk,n_snowgrain_radius),g(kk,n_snowgrain_radius)
c	enddo
c	close (76)


c Compute the downward bulk extinction coefficient.
      CALL BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen,
     &  z_without_bc,spect_extcoef_snow,bulk_ext_snow,
     &  bulk_ext_snowbc)


c	open (77,file='./output/bulkextcoef.OUT')
c	do kk=1,nz+2
c		write (77,*) kk,bulk_ext_snowbc(kk)
c	enddo
c	close (77)

c Compute the a and r coefficients from knowledge of the
c   surface albedo and the extinction coefficient.
      CALL GETAANDR(nz,bulk_ext_snowbc,a,r,albedo)

c Solve the system of equations.
      CALL SOLVETWOSTREAM(nz,a,r,deltaz,bulk_ext_snowbc,Qsi,rad,
     &  up,down)

c Optional output of down and up streams.  This is useful for getting around
c the variable deltaz problem: Set deltaz to be constant, output down & up, change
c deltaz back to the desired values, then read down & up in to bypass their calculation.
	open (24,file='./output/downup.OUT')
	do kk=1,nz+2
		write (24,*) z_with_bc(kk),down(kk),up(kk),down(kk)-up(kk)
	enddo
	close (24)

c Convert up and down streams to the external grid system.
c Set the upper & lower bdy.
	upext(1)=up(1)
	downext(1)=down(1)

c Find matching values from top down
	kkext=2
	do kk=2,nz+2
		if (z_with_bc(kk).ge.y_crdsext(kkext)) then
c We have a match! Linear interpolate to find value.
			upext(kkext) = up(kk-1) + (up(kk)-up(kk-1)) *
     &	(y_crdsext(kkext)-z_with_bc(kk-1)) / 
     &	(z_with_bc(kk)-z_with_bc(kk-1))
			downext(kkext) = down(kk-1) + (down(kk)-down(kk-1)) *
     &	(y_crdsext(kkext)-z_with_bc(kk-1)) / 
     &	(z_with_bc(kk)-z_with_bc(kk-1))
			kkext=kkext+1
		end if
c break out when external values have been calc'ed
			if (kkext.gt.nzext+2) then
				goto 137
			endif
	end do

137	continue
	upext(nzext+2)=up(nz+2)
	downext(nzext+2)=down(nz+2)

c ----Write out the interpolated up/down-----

c Provide the source terms.
c Build an array of values on the c.v. boundaries.
      xynet(1) = upext(1) - downext(1)
      do j=2,nzext
        xynet(j) = (upext(j) + upext(j+1))/2.0 - (downext(j) + 
     &	  downext(j+1))/2.0
      enddo
      xynet(nzext+1) = upext(nzext+2) - downext(nzext+2)

      do j=1,nzext
c This is dq/dz for q=Q0*exp(-extcoef*z).
c       Sc(j) = Qsip * extcoef * exp(- extcoef * y_crds(j+1))

c This is dq/dz for q=eqn 9 in Liston et al. 1999, where the values
c   are scaled by the ratio of the incoming solar to the solar
c   used to get up and down.
c After further review, mh change this to first version (see J 1119)
c        Sc(j) = - (xynet(j) - xynet(j+1)) / dy_p(j)
ccc        Sc(j) = - Qsi / total_solar * (xynet(j) - xynet(j+1)) / dy_p(j)
	Sc(j) = (xynet(j) - xynet(j+1)) / deltazext(j)
	xydiff(j) = (xynet(j) - xynet(j+1))
	end do

	runnamelen=strlen(runname)

	open (24,file='./output/'//runname(1:runnamelen)// 
     &	'/' // 'downupext.out')
	write (24,*) total_solar
	do kkext=1,nzext+2
c		write (24,*) y_crdsext(kkext),downext(kkext),upext(kkext),
c     &  xydiff(kkext)
		write (24,*) downext(kkext),upext(kkext)
	enddo
	close (24)

c write out the coordinate system used for the T points, excluding boundaries
	open (24,file='./output/'//runname(1:runnamelen)// 
     &	'/' // 'ycrds.out')
	do kkext=1,nzext+2
		write (24,*) y_crdsext(kkext)
	enddo
	close (24)

c	open (24,file='./output/sourceterm.txt')
c	do kkext=1,nzext
c		write (24,*) y_crdsext(kkext+1),Sc(kkext),xydiff(kkext)
c	enddo
c	close (24)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
      open (46,file='./input/mie.dat',form='formatted')
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

cc Jon checks values of ss_coalb and changes 0 to 1e-7
	do i=1,nvalues
	do j=1,nclasses
		if (ss_coalb(i,j) .eq. 0 ) then
			ss_coalb(i,j) = 1e-7
		end if
	enddo
	enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETAANDR(nz,bulk_ext_snowbc,a,r,albedo)

      implicit none

      integer k,nz

	real albedo

      real bulk_ext_snowbc(nz+2)
      real a(nz+2)
      real r(nz+2)

c	albedo = 0.65
      
c Compute the a and r coefficients from knowledge of the
c   albedo and the bulk extinction coefficient.
      do k=1,nz+2
        a(k) = (1.0 - albedo) /
     &    (1.0 + albedo) * bulk_ext_snowbc(k)
        r(k) = 2.0 * albedo * bulk_ext_snowbc(k) /
     &    (1.0 - albedo**2)
c Jon had this incorrect formula
c	r(k) = 2.0 * albedo * bulk_ext_snowbc(k) /
c     &    (albedo**2)
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETZ(z_with_bc,z_without_bc,deltaz,nz)

      implicit none

      integer k,nz

      real deltaz
      real z_without_bc(nz)
      real z_with_bc(nz+2)

c Make a z-depth array, in meters.  These are z-levels, the centers
c   of each grid cell.  z_with_bc includes the top and bottom edges of
c   the top and bottom grid cells. 
ccc ***  For constant dz
	z_with_bc(1)=0.0
	do k=2,nz+1
	z_with_bc(k)=deltaz*(k-1)-deltaz/2.0
	z_without_bc(k-1) = z_with_bc(k)
	enddo
	z_with_bc(nz+2)=deltaz*(nz)

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE BULKEXTCOEF(deltaz,nz,solar,nvalues,dwavelen,
     &  z_without_bc,spect_extcoef_snow,bulk_ext_snow,
     &  bulk_ext_snowbc)

      implicit none

      integer k,kk,nvalues,nz

      real sum1,sum2
	real deltaz
      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)
      real z_without_bc(nz)
      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)
	real default_bulk_ext

	default_bulk_ext=0.0

c Compute the downward bulk extinction coefficient.
      do kk=1,nz
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nvalues
          sum1 = sum1 + solar(k) *
     &      exp(- spect_extcoef_snow(k) * (z_without_bc(kk)+deltaz))
     &      * dwavelen(k)
          sum2 = sum2 + solar(k) *
     &      exp(- spect_extcoef_snow(k) * z_without_bc(kk)) *
     &      dwavelen(k)
      enddo

	if (sum2.ge.1e-27) then
      bulk_ext_snow(kk) = - (1.0 / deltaz) * log(sum1/sum2)
	endif

	if ((sum1.le. 1e-20).or.(sum2.le.1e-20)) then
		if (default_bulk_ext .lt. 0.1) then
			default_bulk_ext = bulk_ext_snow(kk)
		endif
		bulk_ext_snow(kk)=default_bulk_ext
	endif

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

      SUBROUTINE BULKEXTCOEF2(deltaz,nz,solar,nvalues,dwavelen,
     &  z_without_bc,spect_extcoef_snow,bulk_ext_snow,
     &  bulk_ext_snowbc)

c MH: This is my attempt to rewrite the calculation of the bulk ext coeff
c in terms of the 'net bulk-extinction coefficient', KGM (Grenfell & Maykut, 1977)
c See Brandt & Warren, 1993, p. 101
      implicit none

      integer k,kk,nvalues,nz

      real sum1,sum2,default_bulk_ext
	real deltaz(nz)
      real spect_extcoef_snow(nvalues)
      real solar(nvalues)
      real dwavelen(nvalues)
      real z_without_bc(nz)
      real bulk_ext_snow(nz)
      real bulk_ext_snowbc(nz+2)
	real specalb(nvalues)

	open (47,file='./input/specalb.dat',form='formatted')
      read (47,302) specalb
  302 format(118f5.3)

	default_bulk_ext=0.0      

c Compute the downward bulk extinction coefficient.
      do kk=1,nz
	if (kk.ge.142) then
c		hi=1
	endif
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nvalues
          sum1 = sum1 + ( spect_extcoef_snow(k)*(1-specalb(k)) ) * 
     &		solar(k) *
     &      exp(- spect_extcoef_snow(k) * z_without_bc(kk) )
     &      * dwavelen(k)
          sum2 = sum2 + ( (1-specalb(k)) ) * 
     &		solar(k) *
     &      exp(- spect_extcoef_snow(k) * z_without_bc(kk) )
     &      * dwavelen(k)
        enddo

	if (sum2.ge.1e-30) then
      bulk_ext_snow(kk) = (sum1/sum2)
	endif

	if ((sum1.le. 1e-25).or.(sum2.le.1e-25)) then
		if (default_bulk_ext .lt. 0.1) then
			default_bulk_ext = bulk_ext_snow(kk)
		endif
		bulk_ext_snow(kk)=default_bulk_ext
	endif




	if ((sum1.eq.0.0).or.(sum2.eq.0.0)) then
		bulk_ext_snow(kk)=6.281
	else
        bulk_ext_snow(kk) = (sum1/sum2)
	endif

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

      SUBROUTINE SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow,
     &  r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses,ro_pure_ice)

      implicit none

      integer k,n_snowgrain_radius,nvalues,nclasses

      real ro_snow,ro_ice,sigma_e,r_snow
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real spect_extcoef_snow(nvalues)

	real ro_pure_ice

      ro_ice = 917.0

      do k=1,nvalues
        sigma_e = 3.0/4.0 * qext(k,n_snowgrain_radius)/r_snow *
     &    ro_snow/ro_pure_ice 
c        sigma_e = 3.0/4.0 * qext(k,n_snowgrain_radius)/r_snow *
c     &    ro_snow/ro_ice
        spect_extcoef_snow(k) = sigma_e *
     &    sqrt(ss_coalb(k,n_snowgrain_radius) -
     &    ss_coalb(k,n_snowgrain_radius) * g(k,n_snowgrain_radius) +
     &    ss_coalb(k,n_snowgrain_radius)**2 * g(k,n_snowgrain_radius))
      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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

cc      do k=1,nvalues
c       print *,k,dwavelen(k),wavelength(k,1)
c       print *,wavelength(k,1)
cc      enddo

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

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
      open (76,file='./input/solar.dat',form='formatted')
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
c      print *, 'total solar radiation = ',total_solar

      close (76)

      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE GETUPDOWN(a,r,deltaz,nz,rad,up,down,Qsi)

      implicit none

      integer k,nz

      real alfa,Qsi
	real deltaz
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
c Do the ends. note: Glen:all denom 1.5, jon:.5,5,.5,5
      down_tmp(1) = (down(2) + 0.5 * down(1)) / 1.5
      down_tmp(nz+2) = (down(nz+1) + 0.5 * down(nz+2)) / 1.5
      up_tmp(1) = (up(2) + 0.5 * up(1)) / 1.5
      up_tmp(nz+2) = (up(nz+1) + 0.5 * up(nz+2)) / 1.5
	down_tmp(1)=down(1)
	up_tmp(1)=up(1)

c MH: I have changed the smoothing. Smoothing and then adjusting by Qsi at surface
c led to increasing the total absorbed rad by a not insignificant amount.
c Intead I smooth all but z=0.  This ensures that we get Qsi at the surface, 
c but does take out a few of the kinks in the first few depths.
c Rewrite the arrays.
      do k=1,nz+2
        down(k) = down_tmp(k)
        up(k) = up_tmp(k)
      enddo
c Now adjust to get back the Qsi at the surface.
c      do k=1,nz+2
c        down(k) = down(k) * Qsi / down_tmp(1)
c        up(k) = up(k) * Qsi / down_tmp(1)
c      enddo

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE SOLVETWOSTREAM(nz,a,r,deltaz,xmu,Qsi,rad,
     &  up,down)

c Note: The matrix is opposite in sign from Schlatter eqn A6

      implicit none

      integer k,nz

      real alfa,xk1,gamma,Qsi
	real deltaz
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
	real tmp1,tmp2,tmp3,tmp4

c Compute matrix diagonal and b coeffs.
      do k=1,nz
		tmp1 = (2.0 + deltaz**2 * xmu(k+1)**2)
		tmp2 = deltaz / (2.0 * r(k+1))
		tmp3 = a(k+1) * (r(k+2) - r(k))
		tmp4 = r(k+1) * (a(k+2) - a(k))
c        A_main(k) = (2.0 + deltaz(k)**2 * xmu(k+1)**2) -
c     &    (deltaz(k) / (2.0 * r(k+1)) * (a(k+1) * (r(k+2) - r(k)) -
c     &    r(k+1) * (a(k+2) - a(k))))
		A_main(k) = tmp1 - (tmp2 * (tmp3 - tmp4))
        b_vector(k) = 0.0
      end do

c Account for the boundary conditions on the main diagonal.
      alfa = 1.0 / (a(1) + r(1))
      xk1 = 1.0 + (r(3) - r(1)) / (4.0 * r(2))
      gamma = xk1 * alfa / (deltaz + alfa)
      A_main(1) = A_main(1) - gamma

c Modify b to account for the boundary conditions.
      b_vector(1) = xk1 * alfa * r(1) * Qsi * deltaz / (deltaz 
     &	+ alfa)
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

      SUBROUTINE SNOWEVENT(model_day,
     &  snow_cover_depth_old,i_yearstart,ro_snow_on_top,ro_water,
     &  Qsi,Qsi_fraction,albedo,rh,windspd)

c This Sub includes mods dealing with snow that Jon has added

      implicit none

	integer conditional_snow,model_day,i_yearstart,k
	real snow_cover_depth_old,ro_snow_on_top,ro_water
c	real stake_depth(25,3)
	real stake_density(1995:2005)
	real Qsi,Qsi_fraction,albedo,rh,windspd

c	read (34,*) conditional_snow

c Include conditional routine to determine the presence of a snow event.
c If the data are available, determine actual thicknesses of the snowfalls
c otherwise use arbitrary value such as 1 cm of snow.
	conditional_snow = 0

	if (Qsi .gt. 10.0) then
		if ((albedo .gt. 0.60) .and. (Qsi_fraction .lt. 0.57)) then
			conditional_snow = 1
		endif
	else
		if ((rh .gt. 73.0) .and. (windspd .lt. 4.0)) then
			conditional_snow = 1
		endif
	endif
	
	if (conditional_snow.eq.1) then
		if ((model_day.le.153).or.(model_day.ge.216)) then
c			snow_cover_depth_old = snow_cover_depth_old + 0.25
		elseif ((model_day.ge.154).and.(model_day.le.215)) then
c			snow_cover_depth_old = snow_cover_depth_old + 1.00
		endif
	endif

c Jon: Include the snow cover depths data obtained during ablation stake
c measurements to more accurately obtain ablation values for the
c meteorological station cell.
c This could be read in from a data file perhaps.

c Year,Day,snow depth data:
c CAA
c	stake_depth(1,:) = (/1998.0, 1.0, 0.0/)
c	stake_depth(2,:) = (/1998.0, 129.0, 12.4/)
c	stake_depth(3,:) = (/1998.0, 197.0, 9.9/)
c	stake_depth(4,:) = (/1999.0, 125.0, 0.0/)
c	stake_depth(5,:) = (/1999.0, 200.0, 0.0/)
c	stake_depth(6,:) = (/2000.0, 135.0, 0.0/)
c	stake_depth(7,:) = (/2000.0, 200.0, 3.9/)
c	stake_depth(8,:) = (/2001.0, 133.0, 0.0/)
c	stake_depth(9,:) = (/2001.0, 200.0, 0.0/)
c	stake_depth(10,:) = (/2002.0, 136.0, 1.6/)
c	stake_depth(11,:) = (/2002.0, 203.0, 0.0/)
c	stake_depth(12,:) = (/2003.0, 134.0, 0.0/)
c	stake_depth(13,:) = (/2003.0, 205.0, 3.2/)
c	stake_depth(14,:) = (/2004.0, 129.0, 0.0/)
c	stake_depth(15,:) = (/2004.0, 200.0, 7.2/)
c
c TAR
c	stake_depth(1,1:3) = (/1995.0, 139.0, 0.0/)
c	stake_depth(2,1:3) = (/1995.0, 172.0, 0.0/)
c	stake_depth(3,1:3) = (/1995.0, 203.0, 0.0/)
c	stake_depth(4,1:3) = (/1996.0, 136.0, 0.0/)
c	stake_depth(5,1:3) = (/1996.0, 156.0, 0.0/)
c	stake_depth(6,1:3) = (/1996.0, 205.0, 0.0/)
c	stake_depth(7,1:3) = (/1997.0, 134.0, 0.0/)
c	stake_depth(8,1:3) = (/1997.0, 207.0, 0.0/)
c	stake_depth(9,1:3) = (/1998.0, 136.0, 0.0/)
c	stake_depth(10,1:3) = (/1998.0, 202.0, 0.0/)
c	stake_depth(11,1:3) = (/1999.0, 126.0, 0.0/)
c	stake_depth(12,1:3) = (/1999.0, 134.0, 0.0/)
c	stake_depth(13,1:3) = (/1999.0, 201.0, 0.0/)
c	stake_depth(14,1:3) = (/2000.0, 139.0, 0.0/)
c	stake_depth(15,1:3) = (/2000.0, 203.0, 0.0/)
c	stake_depth(16,1:3) = (/2001.0, 138.0, 0.0/)
c	stake_depth(17,1:3) = (/2001.0, 207.0, 0.0/)
c	stake_depth(18,1:3) = (/2002.0, 138.0, 3.3/)
c	stake_depth(19,1:3) = (/2002.0, 201.0, 0.0/)
c	stake_depth(20,1:3) = (/2003.0, 137.0, 0.0/)
c	stake_depth(21,1:3) = (/2003.0, 201.0, 0.0/)
c	stake_depth(22,1:3) = (/2004.0, 135.0, 0.0/)
c	stake_depth(23,1:3) = (/2004.0, 203.0, 0.0/)
c	stake_depth(24,1:3) = (/2005.0, 130.0, 0.0/)
c	stake_depth(25,1:3) = (/2005.0, 201.0, 0.0/)


c	do k=1,size(stake_depth,1)
c		if ((i_yearstart .eq. stake_depth(k,1)) .and.
c     &			(model_day .eq. stake_depth(k,2))) then
c			snow_cover_depth_old = stake_depth(k,3)
c			exit
c		endif
c	enddo


c Snow density determined by data gathered during ablation stake
c meausrements.  The ro_snow_on_top variable is used only for ablation
c calculations and does not affect scattering within the ice.  All other
c years with no snow ro_snow = ro_ice.
c CAA
c	stake_density(1998) = 400.0
c	stake_density(1999) = 400.0
c	stake_density(2000) = 165.0
c	stake_density(2001) = 270.0
c	stake_density(2002) = 175.0
c	stake_density(2003) = 265.0
c	stake_density(2004) = 170.0

c TAR  only snow at met station was modelday 138, year 2002, value=300
	do k=1995,2005
		stake_density(k) = 300.0
	enddo

	if (model_day.gt.1) then
		ro_snow_on_top = ro_water
	endif

	if (snow_cover_depth_old.gt.0.0) then                
		  ro_snow_on_top = stake_density(i_yearstart)
	endif



      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	SUBROUTINE DARKENLAYERS(nz,dy_p,up,down,ndarklayers,qsfactor)	
c This subroutine automatically calculates what % of the net Qs to 
c eliminate from the SEB based on how many layers you want to black out.

	implicit none

	integer nz
	real dy_p(nz)
	real up(nz+2)
	real down(nz+2)
	real xynet(nz+2)
	real Sctotal,Scdark,Sc
	integer j,ndarklayers
	real qsfactor

	Sctotal=0.0
	Scdark=0.0

      xynet(1) = up(1) - down(1)
      do j=2,nz
        xynet(j) = (up(j) + up(j+1))/2.0 - (down(j) + 
     &	  down(j+1))/2.0
      enddo
      xynet(nz+1) = up(nz+2) - down(nz+2)

      do j=1,nz
		Sc=(xynet(j) - xynet(j+1)) / dy_p(j) * dy_p(j)
		Sctotal=Sctotal+sc
		if (j.le.ndarklayers) then
			Scdark=Scdark+Sc
		endif
	end do

	qsfactor=Scdark/Sctotal

	return
	end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      INTEGER FUNCTION strlen (st)
      integer i
      character st*(*)
      i = len(st)
      do while (st(i:i) .eq. ' ')
        i = i - 1
      enddo
      strlen = i
      END FUNCTION

