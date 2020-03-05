c mgc compute the spectral extinction coefficient
c mgc NOTE: this is where i would update N
c mgc add ro_snow_z (removed)

c mgc PICK UP HERE ... realized at this point that by making ro_snow_z
c vertically resolved, it makes the extinction coefficient vertical,
c but then that needs to be dealt with in the bulk ext. calculation,
c and we don't have methods developed for that yet

      SUBROUTINE SPECTEXTCOEF(nvalues,n_snowgrain_radius,ro_snow,
     &  r_snow,qext,ss_coalb,g,spect_extcoef_snow,nclasses,ro_snow_z)

      implicit none

      integer k,n_snowgrain_radius,nvalues,nclasses
      real ro_snow,ro_ice,sigma_e,r_snow
      real g(nvalues,nclasses)
      real qext(nvalues,nclasses)
      real ss_coalb(nvalues,nclasses)
      real spect_extcoef_snow(nvalues,JJ)
	  real ro_snow_z(JJ+2)

      ro_ice = 917.0

	  do j=1,JJ
		  ro_z = ro_snow_z(j)
      
	  enddo

      return
      end