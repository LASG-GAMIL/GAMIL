      logical fixYear_ghg   ! true => Ramped gases fixed at specified year.
      logical fixYear_so4   ! true => Ramped gases fixed at specified year.
      logical fixYear_scon  ! true => Ramped gases fixed at specified year.
      common /ramp_l/ fixYear_ghg, fixYear_so4, fixYear_scon

      integer rampYear_ghg   ! ramped gases fixed at this year
      integer rampYear_so4   ! ramped gases fixed at this year
      integer rampYear_scon  ! ramped gases fixed at this year
      common /ramp_i/ rampYear_ghg, rampYear_so4, rampYear_scon


