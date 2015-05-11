#-----------------------------------------------------------------------------------------------
#
# GAMIL_config.pm
#
# Object to read in configure config_cache.xml 
# configuration files.
#
# Methods:
#
# new ---------------- Constructor
# read_config_cache -- Read in the config_cache.xml file.
# exists ------------- Check if a configuration name exists
# setcfg ------------- Set configuration value
# unsetcfg ----------- unset configuration value
# cfg ---------------- Return configuration value
#
# Date        Author                  Modification
# ---------------------------------------------------------------------------
# 23.Apr.2002 Erik Kluzek             Original version
#-----------------------------------------------------------------------------------------------

package GAMIL_config;

use strict;
#use warnings;
#use diagnostics;
use XML::Lite;

sub new {
#
# Constructor
#
  my $class = shift;

  my $self = {};
  my @cfg_vars = (
         "GAMILROOT",         # GAMIL root location
         "MODEL_BLDDIR",    # Model build location
         "MODEL_EXEDIR",    # Model execution location
         "MODEL_CFGDIR",    # Model configuration/build location
         "EXENAME",         # Executable name
         "MODEL_MODDIR",    # Directory with modified code
         "COUP_CSM",        # If running coupled or not
         "DYNAMICS",        # Dynamics to use
         "ICEMODEL",        # Sea-Ice model to use
         "PHYSICS",         # Physics directory to use
         "LANDMODEL",       # Land-model to use
         "OCEANMODEL",      # Ocean-model to use
         "PERGRO",          # Error growth option
         "RESOLUTION",      # Horizontal resolution
         "PLON",            # No. of longitudes
         "PLAT",            # No. of latitudes
         "PLEV",            # No. of vertical levels
         "PCNST",           # No. of advected constituents
         "PNATS",           # No. of non-advected constituents
         "PTRM",            # M truncation
         "PTRN",            # N truncation
         "PTRK",            # K truncation
         "DEBUG",           # Compiler debug setting
         "SPMD",            # SPMD mode (MPI multitasking)
         "SMP",             # Shared Memory Processing mode
         "PCOLS",           # physics grid # of columns
         "TWOD_YZ",         # 2D decomposition for FV dynamics
         "INC_NETCDF",      # NetCDF include
         "LIB_NETCDF",      # NetCDF library
         "INC_MPI",         # MPI include
         "LIB_MPI"          # MPI library
  );
  $self->{CFG_LIST} = \@cfg_vars;
  foreach my $cfg ( @cfg_vars ) {
    $self->{$cfg} = undef;
  }
  bless( $self, $class );
}

sub read_config_cache {
#
# Read in the config_cache.xml file
#
  my $self = shift;
  my $file = shift;

  my $nm = ref($self) . "::read_config_cache";
  if ( ! defined($file) ) {
    die "$nm:: ERROR:: Configuration cache file not given to read_config_cache\n";
  }
  if ( ! -f "$file" ) {
    die "$nm:: ERROR:: Configuration cache file: $file does not exist!";
  }
  my $xml = XML::Lite->new( $file );
  my $root = $xml->root_element();

  # Check for valid root node
  my $name = $root->get_name();
  $name eq "config_bld" or die "$nm:: ERROR:: file $file is not a GAMIL configuration file\n";

  # Get source and build directories
  my $dirs = $xml->elements_by_name( "directories" );
  my %dirs = $dirs->get_attributes();
  my %list = ( GAMILROOT=>"gamil_root", MODEL_BLDDIR=>"gamil_bld", MODEL_EXEDIR=>"gamil_exedir", 
               MODEL_MODDIR=>"usr_src" );
  foreach my $item ( keys(%list) ) {
    if ( ! defined($dirs{$list{$item}}) ) {
      print "$nm: WARNING:: $list{$item} not included in directory element in $file\n";
      $self->unsetcfg($item);
    } else {
      $self->setcfg($item,  $dirs{$list{$item}} );
    }
  }
  $self->setcfg("MODEL_CFGDIR",  $dirs{gamil_root}."/models/atm/gamil/bld/" );

  # Get packages
  my $pkgs = $xml->elements_by_name( "packages" );
  my %pkgs = $pkgs->get_attributes();
  my %list = ( DYNAMICS=>"dyn", PHYSICS=>"phys", LANDMODEL=>"lnd", OCEANMODEL=>"ocn",
               PERGRO=>"pergro" );
  foreach my $item ( keys(%list) ) {
    if ( ! defined($pkgs{$list{$item}}) ) {
      print "$nm: WARNING:: $list{$item} not included in package element in $file\n";
      $self->unsetcfg($item);
    } else {
      $self->setcfg($item,  $pkgs{$list{$item}} );
    }
  }
  $self->setcfg("COUP_CSM",    "FALSE" );
  if ( defined($self->cfg("PERGRO")) ) {
    if ( $self->cfg("PERGRO") ) {
      $self->setcfg("PERGRO",  "TRUE" );
    } else {
      $self->setcfg("PERGRO",  "FALSE" );
    }
  }

  # Get resolution parameters
  my $resparms = $xml->elements_by_name( "resolution" );
  my %resparms = $resparms->get_attributes();
  my %list = ( RESOLUTION=>"res", PLON=>"nlon", PLAT=>"nlat", PLEV=>"nlev", 
               PCNST=>"nadv", PNATS=>"nnadv", PTRM=>"trm", PTRN=>"trn", 
               PTRK=>"trk", PCOLS=>"pcols" );
  foreach my $item ( keys(%list) ) {
    if ( ! defined($resparms{$list{$item}}) ) {
      print "$nm: WARNING:: $list{$item} not included in resolution element in $file\n";
      $self->unsetcfg($item);
    } else {
      $self->setcfg($item,  $resparms{$list{$item}} );
    }
  }
  if ( $self->cfg("DYNAMICS") eq "fv" ) {
    if ( ! defined($resparms{twod_yz}) ) {
      print "$nm: WARNING:: twod_yz not included in resolution element in $file\n";
      $self->unsetcfg("TWOD_YZ");
    } else {
      if ( $resparms{twod_yz} ) {
        $self->setcfg("TWOD_YZ",  "TRUE" );
      } else {
        $self->setcfg("TWOD_YZ",  "FALSE" );
      }
    }
  }

  # Get settings for Makefile (parallelism and library locations)
  my $make = $xml->elements_by_name( "makefile" );
  my %make = $make->get_attributes();
  my %list = ( EXENAME=>'gamil_exe', DEBUG=>"debug", SPMD=>"spmd", SMP=>"smp", 
               INC_NETCDF=>"nc_inc", LIB_NETCDF=>"nc_lib", INC_MPI=>"mpi_inc", 
               LIB_MPI=>"mpi_lib" );
  foreach my $item ( keys(%list) ) {
    if ( ! defined($make{$list{$item}}) ) {
      print "$nm:: WARNING:: $list{$item} not included in makefile element in $file\n";
      $self->unsetcfg($item);
    } else {
      $self->setcfg($item, $make{$list{$item}} );
    }
  }
  if ( defined($self->cfg("DEBUG")) ) {
    if ( $make{debug} ) {
      $self->setcfg("DEBUG",  "TRUE" );
    } else {
      $self->setcfg("DEBUG",  "FALSE" );
    }
  }
  if ( defined($self->cfg("SPMD")) ) {
    if ( $make{spmd} ) {
      $self->setcfg("SPMD",  "TRUE" );
    } else {
      $self->setcfg("SPMD",  "FALSE" );
    }
  }
  if ( defined($self->cfg("SMP")) ) {
    if ( $make{smp} ) {
      $self->setcfg("SMP",  "TRUE" );
    } else {
      $self->setcfg("SMP",  "FALSE" );
    }
  }
}


sub exists {
#
# Check if the given cfg variable is defined here
#
  my $self = shift;
  my $cfg = shift;

  if ( ! defined($cfg) ) {
    print "ERROR:: Configuration variable not given to exists method\n";
  }
  my $list = $self->{'CFG_LIST'};
  my @list = @$list;

  # Make sure this config var defined in list
  my $found = 0;
  foreach my $item ( @list ) {
    if ( $item eq $cfg ) { $found = 1; last; }
  }
  return( $found );
}

sub setcfg {
#
# Set a given cfg value setting
#
  my $self = shift;
  my $cfg = shift;
  my $value = shift;

  if ( defined($cfg) ) {
    if ( ! $self->exists($cfg) ) {
      die ref($self) . "::ERROR:: config variable $cfg not found in this class.\n";
    }
    if ( defined($value) ) {
      $self->{$cfg} = $value;
    } else {
      die ref($self) . "::ERROR:: value not passed into setcfg method for $cfg\n";
    }
  } else {
    die ref($self) . "::ERROR:: config variable not passed into setcfg method\n";
  }
}

sub unsetcfg {
#
# un-Set a given cfg value setting
#
  my $self = shift;
  my $cfg = shift;

  if ( defined($cfg) ) {
    if ( ! $self->exists($cfg) ) {
      die ref($self) . "::ERROR:: config variable $cfg not found in this class.\n";
    }
    $self->{$cfg} = undef;
  } else {
    die ref($self) . "::ERROR:: config variable not passed into setcfg method\n";
  }
}


sub cfg {
#
# Get the value of a given cfg variable
#
  my $self = shift;
  my $cfg = shift;

  if ( defined($cfg) ) {
    if ( ! $self->exists($cfg) ) {
      die \ref($self) . "::ERROR:: config variable $cfg not found in this class.\n";
    }
    my $value = $self->{$cfg};
    return( $value );
  } else {
    die \ref($self) . "::ERROR:: config variable not passed into this method\n";
  }
}

1 # to make use or require happy
