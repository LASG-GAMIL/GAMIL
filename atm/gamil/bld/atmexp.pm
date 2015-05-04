#
#	atmexp.pm			Erik Kluzek
#
#	Perl module to create a namelist for GAMIL.
#
#	Description of methods:
#
#	new ----------------------  Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package atmexp;

use atmlndnl;
@atmexp::ISA = qw(atmlndnl  namelist);

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;
  my $GAMIL_config = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'out'};
  my $printlev = $$optsref{'printlev'};

  my $default_vals = {};

  my $self = $class->SUPER::new( "atmexp", \%main::ATMEXP, $interactive, $file,
                                 "DefaultATMEXPNamelist.xml", $default_vals,
                                 $GAMIL_config, $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  $self->{DYNAMICS}   = $GAMIL_config->cfg("DYNAMICS");  # model dynamics
  $self->{PERGRO}     = $GAMIL_config->cfg("PERGRO");    # PERGRO error-growth option
  $self->{PLEV}       = $GAMIL_config->cfg("PLEV");      # Ocean/sea-ice directory to use (dom/som)
  $self->{PHYSICS}    = $GAMIL_config->cfg("PHYSICS");   # physics directory to use
  $self->{RESOLUTION} = $GAMIL_config->cfg("RESOLUTION");# horizontal resolution

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the ATM namelist variables.

  my $self = shift;
  my $runtype = shift;

  my $class = ref($self);
  my $nm = "$class\:\:set_default_values";

  my $NLref = $self->{'NLREF'};
  my $optsref = $self->{'optsref'};
  my $default_vals = $self->{'default_vals'};
  my $opt;

  # Get the default values from the XML file
  $self->get_default_values;

  # Run type
  if ($runtype eq 'initial') { $opt = 0; }
  elsif ($runtype eq 'restart') { $opt = 1; }
  elsif ($runtype eq 'branch') { $opt = 3; }
  $NLref->{'nsrest'} = $opt;

  # Check that "nrevsn" is set if this is a branch simulation
  if ($NLref->{'nsrest'}==3 and !defined($NLref->{'nrevsn'})) {
      if ( $self->do_interactive ) {
	  print "Enter absolute pathname for GAMIL master restart file from which to branch: ";
	  $opt = <>; chomp $opt;
	  $NLref->{'nrevsn'} = namelist::quote_string($opt);
      } else {
	  die "ERROR: The GAMIL master restart file must be specified for a branch\n".
	  "       run.  Set the namelist variable NREVSN to the absolute\n".
	  "       pathname for this dataset.\n".
	  "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option.\n";
      }
  }

  # Case name
  if (defined($optsref->{'case'})) {
      $opt = $optsref->{'case'};
  } elsif (defined($NLref->{'caseid'})) {
      $opt = $NLref->{'caseid'};
  } else {
      $opt = $default_vals->{'caseid'};
  }
  my $case = $opt;
  if ( $case !~ /([^\/]{1,32})([^\/]*)/ ) {
    die "$nm Bad casename: $case\n";
  }
  # If casename too long
  if ( $2 ne "" ) {
    $case = $1;
    print "\n\nWARNING:: Truncating caseid in namelist from $opt to $case\n\n" if ($self->{'printlev'});
  }
  $NLref->{'caseid'} = namelist::quote_string($case);

  # Length of simulation
  unless ( defined($NLref->{'nelapse'}) or defined($NLref->{'nestep'}) or
	   defined($NLref->{'stop_ymd'})) {
    $NLref->{'nelapse'} = $default_vals->{'nelapse'};
  }

  # Orbit
  unless ( defined($NLref->{'obliq'}) and defined($NLref->{'eccen'}) and
	   defined($NLref->{'mvelp'}) ) {
      unless ( defined($NLref->{'iyear_ad'}) ) {
	  $NLref->{'iyear_ad'} = $default_vals->{'iyear_ad'};
      }
  }

  # Diffusion for Eulerian dynamics
  unless (defined($NLref->{'dif2'})) {
      if ( defined($default_vals->{'dif2'}) ) {
	  $NLref->{'dif2'} = $default_vals->{'dif2'};
      }
  }
  unless (defined($NLref->{'dif4'})) {
      if ( defined($default_vals->{'dif4'}) ) {
	  $NLref->{'dif4'} = $default_vals->{'dif4'};
      }
  }

  # Timestep size
  unless (defined($NLref->{'dtime'})) {
      if ( defined($default_vals->{'dtime'}) ) {
	  $NLref->{'dtime'} = $default_vals->{'dtime'};
      }
  }

  # reset_csim_iceprops 
  unless (defined($NLref->{'reset_csim_iceprops'})) {
      if ( defined($default_vals->{'reset_csim_iceprops'}) ) {
	  $NLref->{'reset_csim_iceprops'} = $default_vals->{'reset_csim_iceprops'};
      }
  }

  # PERGRO settings
  if ($self->{PERGRO} eq 'TRUE') {
      unless (defined($NLref->{'prognostic_icesnow'})) {
	  $NLref->{'prognostic_icesnow'} = ".false.";
      }
  }

  # Root directory for default initial and boundary datasets
  my $rootdir;
  if (defined($optsref->{'csmdata'})) {
      $rootdir = $optsref->{'csmdata'};
  } elsif (defined $ENV{'CSMDATA'}) {
      $rootdir = $ENV{'CSMDATA'};
  } else {
      $rootdir = $default_vals->{'csmdata'};
  }
  my $am_datdir = "$rootdir/atm/gamil";

  # Initial conditions
  if ( defined($NLref->{'ncdata'}) ) {
      $opt = $NLref->{'ncdata'};
  } elsif ( defined($default_vals->{'ncdata'}) ) {
      $opt = "$am_datdir/$default_vals->{'ncdata'}";
  } else {
      if ( $self->do_interactive ) {
	  print "Enter absolute pathname for initial dataset: ";
	  $opt = <>; chomp $opt;
      } else {
	  die "ERROR: The initial dataset must be specified (no suitable\n".
	  "       default value was found).  Set the namelist variable\n".
	  "       NCDATA to the absolute pathname for this dataset.\n".
	  "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option.\n";
      }
  }
  $NLref->{'ncdata'} = namelist::quote_string($opt);
  $self->checkinputfile('ncdata') if $optsref->{'test'};

  # Absorptivity and emissivity data
  if ( defined($NLref->{'absems_data'}) ) {
      $opt = $NLref->{'absems_data'};
  } else {
      $opt = "$am_datdir/$default_vals->{'absems_data'}";
  }
  $NLref->{'absems_data'} = namelist::quote_string($opt);
  $self->checkinputfile('absems_data') if $optsref->{'test'};


  # SST dataset
  if ( defined($NLref->{'bndtvs'}) ) {
      $opt = $NLref->{'bndtvs'};
  } elsif ( defined($default_vals->{'bndtvs'}) ) {
      $opt = "$am_datdir/$default_vals->{'bndtvs'}";
  } else {
      if ( $self->do_interactive ) {
	  print "Enter absolute pathname for SST dataset: ";
	  $opt = <>; chomp $opt;
      } else {
	  die "ERROR: The SST dataset must be specified (no suitable\n".
	  "       default value was found).  Set the namelist variable\n".
	  "       BNDTVS to the absolute pathname for this dataset.\n".
	  "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option.\n";
      }
  }
  $NLref->{'bndtvs'} = namelist::quote_string($opt);
  $self->checkinputfile('bndtvs') if $optsref->{'test'};

  # Ozone dataset
  if ( defined($NLref->{'bndtvo'}) ) {
      $opt = $NLref->{'bndtvo'};
  } else {
      $opt = "$am_datdir/$default_vals->{'bndtvo'}";
  }
  $NLref->{'bndtvo'} = namelist::quote_string($opt);
  $self->checkinputfile('bndtvo') if $optsref->{'test'};


  # Greenhouse gas chemistry dataset
  if ( defined($NLref->{'trace_gas'}) && $NLref->{'trace_gas'} =~ /\.true\./i ) {
      if ( defined($NLref->{'bndtvg'}) ) {
	  $opt = $NLref->{'bndtvg'};
      } else {
	  $opt = "$am_datdir/$default_vals->{'bndtvg'}";
      }
      $NLref->{'bndtvg'} = namelist::quote_string($opt);
      $self->checkinputfile('bndtvg') if $optsref->{'test'};
  }

}

#============================================================================

sub print_hash {
    my %h = @_;
    my ($k, $v);
    while ( ($k,$v) = each %h ) { print "$k => $v\n"; }
}


1   # to make use or require happy
