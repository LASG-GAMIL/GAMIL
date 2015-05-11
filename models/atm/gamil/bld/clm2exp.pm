#
#	clm2exp.pm			Erik Kluzek
#
#	Perl module to create a namelist for CLM2.
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	set_output_values --------  Set output values based on precedence of the various input
#                                   values and ensure that a valid namelist is produced.
#
#-----------------------------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package clm2exp;

use atmlndnl;
@clm2exp::ISA = qw(atmlndnl  namelist);

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

  my $self = $class->SUPER::new( "clmexp", \%main::CLMEXP, $interactive, $file, 
                                 "DefaultCLMEXPNamelist.xml", $default_vals,
                                 $GAMIL_config, $printlev );

  $self->{'printlev'} = $printlev;
  $self->{'optsref'}  = $optsref;

  $self->{DYNAMICS}   = $GAMIL_config->cfg("DYNAMICS");       # Dynamics to use
  $self->{RESOLUTION} = $GAMIL_config->cfg("RESOLUTION");     # horizontal resolution
  $self->{OCEANMODEL} = $GAMIL_config->cfg("OCEANMODEL");     # Ocean using with GAMIL (dom or som)

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the CLM2 namelist variables.

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

  # Check that "nrevsn" is set if this is a branch simulation
  if ($runtype eq 'branch' and !defined($NLref->{'nrevsn'})) {
      if ( $self->do_interactive ) {
	  print "Enter absolute pathname for CLM2 master restart file from which to branch: ";
	  $opt = <>; chomp $opt;
	  $NLref->{'nrevsn'} = namelist::quote_string($opt);
      } else {
	  die "ERROR: The CLM2 master restart file must be specified for a branch\n".
	  "       run.  Set the namelist variable NREVSN to the absolute\n".
	  "       pathname for this dataset.\n".
	  "       This can be done on the command-line using the -namelist\n".
          "       option or in an input namelist file that is specified\n".
          "       using the -infile option.\n";
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
  my $datdir = "$rootdir/lnd/clm2";

  # Plant function types.
  unless (defined($NLref->{'fpftcon'})) {
      $NLref->{'fpftcon'} = namelist::quote_string("$datdir/$default_vals->{'fpftcon'}");
  }
  $self->checkinputfile('fpftcon') if $optsref->{'test'};

  # Initial conditions
  unless ( defined($NLref->{'finidat'}) ) {
      if ( defined($default_vals->{'finidat'}) ) {
	  $NLref->{'finidat'} = namelist::quote_string("$datdir/$default_vals->{'finidat'}");
      }
  }
  if ( defined $NLref->{'finidat'} and $optsref->{'test'} ) {
      $self->checkinputfile('finidat');
  }

  # Surface datasets
  # Look for default unless user has set a value for fsurdat.  I
  unless ( defined($NLref->{'fsurdat'}) ) {
      if ( defined($default_vals->{'fsurdat'}) ) {
	  $NLref->{'fsurdat'} = namelist::quote_string("$datdir/$default_vals->{'fsurdat'}");
      }
  }
  # If user has set a "blank" value for fsurdat then remove the key so that 
  # setting the test option doesn't look for it, and the high resolution datasets
  # will be included in the namelist.
  if ($NLref->{'fsurdat'} =~ /'\s*'/) {
      delete $NLref->{'fsurdat'};
  }
  if ( defined $NLref->{'fsurdat'} and $optsref->{'test'} ) {
      $self->checkinputfile('fsurdat');
  }

  # High resolution surface datasets
  unless ( defined($NLref->{'fsurdat'}) ) {
      $NLref->{'mksrf_fvegtyp'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fvegtyp'}");
      $self->checkinputfile('mksrf_fvegtyp') if $optsref->{'test'};

      $NLref->{'mksrf_fsoitex'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fsoitex'}");
      $self->checkinputfile('mksrf_fsoitex') if $optsref->{'test'};

      $NLref->{'mksrf_fsoicol'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fsoicol'}");
      $self->checkinputfile('mksrf_fsoicol') if $optsref->{'test'};

      $NLref->{'mksrf_flanwat'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_flanwat'}");
      $self->checkinputfile('mksrf_flanwat') if $optsref->{'test'};

      $NLref->{'mksrf_furban'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_furban'}");
      $self->checkinputfile('mksrf_furban') if $optsref->{'test'};

      $NLref->{'mksrf_fglacier'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_fglacier'}");
      $self->checkinputfile('mksrf_fglacier') if $optsref->{'test'};

      $NLref->{'mksrf_flai'} = namelist::quote_string("$datdir/$default_vals->{'mksrf_flai'}");
      $self->checkinputfile('mksrf_flai') if $optsref->{'test'};
  }
}

#============================================================================

1   # to make use or require happy
