#
#	atmlndnl.pm			Erik Kluzek
#
#	Perl module to deal with namelists for the atmosphere and land.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#	setcfg -------------------- Set a given configuration value setting.
#	cfg ----------------------- Return a configuration value.
#	get_default_values -------- Set the default values from the Default XML file.
#	exists -------------------- Check if a cfg value exists.
#	checkinputfile ------------ Check if the given input file exists.
#	check --------------------- Check that the namelist is valid.
#
#	$Id: atmlndnl.pm,v 1.1.6.1 2002/05/13 17:21:23 eaton Exp $
#
use strict;
#use diagnostics;
use Cwd;
use XML::Lite;

package atmlndnl;
use namelist;
@atmlndnl::ISA = "namelist";
#
# Extend the namelist class to have a method to check that input files
# exist either on local disk, full-path given, or on Mass Store.
#
sub new {
#
# Constructor
#
  my $class = shift;
  my $name  = shift;
  my $refNL = shift;
  my $interactive = shift;
  my $file  = shift;
  my $defaults_file = shift;
  my $default_vals = shift;
  my $CFG = shift;
  my $printlev = shift;

  if ( ! defined($file) ) {
    die "ERROR($class): atmlndnl constructor was not sent the namelist filename\n";
  }

  my $self = $class->SUPER::new( $name, "$file", $refNL, $printlev );

  $self->{'DEFAULTS_FILE'} = $defaults_file;       # XML File with default values

  if ( ($interactive != 0) && ($interactive != 1) ) {
    die "ERROR($class): interactive option passed in to new was not valid: $interactive\n";
  }
  $self->{'INTERACTIVE'} = $interactive;           # Interactive mode (0 or 1)
  if ( ref($CFG) ne "GAMIL_config" ) {
    die "ERROR($class): Object sent to atmlndnl constructor not a GAMIL_config object\n";
  }

  $self->{CFG} = $CFG;        # GAMIL Configuration object
  $self->{MODEL_EXEDIR} = $CFG->cfg( "MODEL_EXEDIR" );    # Model execution directory
  $self->{MODEL_CFGDIR} = $CFG->cfg( "MODEL_CFGDIR" );    # Location of configuration files

  $self->{'default_vals'} = $default_vals;

  bless( $self, $class );
  return( $self );
}

#============================================================================

sub get_default_values {
#
#  Parse the Default XML file that gives most of the default settings
#  for different conditions (different resolutions, dynamics etcetera).
#
  my $self = shift;
  my $class = ref($self);
  my $nm = "$class\:\:get_default_values";

  my $EXPNLref = $self->{'default_vals'};
  my $MODEL_CFGDIR = $self->{"MODEL_CFGDIR"};
  my $file;
  if ( defined($MODEL_CFGDIR) ) {
    $file = $MODEL_CFGDIR . "/" . $self->{'DEFAULTS_FILE'};
  } else {
    $file = $self->{'DEFAULTS_FILE'};
  }
  print "($nm) Read: $file\n" if ($self->{'printlev'}>2);
  my $xml = XML::Lite->new( $file );
  if ( ! defined($xml) ) {
    die "ERROR($nm): Trouble opening or reading $file\n";
  }
  #
  # Find the namelist element for this namelist
  #
  my $elm = $xml->root_element( );
  my $namelist = $self->{'NAME'};
  my @list = $xml->elements_by_name( $namelist );
  if ( $#list < 0 ) {
    die "ERROR($nm): could not find the main $namelist namelist element in $file\n";
  }
  if ( $#list != 0 ) {
    die "ERROR($nm): $namelist namelist element in $file is duplicated, there should only be one\n";
  }
  #
  # Go through the sub-elements to the namelist element
  #
  $elm = $list[0];
  my @children = $elm->get_children();
  if ( $#children < 0 ) {
    die "ERROR($nm): There are no sub-elements to the $namelist element in $file\n";
  }
  foreach my $child ( @children ) {
    #
    # Get the attributes for each namelist element
    # The attributes describe either config settings that need to match
    # or other namelist elements that need to match
    #
    my %atts = $child->get_attributes;
    # Name of element, and it's associated value
    my $name = $child->get_name();
    my $value =  $child->get_text();
    $value =~ s/\n//g;   # Get rid of extra returns 
    # Expand the internal variables that might be in the string
    $value = $self->expand_vars_in_string( $value );
    my @keys = keys(%atts);
    my $match = 1;
    if ( $#keys >= 0 ) {
      #
      # Check that all values match the appropriate settings
      #
      foreach my $key ( @keys ) {
        # For config variables
        if ( $self->{CFG}->exists($key) && defined($self->{CFG}->cfg($key)) && 
             ($self->{CFG}->cfg($key) !~ /$atts{$key}/ ) ) {
          $match = 0;
          last;
        }
        # For namelist items
        if ( exists($$EXPNLref{$key}) && defined($$EXPNLref{$key}) &&
             ($$EXPNLref{$key} !~ /$atts{$key}/ ) ) {
          $match = 0;
          last;
        }
      }
    }
    # If match all attributes, and value isn't currently set
    if ( $match && ( ! exists($$EXPNLref{$name}) || 
      ! defined($$EXPNLref{$name}) ) ) {
      print "Set default value for: $name = $value\n" if ($self->{'printlev'}>2);
      $$EXPNLref{$name} = $value;
    }
  }
}

#============================================================================

sub do_interactive {
#
# Return true if interactive option set
#
  my $self = shift;

  my $value = $self->{INTERACTIVE};
  return( $value );
}

#============================================================================

sub checkinputfile {

# Check that the namelist value for an initial or boundary datasets is
# properly quoted.  Then check that the file exists on local filesystem.
# If the file is not found by looking at the full filepath, check for it in
# the directory where the GAMIL executable was created.

  my $self = shift;
  my $item = shift;

  my $class = ref($self);
  my $nm = "$class\:\:checkinputfile";

  my $EXPNLref = $self->{'NLREF'};
  my %EXPNL = %$EXPNLref;
  my $name = $EXPNL{$item};

  # check for quoting
  if ( $name !~ /["'](.*)['"]/ ) {
    die "$nm: $item needs quotes around filename: value = $name";
  }
  my $infile = $1;

  my $found_message = "Found $item dataset on local disk.";

  # check full pathname
  if ( -f $infile ) { 
      print "$found_message\n" if ($self->{'printlev'}>1);
      return;
  }

  # check for file in directory containing GAMIL executable
  $infile =~ /([^\/]+$)/;     # strip filename from the path
  my $file = $1;
  my $MODEL_EXEDIR = $self->{'MODEL_EXEDIR'};
  if ( defined($MODEL_EXEDIR) ) {
      if ( -f "$MODEL_EXEDIR/$file" ) { 
	  print "$found_message\n" if ($self->{'printlev'}>1);
	  return;
      }
  }

  print "Warning($nm): $item dataset $infile not found on local disk\n".
        "This dataset must be copied or linked to the run directory.\n";
}

#============================================================================

sub expand_vars_in_string {
#
# Expand any internal variables that are in a string
#
  my $self = shift;
  my $value = shift;

  while ( $value =~ /^(.*)\${*([a-zA-Z_]+[a-zA-Z0-9_]{0,19})}*(.*)$/ ) {
    my $var = $2;
    my $lead = $1;
    my $tail = $3;
    my $var_value;
    # If internal variable exists for this variable name
    if ( exists($self->{$var}) ) {
      $var_value = $self->{$var};
    } else {
      die "ERROR:: Internal variable $var needed in setting this value $value";
    }
    $value = "${lead}${var_value}${tail}";
  }
  return( $value );
}

#============================================================================

1   # to make use or require happy
