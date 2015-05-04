#------------------------------------------------------------------------
#
#	mprun2d.pm			Erik Kluzek
#
#	Perl module to create a namelist for mprun2d.
#
#	Description of methods:
#
#	new ------------------ Constructor
#	set_output_values ---- Set output values based on precedence of the various input
#                              values and ensure that a valid namelist is produced.
#
#------------------------------------------------------------------------

use strict;
#use diagnostics;
use Cwd;

package mprun2d;

@mprun2d::ISA = qw(namelist);
use namelist;

#============================================================================

sub new {
#
# Constructor
#
  my $class = shift;
  my $optsref = shift;

  my $interactive = $$optsref{'interactive'};
  my $file = $$optsref{'out'};
  my $printlev = $$optsref{'printlev'};

  my $self = $class->SUPER::new( "mprun2d", $file, \%main::MPRUN2D, $printlev );

  $self->{'interactive'}  = $interactive;
  bless( $self, $class );
  return( $self );
}

#============================================================================

sub set_output_values {

# Set the MPRUN2D namelist variables.

  my $self = shift;
  my $runtype = shift;

  my $class = ref($self);
  my $nm = "$class\:\:set_default_values";

  my $NLref = $self->{'NLREF'};
  my $opt;

  # npr_yz
  if ( ! defined($NLref->{'npr_yz'}) ) {
      if ( $self->{'interactive'} ) {
         print "npr_yz not set in namelist\n";
         print "Enter 2D parallel decomposition (Y,Z,X',Y'): ";
         my $opt = <>; chomp $opt;
         $NLref->{'npr_yz'} = $opt;
      } else {
         die "ERROR($nm):: npr_yz not set in namelist\n";
      }
  }
  my @npr_y = split( /,/, $NLref->{'npr_yz'} );
  if ( $#npr_y != 3 ) {
     die "ERROR($nm):: Wrong number of values for the npr_yz array - there must be four values.\n";
  }
  foreach my $i ( @npr_y ) {
    if ( $i !~ /[0-9]+/ ) {
       die "ERROR($nm):: The values in the npr_yz array must be integers\n";
    }
  }
  if ( $npr_y[0]*$npr_y[1] != $npr_y[2]*$npr_y[3] ) {
     die "ERROR($nm):: product of the first two values in npr_yz must equal the product of the last two\n";
  }
  my $SPMD_NODES = $npr_y[0]*$npr_y[1];
  my $name = $self->{NAME};
  print "$name: namelist setup for $SPMD_NODES SPMD tasks\n"  if ($self->{'printlev'}>0);
}

#============================================================================

1   # to make use or require happy
