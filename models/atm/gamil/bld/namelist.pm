#
#	namelist.pm			Erik Kluzek
#
#	Perl module to deal with FORTRAN namelists.
#
#------------------------------------------------------------------------
#
#	Description of methods:
#
#	new ----------------------- Constructor
#	change -------------------- Interactively change values in the namelist
#	checkstring --------------- Check that a string namelist item is handled correctly.
#	parse --------------------- Parse a namelist file into a Perl associative array.
#	print --------------------- Print namelist to screen.
#	Write --------------------- Write the namelist out.
#	convert_case -------------- Convert the keys to lowercase.
#	quote_string -------------- Make sure a string has quotes around it.
#
#	Methods meant to be private to namelist:
#
#	split_namelist_value ------ Split namelist value up if to long.
#	setkeypair ---------------- Set a particular key-value pair into the namelist hash.
#	parse_next ---------------- Parse the next element in a line of a namelist.
#
#	$Id: namelist.pm,v 1.1.6.2 2002/07/31 19:55:16 erik Exp $
#
use strict;
#use diagnostics;
use Cwd;

package namelist;

# Some data to make global within this package, but local to inside it
#
# Perl expression to match a fortran variable
#
$namelist::varmatch = "[A-Za-z_]+[A-Za-z0-9_%()]*";    # % for derived types, () for arrays
#
# Perl pattern to match the value for a fortran constant
#
# Match for logical data
$namelist::vallogical = "\\.[Tt][Rr][Uu][Ee]\\.|\\.[Ff][Aa][Ll][Ss][Ee]\\.";
# Match for integer data
$namelist::valint = "[+-]?[0-9]+";
# Match for real data
# "_" are for f90 precision specification
$namelist::valreal = "$namelist::valint\\.?[0-9]*[EedDqQ]?[0-9+-]*_?[0-9a-z_]*";
# Match for string data
# One problem with below is strings that have \" or \' in them
$namelist::valstring = '\'[^\']+\'|"[^"]+"';
# Match for complex data
$namelist::valcomplex = "\\($namelist::valreal,$namelist::valreal\\)";
# Match for all valid data-types: integer, real, complex, logical, or string data
$namelist::valall   = 
"$namelist::vallogical|$namelist::valstring|$namelist::valreal|$namelist::valint|$namelist::valcomplex";
# Match for logical data with repeater on front
$namelist::vallogical_repeat = "${namelist::valint}\\*\.[Tt][Rr][Uu][Ee]\\.|${namelist::valint}\\*\.[Ff][Aa][Ll][Ss][Ee]\\.";
# Match for string data with repeater on front
$namelist::valstring_repeat = ${namelist::valint}.'\\*\'[^\']+\'|'.${namelist::valint}.'\\*\"[^\"]+\"';
# Match for all valid data-types with repeater: integer, real, complex, logical, or string data
$namelist::valrepeat = 
"$namelist::vallogical_repeat|$namelist::valstring_repeat|${namelist::valint}\\*$namelist::valreal|${namelist::valint}\\*$namelist::valint|${namelist::valint}\\*$namelist::valcomplex";
# Match for all valid data-types with or without numberic repeater at the lead
$namelist::valmatch = "$namelist::valrepeat|$namelist::valall";
# Same as above when a match isn't required
$namelist::nrvalmatch = $namelist::valmatch. "||";

#
# This script takes the %main::CCMEXP associative array and stores the keys in
# Lower case to the following lowercase copies. It uses the values passed in
# and sets needed default values based on configuration variables. Then it writes out
# a namelist according to the corresponding resultant associative array.
#

sub new {
#
# Constructor: usage: my $nl = namelist->new( "CCMEXP", "nl.initial", %main::CCMEXP );
#
  my $class = shift;
  my $name  = shift;     # Name of new namelist
  my $file  = shift;     # Namelist filename
  my $NLref = shift;     # Associated namelist hash
  my $printlev = shift;  # Debug print level

  my $self  = {};
  if ( ! defined($name) ) {
    die "ERROR:: name not given to namelist constructor\n";
  }
  if ( ! defined($file) ) {
    die "ERROR:: filename not given to namelist constructor\n";
  }
  if ( ! defined($NLref) || $NLref !~ /HASH/ ) {
    die "ERROR:: reference to namelist associative array not given to namelist constructor\n";
  }
  $self->{'FILENAME'} = $file;                    # Filename of output namelist
  $self->{'NLREF'}    = $NLref;                   # Reference to namelist
  $self->{'NAME'}     = $name;                    # the name of the namelist
  $self->{'printlev'} = $printlev;
  $self->{'VAR'}      = undef;                    # Variable name when parsing
  $self->{'VALUE'}    = undef;                    # Variable value when parsing
  bless( $self, $class );
  return( $self );
}

#============================================================================

sub change {
#
# Make changes to the namelist by prompting at the command-line
#
  my $self = shift;

  my $nm = ref($self) . "::change";
  my $ref = $self->{'NLREF'};
  my $name = $self->{'NAME'};
  print "Here is the $name namelist:\n";
  $self->print;
  my $ans;
  do {
    print "Do you want to add or change any settings? (y/n):";
    $ans = <>; 
    if ( $ans =~ /[yY][Ee]*[sS]*/ ) {
      print "Enter changes to namelist: (key = value) (enter a / to stop data entry)\n";
      print "(Assumes standard F90 namelist format on one line for each keyword pair)\n";
      print "(Be sure and put \' around string values)\n";
      my $bad = "$nm:: Warning:: bad input: enter as: key = value :: key and value should conform to f90 rules\n";
      my $expect = "variable";
      my $line;
      while ( ($expect ne "end") && defined($line = <>) && ($line =~ /./) ) {
        #
        # Parse the line in standard namelist format
        #
        my $err;
        while ( defined($line) && ($line =~ /./) && ($expect ne "end") && ($expect ne "error") ) {
          $err = $self->parse_next( \$line, \$expect );
        }
        if ( $expect eq "error" ) {
          print "WARNING:: $err\n";
          print "$bad"; 
          next; 
        }
      }
      $self->convert_case;
      print "Ok, here is the new $name namelist:\n";
      $self->print;
    }
  } while( $ans =~ /[yY][Ee]*[sS]*/ );
}

#============================================================================

sub checkstring {
#
# Check that a string namelist item is handled correctly
#
  my $self = shift;
  my $item = shift;

  my $EXPNLref = $self->{'NLREF'};
  my %EXPNL = %$EXPNLref;
  my $name = $EXPNL{$item};
  if ( $name !~ /\'(.+)\'/ ) {
    die "$item needs \' around the value";
  }
}

#============================================================================

sub split_namelist_value {
#
# Return a namelist value split up if longer than 70 characters
#
  my $self = shift;
  my $value = shift;

  my $nm = ref($self) . "::split_namelist_value";
  if ( length($value) > 70 ) {
    my $originalvalue = $value;
    my $expect = "value";
    my @list;
    while ( $value =~ /./ ) { 
      my $err = $self->parse_next( \$value, \$expect ) ;
      if ( $expect eq "error" ) { die "$nm::ERROR::$err"; }
      push( @list, $self->{VALUE} );
      $expect = "value";
    }
    my $numberonline = ( 70*($#list+1) ) / length($originalvalue);
    my $i = 0;
    $value = shift @list;
    foreach my $item ( @list ) {
      $value = $value . ", $item";
      if ( ++$i >= $numberonline ) {
        $value = $value . "\n         ";
        $i = 0;
      }
    }
  }
  return( $value );
}

#============================================================================

sub print {
#
# Print the namelist out
#
  my $self = shift;

  my $ref = $self->{'NLREF'};
  my $key;
  my %namelist = %$ref;
  foreach $key ( sort( keys(%namelist) ) ) {
    if ( defined($namelist{$key}) ) {
      my $value = $self->split_namelist_value( $namelist{$key} );
      print " $key = $value\n";
    }
  }
}

#============================================================================

sub Write {
#
# Write out the namelist based on values set in the associative
# arrays
#
  my $self = shift;
  my $append = shift;

  my $ref = $self->{'NLREF'};
  my %namelist = %$ref;
  my $name = $self->{'NAME'};
  my $file = $self->{'FILENAME'};
  if ( defined($append) && $append =~ /Append/i ) {
    open( OUT, ">>$file" ) || die "Can not open namelist file: $file";
  } else {
    if ( -f $file ) { unlink( $file ); }
    open( OUT, ">$file" ) || die "Can not open namelist file: $file";
  }
  print OUT "&$name\n";
  my $key;
  foreach $key ( sort( keys(%namelist) ) ) {
    if ( defined($namelist{$key}) ) {
      my $value = $self->split_namelist_value( $namelist{$key} );
      print OUT " $key\t\t= $value\n";
    }
  }
  print OUT "/\n";
  close( OUT );
}

#============================================================================

sub convert_case {
#
# Convert the case of the keys in the main associative arrays to lowercase.
# Also terminate if there are two keys with the same name but different case.
#
  my $self = shift;

  my $class = ref($self);
  my $nm = "$class\:\:convert_case";

  my $ref = $self->{'NLREF'};
  my $key;
  foreach $key ( keys(%$ref) ) {
    if ( defined($$ref{$key}) ) {
      my $lckey = $key;
      $lckey =~ tr/[A-Z]/[a-z]/;
      my $value = $$ref{$key};
      if ( $key ne $lckey && defined($$ref{$lckey}) ) {
        print "$lckey already defined\n";
        die "$nm: Fix your namelist so that two definitions of $lckey do not exist";
      }
      $$ref{$key}   = undef;
      $$ref{$lckey} = $value;
    }
  }
}

#============================================================================

sub parse {
#
# Parse the namelist from a file
#
  my $self = shift;
  my $filename = shift;

  my $class = ref($self);
  my $nm = "$class\:\:parse";

  my $name = $self->{'NAME'};
  if ( ! defined( $filename ) ) {
    die "ERROR($nm): Namelist filename not passed to parse method\n";
  }
  open( NAMELIST, "<$filename") || die "ERROR($nm): Can not open namelist: $filename\n";
  print "Parse namelist: $name from file: $filename\n" if ($self->{'printlev'}>2);
  #
  # Find the designator for this namelist
  #
  my $found = undef;
  my $line;
  while ( defined($_ = <NAMELIST>) && (/./)  ) {
    if ( /[\$\&]$name(.*?)$/i ) {
      $line = $1;
      $found = 1;
      last;
    }
  }
  if ( ! defined($found) ) {
    print "WARNING($nm): did not find the correct namelist: $name in file: $filename\n" if ($self->{'printlev'}>2);
    return;
  }
  my $expect = "variable";    # First item expected in a namelist is a variable
  goto LINE;
  #
  # Loop over each line in the namelist
  #
NEXT: while ( defined($line = <NAMELIST>) && (/./)  ) {
    #
    # Loop over each item type in each line
    #
LINE: while ( defined($line) && ($line =~ /./) ) {
      my $err = $self->parse_next( \$line, \$expect );
      if ( $expect eq "error" ) { die "$nm::ERROR::$err"; }
      if ( $expect eq "end" ) {
        last LINE;
      } 
    }
    if ( $expect eq "end" ) { last; }
  }
  close( NAMELIST );
  $self->convert_case;
}

#============================================================================

sub setkeypair {
#
# Set the keyword pair
# Loads variable and complementary value into namelist hash from
# VAR and VALUE hash elements.
# Called from parse_next when a variable assignment is complete.
#
  my $self = shift;

  my $nm = ref($self) . "::setkeypair";
  if ( defined( $self->{'VAR'} ) ) {
    my $ref = $self->{'NLREF'};
    my $var = $self->{'VAR'};
    my $val = $self->{'VALUE'};
    if ( ! defined($val) ) {
      die "ERROR::($nm) Value not defined for variable: $var\n";
    }
    $$ref{$var} = $val;
    $self->{'VAR'}   = undef;
    $self->{'VALUE'} = undef;
  }
}

#============================================================================

sub parse_next {
#
# Parse the next item in the line
# parse_next( \$line, \$expect )
# Loads values and variables into VALUE and VAR hash element.
#
# Returns information on an error, if $expect changes to "error"
# otherwise returns "success".
#
  my $self = shift;
  my $line = shift;        # Line read from file
  my $expect = shift;      # Type of item you expect next 
                           # ("variable", "varorvalue", "=", or "value")

  my $class = ref($self);
  my $nm = "$class\:\:parse_next";

  my $err;
  $_ = $$line;
  # Blank line, return and continue
  if ( /^\s*$/ ) {
    $$line = undef;
    return( "success" );
  }
  #
  # Switch based on what type of item you expect
  #
  SWITCH: {
    # Expect a variable
    (($$expect eq "variable") || ($$expect eq "varorvalue")) && do {
       # End-designator (F90 form "/" and non-standard F77 forms (&end) )
       if ( /^\s*\// || /^\s*[\$\&]end/i ) {
         $$line = undef;
         $self->setkeypair;
         $$expect = "end";
         return( "success" );
       }
       # variable
       if ( /^\s*,?\s*($namelist::varmatch)(.*?)$/ ) {
         $$line = $2;
         $$expect = "=";
         $self->setkeypair;
         $self->{'VAR'} = $1;
       } elsif ( $$expect ne "varorvalue" ) {
         $err = "ERROR($nm): expect a variable instead got: $_\n";
         $$expect = "error";
         return( $err );
       # value
       } elsif ( $$expect eq "varorvalue"
         &&   /^\s*([\s,]*)($namelist::nrvalmatch)([\s,]*)(.*?)$/ ) {
         $$line = $4;
         $$expect = "varorvalue";
         my $val = $2;
         my $repeat = undef;
         if ( $val =~ /($namelist::valint)\*($namelist::valall)/ ) {
           $repeat = $1;
           $val = $2;
         }
         $self->{'VALUE'} = $self->{'VALUE'} . ",$val";
         if ( defined($repeat) ) {
           for( my $i = 0; $i < ($repeat-1); $i++ ) {
             $self->{'VALUE'} = $self->{'VALUE'} . ",$val";
           }
         }
         # Comments, only can follow a value
         if ( $$line =~ /^([\s,])*![^!]*$/ ) {
           $$line = undef;
         }
       } else {
         $err = "ERROR($nm): expect a F90 namelist constant or variable instead got: $_\n";
         $$expect = "error";
         return( $err );
       }
       last SWITCH;
    };
    # Expect a "="
    ($$expect eq "=") && do {
       if ( /^\s*=(.*?)$/ ) {
         $$line = $1;
         $$expect = "value";
       } else {
         $err = "ERROR($nm): expect a equal \'=\' sign instead got: $_\n";
         $$expect = "error";
         return( $err );
       }
       last SWITCH;
    };
    # Expect a value
    ($$expect eq "value") && do {
       # value
       if ( /^\s*(${namelist::valmatch})([\s,]*)(.*?)$/ ) {
         $$line = $3;
         $$expect = "varorvalue";
         my $val = $1;
         my $repeat = undef;
         if ( $val =~ /(${namelist::valint})\*($namelist::valall)/ ) {
           $repeat = $1;
           $val = $2;
         }
         $self->{'VALUE'} = "$val";
         if ( defined($repeat) ) {
           for( my $i = 0; $i < ($repeat-1); $i++ ) {
             $self->{'VALUE'} = $self->{'VALUE'} . ",$val";
           }
         }
         # FORTRAN only allows comments after values
         if ( $$line =~ /^\s*![^!]*$/ ) {
           $$line = undef;
         }
       } else {
         $err = "ERROR($nm): expect a F90 constant for a namelist instead got: $_\n";
         $$expect = "error";
         return( $err );
       }
       last SWITCH;
    };
    # default
    $err = "ERROR($nm): Bad type to expect: $$expect\n";
    $$expect = "error";
    return( $err );
  }
}

#============================================================================

# Quoting should be done in the Write method rather
# than when string values are added to the namelist hash.
# But the namelist variable type isn't known in the Write method.
sub quote_string {
    my $str = shift;
    $str =~ s/^\s+//;
    $str =~ s/\s+$//;
    unless ($str =~ /^['"]/) {        #"'
        $str = "\'$str\'";
    }
    return $str;
}

#============================================================================

1   # to make use or require happy
