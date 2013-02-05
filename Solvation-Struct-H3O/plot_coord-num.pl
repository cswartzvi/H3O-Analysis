#! /usr/bin/perl -w
#
#

use strict;
use warnings;

use File::Basename;
use Term::ReadLine;

# Check command line arguments, open filename
my $filename = undef;
if ($#ARGV != 0){
   print (" ERROR: Filename missing from command line.");
} 
else {
   $filename = $ARGV[0];
   open INFILE, "<", "$filename" or die " ERROR: Can't Open $filename!!";
}

while (<INFILE>) {

      my @temp = split /\s+\n?/, $_; 

      if ( $temp[2] =~ /^[^PT]/) {
         print (" $temp[1] $temp[3] $temp[4] \n"); 
      }
}


#Takes User prompt and reads in variable
#sub prompt {
#   #print("@_:");
#
#   my $Term = new Term::ReadLine 'Installer';
#
#   my $val = $Term->readline("@_:");
#   #my $val = <STDIN> ;
#   chomp($val);
#
#   return $val;
#}
