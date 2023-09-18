#!/usr/bin/env perl 
#===============================================================================
=pod

=head2

         FILE: translate_fasta_headers.pl

        USAGE: ./translate_fasta_headers.pl [options] <file>

               # From long to short labels:
               ./translate_fasta_headers.pl --out=short.fas long.fas

               # And back, using a translation table:
               ./translate_fasta_headers.pl --tabfile=short.fas.translation.tab short.fas

               # Slightly shorter/different version:
               ./translate_fasta_headers.pl long.fas > short.fas
               ./translate_fasta_headers.pl -t long.fas.translation.tab short.fas > back.fas

               # Use your own prefix:
               ./translate_fasta_headers.pl --prefix='Own_' long.fas > short.fas


  DESCRIPTION: Replace fasta headers with headers taken from tab delimited file.
               If no tab file is given, the (potentially long) fasta headers are replaced 
               by short labels "Seq_1", "Seq_2", etc, and the short and original headers
               are printed to a translation file.

      OPTIONS: -t, --tabfile=<filename> -- Specify tab-separated translation file with
                                           unique "short" labels to the left, and "long"
                                           names to the right. Translation will be from
                                           left to right.

               -i, --in=<filename>      -- Specify name of fasta file. Can be skipped as
                                           script reads files from STDIN.

               -o, --out=<filename>     -- Specify output file for the fasta sequences.
                                           Note: If --out=<filename> is specified, the
                                           translation file will be named
                                           <filename>.translation.tab. This simplifies
                                           back translation. If '--out' is not used,
                                           the translation file will be named after
                                           the infile!

               -n, --notab              -- Do not create a translation file.

               -p, --prefix=<string>    -- Prefix for short label. Defaults to 'Seq_'.

               -f, --forceorder         -- [NOT IMPLEMENTED] translate in order of
                                           appearance in the fasta file, and use
                                           the same order as in the tabfile - without
                                           rigid checking of the names! This allows
                                           non-unique labels in the left column.

                -h, --help              -- Show this help text and quit.

 REQUIREMENTS: ---

         BUGS: ---

        NOTES: ---

       AUTHOR: Johan.Nylander\@nbis.se 

      COMPANY: NBIS/NRM

      VERSION: 1.0.1

      CREATED: 03/13/2013 01:52:28 PM

     REVISION: 01/10/2018 12:59:31 PM

         TODO: Handle non-unique values in the left tabfile column
               (can't use hash):
               Test if values in translation table are unique. If so,
               use read_tabfile. If not, read into two arrays, check
               for same lengths, and then use an iterator while reading
               the infile. Warn if number of sequences doesn't match 
               number of entries in the tab file. Plus give a warning
               that labels where not unique. Use the array approach
               when '--forceorder'.

      LICENSE: Copyright (c) 2019-2020 Johan Nylander
               
               Permission is hereby granted, free of charge, to any person
               obtaining a copy of this software and associated documentation
               files (the "Software"), to deal in the Software without
               restriction, including without limitation the rights to use,
               copy, modify, merge, publish, distribute, sublicense, and/or
               sell copies of the Software, and to permit persons to whom the
               Software is furnished to do so, subject to the following
               conditions:
               
               The above copyright notice and this permission notice shall be
               included in all copies or substantial portions of the Software.
               
               THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
               EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
               OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
               NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
               HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
               WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
               FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
               OTHER DEALINGS IN THE SOFTWARE.
 

=cut
#===============================================================================

use strict;
use warnings;
use Getopt::Long;

## Globals
my $prefix       = 'Seq_'; # Prefix for short names
my $tabfile      = q{};
my $in           = q{};
my $out          = q{};
my $notab        = q{};
my $help         = q{};
my %header_hash  = ();     # Key: short, value: long.
my @short_array  = ();     # Array with short headers.
my @long_array   = ();     # Array with long headers.
my $PRINT;                 # Print file handle. Using the typeglob notation below
my $IN;                    # in order to use STDOUT as a variable.
#my $forceorder   = q{};

## If no arguments
exec("perldoc", $0) unless (@ARGV);

## Get args
my $r = GetOptions("tabfile|translationfile=s" => \$tabfile,
                   "out=s"                     => \$out,
                   "notab"                     => \$notab,
                   "in=s"                      => \$in,
                   "prefix=s"                  => \$prefix,
                   #"forceorder"                => \$forceorder,
                   "help"                      => sub { exec("perldoc", $0); exit(0); },
                  );

## If tabfile, read it assuming "short" labels in the left column
if ($tabfile) {
    %header_hash = read_tabfile($tabfile);
}

## If outputfile
if ($out) {
    open ($PRINT, '>', $out) or die "$0 : Failed to open output file $out : $!\n\n";
}
else {
    $PRINT = *STDOUT; # Using the typeglob notation in order to use STDOUT as a variable
}

## If infile
if ($in) {
    read_infile($in, $prefix);
}
else {
    while (my $infile = shift(@ARGV)) {
        read_infile($infile, $prefix);
    }
}

exit(0);


#===  FUNCTION  ================================================================
#         NAME:  read_infile
#      VERSION:  01/10/2018 01:03:55 PM
#  DESCRIPTION:  Reads a tab separated file and returns a hash. Expects all values
#                in left column ("short") to be unique!
#   PARAMETERS:  filename, prefix
#      RETURNS:  hash: key:short, value:long
#         TODO: 
#===============================================================================
sub read_infile {

    my ($file, $shortlabel) = (@_);
    my @in_headers_array    = ();
    my $counter             = 1;
    my $OUTTAB;
    my $outtabfile;

    ## Name translation.file from out if available. Saves time with back translation...
    if ($out) {
        $outtabfile = $out . '.translation.tab';
    }
    else {
        $outtabfile = $file . '.translation.tab';
    }

    open my $IN, "<", $file or die "Could not open $file for reading : $! \n";

    ## Open tab translation outfile or not
    if ($tabfile) {
        # not open
    }
    elsif ($notab) {
        # not open
    }
    else {
        if ( -e $outtabfile ) {
            die "\nCovardly refusing to overwrite existing $outtabfile . Rename files of use \'--notab\'.\n\n";
        }
        else {
            open $OUTTAB, ">", $outtabfile or die "Could not open $outtabfile for writing : $! \n";
        }
    }

    ## Read fasta file
    while(<$IN>) {
        my $line = $_;
        chomp($line);
        if ($line =~ /^\s*>/) {
            $line =~ s/^\s*>//;
            $line = trim_white_space($line);
            if ($tabfile) {
                if (exists $header_hash{$line}) {
                    print $PRINT '>', $header_hash{$line}, "\n";
                }
                else {
                    die "\nLabel \"$line\" in $file not found in $tabfile.\n";
                }
            }
            else {
                my $header = $shortlabel . $counter;
                print $PRINT '>', $header, "\n";
                print($OUTTAB $header, "\t", $line, "\n") unless $notab;
                $counter++;
            }
        }
        else {
            print $PRINT $line, "\n";
        }
    }
    close($IN);
    #close($OUTTAB); # Let garbage collector do some work

} # end of read_infile


#===  FUNCTION  ================================================================
#         NAME:  read_tabfile
#      VERSION:  03/14/2013 09:42:47 PM
#  DESCRIPTION:  read tabfile and return hash
#   PARAMETERS:  filename
#      RETURNS:  hash (key:short, value:long)
#         TODO:  ????
#===============================================================================
sub read_tabfile {

    my ($file) = @_; 

    my %hash = (); # Key: short, value: long

    open my $TAB, "<", $file or die "Could not open $file for reading : $! \n";
    while(<$TAB>) {
        chomp;
        next if (/^\s*$/);
        my ($short, $long) = split /\t/, $_;
        $short = trim_white_space($short);
        $long = trim_white_space($long);
        if ($hash{$short}++) {
            die "Warning, the values in the left column of the tabfile must be unique (I found more than one $short)\n";
        }
        else {
            $hash{$short} = $long;
        }
    }
    close($TAB);

    return(%hash);

} # end of read_tabfile


#===  FUNCTION  ================================================================
#         NAME:  trim_white_space
#      VERSION:  03/14/2013 09:44:02 PM
#  DESCRIPTION:  trim white space from both ends of string 
#   PARAMETERS:  string
#      RETURNS:  string
#         TODO:  ????
#===============================================================================
sub trim_white_space {

	my ($a) = @_;

	$a =~ s/^\s+//;
	$a =~ s/\s+$//;

	return($a);

} # end of trim_white_space

