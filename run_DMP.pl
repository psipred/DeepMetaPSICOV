#!/usr/bin/perl

# Script to run DeepMetaPSICOV (DMP)

# Further modifications SMK June/July 2018:
#    Cleanup dead code
#    modify to start from a user-provided alignment; exclude HHblits-pdb70, HHblits-uniclust30, and jack_hhblits steps
#    skip some steps if output files are present - this may or may not be what a user wants; doing it to save time in case of failures.
#    All errors now produce an exit status of 1 rather than 0. Exit status 0 is only produced on successful runs.
#    cleanup to remove email-related stuff
# TODO:
#    add a --force option to redo all steps even if output files exist
#    add option for user to supply psiblast results(?) - should be just the .chk or .mtx if available
#        interactions with domain parsing
#    add option to cleanup intermediate files; some can be large.
#    check provided db paths before proceeding? check all executables are in place and executable?

use warnings;
use Digest::MD5 qw(md5_hex); 
use BSD::Resource; # SMK for setpriority and setrlimit. Seems to be the only additional Perl dependency, at least on CAMP.
use IO::Handle;
use Text::Wrap;
use File::Basename;
use Getopt::Long qw(GetOptions);

Getopt::Long::Configure qw(gnu_getopt);

setpriority 0, 0, 1;

setrlimit(RLIMIT_AS, RLIM_INFINITY, RLIM_INFINITY);  # (virtual) address space in bytes
setrlimit(RLIMIT_DATA, RLIM_INFINITY, RLIM_INFINITY);  # data size in bytes
setrlimit(RLIMIT_STACK, RLIM_INFINITY, RLIM_INFINITY);  # stack size in bytes (important)

$DMP_VERSION = "1.0.0";
$maindir = dirname(__FILE__);

####################### USER VARIABLES #######################

# Path to python executable to use.
# currently need python 2 or 3, preferably anaconda/miniconda, with pytorch >= 0.3.0 and < 0.4.0
# TODO we will eventually move this to pytorch 1.0 as DMPfold uses the same
$python = "/home/camp/kandats/working/DMP/miniconda2/bin/python";

# Enviromment variables for BLAST and HH-suite
$ENV{'BLASTMAT'} = "$maindir/data/blast";
$ENV{'BLASTDB'} = "$maindir/data/blast";
$ENV{'HHLIB'} = "$maindir/hh-suite-patched";

# Location of HH-suite binaries
$hhbindir = "$maindir/hh-suite-patched/bin";

# Number of threads to use for multithread-capable programs
$psiblast_threads=4;
$hhblits_threads=4;
$psicov_threads=6;
$ccmpred_threads=6;
$freecontact_threads=8;

# Timeouts for PSICOV and CCMpred in seconds. Increase if PSICOV/CCMpred frequently timeout.
# Please also see the man pages for 'timeout'.
$psicov_timeout=86400;
$ccmpred_timeout=86400;

# Database paths:

# HHblits UniClust30 database
$hhblits_db = "/home/camp/kandats/working/hhblitsdb/uniclust30_2017_10/uniclust30_2017_10";

# Legacy BLAST-formatted database for PSI-BLAST (usually UniRef90)
$psiblast_db = "$maindir/data/blast/nr";

################### END OF USER VARIABLES ###################


########################################################
##### Subroutines
#####

sub MailError {
  my $msg = $_[0];

  print STDERR $msg;

  if (open(LOGF, ">>$maindir/DMP.log")) {
      flock(LOGF, 2) or die;
      $date = localtime;
      print LOGF $date," : ERROR $jobid\n";
      close(LOGF);
  }
  
}

sub RunDeepMetaPSICOV
{
    my($ffname, $ofname, $mtx_fname) = @_;  # fasta filename, output filename, psiblast mtx filename
    my($tempdir, $prefix, $suffix) = fileparse($ffname, qr/\.[^.]*/);

    if ($mtx_fname eq ''){
	$mtx_fname = "$tempdir/$prefix.mtx";
    } elsif (! -s $mtx_fname) {
	print STDERR "WARNING: Supplied PSI-BLAST PSSM file $mtx_fname not found or is empty. Running PSI-BLAST to make new PSSM.\n";
	$mtx_fname = "$tempdir/$prefix.mtx";
    }

    if (! -s $mtx_fname){

	# check for existing (finished) psiblast outputs. The check is a bit trickier than all subsequent output file checks
	$run_psiblast=1;

	if (-s "$tempdir/$prefix.blast" && -s "$tempdir/$prefix.chk")
	{
	    # compare cols 1-3 of the last 7 lines of the blast file to a sample prepared earlier from a successful run
	    @cmd=("bash", "-c", "diff -q $maindir/data/psiblast-tail-check <(tail -n 7 $tempdir/$prefix.blast | cut -c 1-3) >/dev/null");

	    if(system(@cmd) == 0)
	    {
		$run_psiblast=0;
		print "PSI-BLAST outputs found at $tempdir/$prefix.chk and $tempdir/$prefix.blast. Will use these to make PSSM.\n";
	    }
	}

	if ($run_psiblast)
	{
	    print "Running PSI-BLAST...\n";
	    if (system("$bindir/blastpgp -a $psiblast_threads -b 0 -v 2000 -j 3 -h 0.001 -e 0.001 -d $psiblast_db -i $ffname -C $tempdir/$prefix.chk > $tempdir/$prefix.blast") != 0)
	    {
		&MailError("DMP ERROR 01 - please report error to psipred\@cs.ucl.ac.uk\n");
		exit 1;
	    }
	}
	open(MMOUT, ">$tempdir/$prefix.pn") or die;
	print MMOUT "$prefix.chk";
	close(MMOUT);
	open(MMOUT, ">$tempdir/$prefix.sn") or die;
	print MMOUT "$ffname";
	close(MMOUT);

	print "Preparing PSSM from PSI-BLAST result...\n";
	if (system("$bindir/makemat -P $tempdir/$prefix > /dev/null") != 0)
	{
	    &MailError("DMP ERROR 02 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    }
    else
    {
	print "Using existing PSI-BLAST PSSM at $mtx_fname.\n";
    }

    if ($no_aln) { # $no_aln is set outside this function
	print "Running HHblits...\n";
	if (system("$hhbindir/hhblits -i $ffname -n 3 -e 0.001 -d $hhblits_db -cpu $hhblits_threads -oa3m $tempdir/$prefix.a3m -diff inf -cov 50 -id 99 > $tempdir/$prefix.hhblog 2>&1") != 0)
	{
	    &MailError("DMP ERROR 03 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
	
	if (system("grep -v '^>' $tempdir/$prefix.a3m | sed 's/[a-z]//g' > $tempdir/$prefix.hhbaln") != 0)
	{
	    &MailError("DMP ERROR 04 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}

	$seqlen = `head -1 $tempdir/$prefix.hhbaln | wc -c`;
	$naln_hhblits = `cat $tempdir/$prefix.hhbaln | wc -l`;
	
	if ($naln_hhblits < 10 * $seqlen)
	{
	    # Try to find more sequences via jackhmmer/UNIREF100
	    print "HHblits alignment has fewer than 10L sequences. Running jack_hhblits to find more...\n";
	    if (system("$bindir/jack_hhblits $prefix $bindir $tempdir $jhhb_db > $tempdir/$prefix.jacklog 2> $tempdir/$prefix.jackerr") != 0)
	    {
		&MailError("DMP ERROR 05 - please report error to psipred\@cs.ucl.ac.uk\n");
		exit 1;
	    }
	    
	    $naln_jack = `cat $tempdir/$prefix.jackaln | wc -l`;
	}
	else
	{
	    print STDERR "WARNING: alignment has fewer than 10L sequences.\n";
	    $naln_jack = 0;
	}
	
	if ($naln_jack > $naln_hhblits)
	{
	    system("cp -f $tempdir/$prefix.jackaln $tempdir/$prefix.aln")
	}
	else
	{
	    system("cp -f $tempdir/$prefix.hhbaln $tempdir/$prefix.aln")
	}
    } # if ($no_aln)
    
    if (! -s "$tempdir/$prefix.ss") {
	print "Running PSIPRED pass 1...\n";
	if (system("$bindir/psipred $tempdir/$prefix.mtx $maindir/data/psipred4/weights.dat $maindir/data/psipred/weights.dat2 $maindir/data/psipred/weights.dat3 > $tempdir/$prefix.ss") != 0)
	{
	    &MailError("DMP ERROR 06 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    }

    if (! -s "$tempdir/$prefix.ss2") {
	print "Running PSIPRED pass 2...\n";
	if (system("$bindir/psipass2 $maindir/data/psipred4/weights_p2.dat 1 1.0 1.0 $tempdir/$prefix.ss2 $tempdir/$prefix.ss > /dev/null") != 0)
	{
	    &MailError("DMP ERROR 07 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    }

    if (! -s "$tempdir/$prefix.solv") {
	print "Running SOLVPRED...\n";
	if (system("$bindir/solvpred $tempdir/$prefix.mtx $maindir/data/weights_solv.dat > $tempdir/$prefix.solv") != 0)
	{
	    &MailError("DMP ERROR 08 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    }

    if ((! -s "$tempdir/$prefix.colstats") || (! -s "$tempdir/$prefix.pairstats") ) {
	print "Calculating alignment statistics...\n";
	if (system("$bindir/alnstats $tempdir/$prefix.aln $tempdir/$prefix.colstats $tempdir/$prefix.pairstats > /dev/null") != 0)
	{
	    &MailError("DMP ERROR 09 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    }
    
    $naln = `cat $tempdir/$prefix.aln | wc -l`;

    if ($naln >= 5)
    {
	# Run contact predictors if they haven't timed out and the output file doesn't exist

	if ((! -e "$tempdir/$prefix.psicov.timeout") && (! -s "$tempdir/$prefix.psicov")){	 
	    print "Running PSICOV...\n";

	    $psicovret = system("timeout $psicov_timeout $bindir/psicov -z $psicov_threads -o -d 0.03 $tempdir/$prefix.aln > $tempdir/$prefix.psicov 2>&1");
	    
	    # handle timeouts and failures differently
	    if ($psicovret == 124){
		print "PSICOV timed out after $psicov_timeout seconds. This is not an error, but you can try increasing this time limit.\n";
		open(PT, ">$tempdir/$prefix.psicov.timeout");
		close PT;
	    }
	    if ($psicovret != 0 && $psicovret >> 8 != 124)
	    {
		&MailError("DMP ERROR 10 ($psicovret) - please report error to psipred\@cs.ucl.ac.uk\n");
		exit 1;
	    }
	}

	if ((! -e "$tempdir/$prefix.ccmpred.timeout") && (! -s "$tempdir/$prefix.ccmpred")){	 
	    print "Running CCMpred...\n";
	    $ccmret = system("timeout $ccmpred_timeout $bindir/ccmpred -t $ccmpred_threads $tempdir/$prefix.aln $tempdir/$prefix.ccmpred > /dev/null 2>&1");
	    if ($ccmret == 124){
		print "CCMpred timed out after $ccmpred_timeout seconds. This is not an error, but you can try increasing this time limit.\n";
		open(CT, ">$prefix.ccmpred.timeout");
		close CT;
	    }
	    
	    if ($ccmret != 0 && $ccmret >> 8 != 124)
	    {
		&MailError("DMP ERROR 11 ($ccmret) - please report error to psipred\@cs.ucl.ac.uk\n");
		exit 1;
	    }
	}

	print "Running FreeContact...\n";
	if (system("$bindir/freecontact -a $freecontact_threads < $tempdir/$prefix.aln > $tempdir/$prefix.evfold") != 0)
	{
	    &MailError("DMP ERROR 12 - please report error to psipred\@cs.ucl.ac.uk\n");
	    exit 1;
	}
    } 

    if (system("touch $tempdir/$prefix.psicov $tempdir/$prefix.evfold $tempdir/$prefix.ccmpred") != 0)
    {
	&MailError("DMP ERROR 13 - please report error to psipred\@cs.ucl.ac.uk\n");
	exit 1;
    }
    
    print "Preparing DMP inputs (1/2)...\n";   
    if (system("$bindir/deepmetapsicov_makepredmap $tempdir/$prefix.colstats $tempdir/$prefix.pairstats $tempdir/$prefix.psicov $tempdir/$prefix.evfold $tempdir/$prefix.ccmpred $tempdir/$prefix.ss2 $tempdir/$prefix.solv $tempdir/$prefix.deepmetapsicov.map $tempdir/$prefix.deepmetapsicov.fix > /dev/null 2>&1") != 0)
    {
	&MailError("DMP ERROR 14 - please report error to psipred\@cs.ucl.ac.uk\n");
	exit 1;
    }
    print "Preparing DMP inputs (2/2)...\n";       
    if (system("$bindir/cov21stats $tempdir/$prefix.aln $tempdir/$prefix.deepmetapsicov.21c > /dev/null 2>&1") != 0)
    {
	&MailError("DMP ERROR 15 - please report error to psipred\@cs.ucl.ac.uk\n");
	exit 1;
    }
    print "Predicting DMP contacts...\n";   
    if (system("$python $deepdir/pytorch_metacov_pred.py $tempdir/$prefix.deepmetapsicov.21c $tempdir/$prefix.deepmetapsicov.map > $ofname 2>&1") != 0)
    {
	&MailError("DMP ERROR 16 - please report error to psipred\@cs.ucl.ac.uk\n");
	exit 1;
    }
}    

###################
##      MAIN
###################

$longusage = <<'END';

required arguments:
  -i fasta_file, --input-fasta fasta_file
                        FASTA-formatted file with target sequence.

optional arguments:
  -h, --help            show this help message and exit.

  -a alignment_file, --input-aln alignment_file
                        Path to input alignment in PSICOV format.
                        If supplied, alignment generation steps will be skipped unless --force is also specified.

  -o contact_file, --output-file contact_file
                        Filename for output contacts in CASP format.
                        The default is to create a file with the extension '.deepmetapsicov.final.con' and the basename of the input file.

END
# options not implemented yet. Some are not implemented because it's unclear how they'd work under all conditions
#  -p, --parse-domains   Enable automatic domain parsing. Requires HHblits-formatted PDB70 database.
#                        By default, no domain parsing is attempted.
#  -j, --jack-hhblits    Enable the use of jack_hhblits for more intensive sequence searching.
#                        Requires additional programs and an indexed large sequence database e.g. UniRef100.
#                        See README for more details. Disabled by default.

#  -d, --deeper          Use double-depth DMP model with 36 residual blocks.
#                        By default, a model with 18 residual blocks is used.

#  --force               Force generation of ALL intermediate files, even if these already exist.
#                        By default, existing intermediate files with the correct path and filename are reused.
#                        Contact prediction is always attempted regardless of --force.

#  -m mtx_file, --psiblast-mtx mtx_file
#                        Path to legacy-format PSI-BLAST PSSM file generated by 'makemat -P' for target sequence (extension .mtx).
#                        If supplied, PSI-BLAST will not be run.

#  -c, --cleanup         Remove all intermediate files except DMP predictions. By default, all intermediate files are retained.

#  --dry-run             Show all commands run in the DMP pipeline but do not execute them. This is useful for troubleshooting.
#                        The default is to not do a dry run.

#  -r, --casp-format     Add header and footer lines to output so that it is fully compliant with the CASP RR format.

sub splash{
    print "DeepMetaPSICOV (DMP) v${DMP_VERSION} by Shaun M. Kandathil and David T. Jones\n";
}

sub usage{
    print STDERR 'Usage: perl ' . __FILE__ . " [-h] -i fasta_file [-a alignment_file] [-o contact_file]\n";
}

sub long_usage{
    usage();
    print $longusage;
}

###########
# Important paths
$bindir = "$maindir/bin";
$deepdir = "$maindir/deepmetapsicov_consens";

# default values for command-line options
$fasta_path = '';
$aln_path = '';
$output_file = '';
$help = 0;
#$parse_domains = 0;
#$go_deeper = 0;
#$jack_hhblits = 0;

# TODO other options not yet implemented:
$force = 0;
$cleanup = 0;
# $dry_run = 1; # TODO change to 0 before release
 $psiblast_mtx = '';
# $casp_format = 0;

# Parse command-line options

    GetOptions('i|input-fasta=s' => \$fasta_path,
	       'h|help' => \$help,
	       'a|input-aln=s' => \$aln_path,
	       'o|output-file=s' => \$output_file
#	       'p|psiblast-mtx=s' => \$psiblast_mtx,
#	       'force' => \$force,
#	       'dry-run' => \$dry_run,
#	       'c|cleanup' => \$cleanup,
#	       'r|casp-format' => \$casp_format,
#	       'p|parse-domains' => \$parse_domains,
#	       'd|deeper' => \$go_deeper,
#	       'j|jack-hhblits' => \$jack_hhblits
    );

# check options & combinations
if ($help){
    splash();
    long_usage();
    exit 0;
}

if ($fasta_path eq '') {
    print STDERR "Error: You must supply a path to a target sequence file in FASTA format using -i.\n";
    usage();
    exit 1;
} elsif (! -s $fasta_path) {
    &MailError("Error: supplied FASTA file $fasta_path does not exist or is empty.");
}

# get $jobid and $tempdir
($jobid, $tempdir, $suffix) = fileparse($fasta_path, qr/\.[^.]*/);

$no_aln = 0;
if ($aln_path eq ''){
    $no_aln = 1;
} elsif (! -s $aln_path){
    &MailError("Error: supplied alignment $aln_path does not exist or is empty.");
    exit 1;
}

if ($no_aln){
    $aln_path = "$tempdir/$prefix.aln";
}

## TODO if aln is also given, need to check that the number of columns is the same as sequence length, and that format is correct

open(SEQFILE, "$fasta_path") or die;
$seq = do { local $/; <SEQFILE> };
$seq =~ s/>.*\n//;
$seq =~ s/[^A-za-z]//g;
close SEQFILE;

$nres = length($seq);

# if ($go_deeper) {
#     $deepdir="$maindir/deeper-DMP"
# }

# undef $suffix;

if ($output_file eq ''){
    $output_file = "$tempdir/$jobid.deepmetapsicov.con"
}

if (open(LOGF, ">>$maindir/DMP.log")) {
    flock(LOGF, 2) or die;
    $date = localtime;
    print LOGF $date," : RUNNING $jobid\n";
    close(LOGF);
}


#print "** Running DMP on target sequence...\n";
&RunDeepMetaPSICOV($fasta_path, $output_file, $psiblast_mtx, $aln_path);

# Additional stuff for CASP format compliance, emailing results etc

# $message = "PFRMAT RR\nTARGET " . $subject . "\nAUTHOR DMP\nMODEL 1\n";
# $subjectline = $subject . " DeepMetaPSICOV Results";
# open(MAILOUT, "| /usr/bin/Mail -s '$subjectline' $email") or die;
# print MAILOUT $message;

# $Text::Wrap::columns = 50;
# open(SEQFILE, "$fasta_file");
# <SEQFILE>;
# print MAILOUT Text::Wrap::fill( '', '', join '', <SEQFILE>);
# close SEQFILE;
# print MAILOUT "\n";

# $results = `cat $output_file | sort -g -r -k 5 | head -10000`;
# print MAILOUT $results, "END\n";

# close(MAILOUT);

if (open(LOGF, ">>$maindir/DMP.log")) {
    flock(LOGF, 2) or die;
    $date = localtime;
    print LOGF $date," : COMPLETED $jobid\n";
    close(LOGF);
}

# TODO cleanup
# rm -rf *.ss* *.hh* *.pdbhhr *.aln *.solv *.psicov *.ccmpred *.evfold *.mtx *.blast *.chk *.pairstats *.colstats *.21c *.map __pycache__ *.timeout *.pn *.sn

print "Done. Predicted contacts can be found in $output_file \n";
exit 0;
