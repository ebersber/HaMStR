#!/usr/bin/perl
use strict;
use Getopt::Long;
use lib '/Users/ingo/src/hamstr/hamstr.v12/bin/lib';
use lib '/Users/ingo/src/hamstr/hamstr.v12/bin/lib/Bio';
use Bio::SearchIO;
use Bio::Search::Hit::BlastHit;
use run_genewise_hamstr;

# PROGRAMNAME: hamstrsearch_local.pl

# Copyright (C) 2009 INGO EBERSBERGER, ingo.ebersberger@univie.ac.at
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 3 of the License
# or any later version.

# This program is distributed in the hope that it will be useful
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program; If not, see http://www.gnu.org/licenses
 
# PROGRAM DESCRIPTION: 

# DATE: Wed Dec 19 10:41:09 CEST 2007
# Date last modified: 
	##23. 07. 2010: found a bug in the extraction of the
	## hmm hit sequence from the sequnence_file. A end-of-line char was missing.
	
	##09.08.2010: added the option to choose the new blastp program from ncbi. Just comment
	##out line 45 in the script and uncomment line 46. Note, in order to make this work I have
	##to slightly modify the blast output since otherwise it will not be parsed by the Bioperl
	##Blast parser. Currently this is a pretty dirty $sedprog hack. It will also take care of removin
        ##the string lcl| that is added in some instances to the hit id in the blast output.
        
        ## I added the option that one can now provide a comma-separated string of phmm names as an
        ## argument for the option -hmm 

	## 08.03.2011: 
	## 1) BUG-FIX: Hamstr will now remove automatically newlines from the input sequence file
	## 2) BUG-FIX: The sequence header remains now the same whether or not the flag -representative
	## has been chosen.
	
	## 10.04.2011
	## 1) added some information to the log file.

	## 20.05.2011
	## 1) BUG-FIX: The grep for the EST sequence in the sub-routine predictORF received also a hit when
	## the search pattern was only a substring of the EST sequence identifier. In some cases the wrong EST
	## was then used to predict the ORF. This has been fixed. 

	## 30.05.2011
	## 1) Extension: a command line option -longhead has been added. The user can now specify that the
	## full sequence id including whitespaces will considered throughout the hamstr search. Note, the
	## whitespaces will be replaced by the string specified in the variabel $idsep.
	## 2) Modification from the bug fix from 20.05.2011. In the grep for the original EST it is no longer
	## necessary that the search string and the EST id are identical over their entire length. Instead the
	## search string may be a prefix of the EST id ending with a whitespace.
	
	## 27.06.2011
	## 1) Extension: I added the option to run a true reciprocal best hit search. Only the best hit from the
	## hmmer search is used to check for reciprocity.

	## 06.12.2011
	## 1) Extension: I added the option -hit_limit to set the number of hmmsearch hits that HaMStR uses for
	## the re-blast.

	## 10.02.2012
	## 1) Extension: I added checks for the appropriate hmmsearch version (HMMER 3) and for genewise and
	## its environmental variable WISECONFIGDIR.
	## 2) Bug fix in the -rbh option.

	## 11.09.2012
	## 1) Bug fix: -hitlimit, even if not set explicitely has been invoked resulting in a more stringent
	## behaviour of HaMStR. This has been fixed resulting in longer run-times.
	## 18.12.2012
	## 1) Bug fix: There was a bug in the CDS extraction for reverse complemented 
	## sequences. A new line was moved to the beginning of the sequence
	## leading to index errors.

    ## 18.12.2013
    ## 1) Bug fix: I have now adapted the script such that it no longer requires the default directory structure
######################## start main #############################
my $version = "hamstrsearch_server-hmmer3.v12.pl\n";
print "$version\n";
### EDIT THE FOLLOWING LINES TO CUSTOMIZE YOUR SCRIPT
my $prog = 'hmmsearch'; #program for the hmm search
my $path =  '/Users/ingo/src/hamstr/hamstr.v12';
my $sedprog = 'sed';
my $grepprog = 'grep';
my $alignmentprog = 'clustalw2';
########## EDIT THE FOLLOWING TWO LINES TO CHOOSE YOUR BLAST PROGRAM ##########
my $blast_prog = 'blastall';
#my $blast_prog = 'blastp';
###############################################################################

my $hmmpath = "$path/core_orthologs"; #path where the hmms are located
my $blastpath = "$path/blast_dir"; #path to the blast-dbs
my $tmpdir = 'tmp';
my $eval = 1; # eval cutoff for the hmm search
my $idsep = '__'; #character used to replace whitespaces in the sequence header with (flag -longhead)

my $hmm_dir = 'hmm_dir';
my $fa_dir  = 'fa_dir';
##############################
my $pid = $$;
my $help;
my $seq2store_file='';
my $cds2store_file='';
my $hmm;
my @hmms;
my $fa;
my $fafile;
my @seqs2store;
my @cds2store;
my $dbpath = '.';	
my $ep2eg;
my $estfile;
my $aln;
my $idfile;
my $taxon_check = 0;
my $hmmset;
my $show_coreortholog_sets;
my $hmmsearch_dir;
my $dbfile = ''; # the file hmmsearch is run against
my $dbfile_short;
my $taxon_file;
my $refspec_string;
my @refspec = qw();
my @primer_taxa;
my $refspec_name = '';
my $taxon_global;
my $fileobj;
my $fa_dir_neu = '';
my $gwrefprot;
my $seqtype;
my $align;
my $rep;
my $estflag;
my $proteinflag;
my $refseq;
my $strict;
my $relaxed;
my $refspec_final = '';
my $concat;
my $seqs2store_file;
my $append;
my $longhead;
my $check = 1;
my @log = qw();
my $filter = 'T';
my $bhh;
my $hitlimit;

#####determine the hostname#######
push @log, "VERSION\n\t$version\n";
my $hostname = `hostname`;
chomp $hostname;
print "hostname is $hostname\n";
push @log, "HOSTNAME\n\t$hostname\n";
#################################
if (@ARGV==0) {
	$help = 1;
}
## help message
my $helpmessage = "
This program is freely distributed under a GPL. See -version for more info
Copyright (c) GRL limited: portions of the code are from separate copyrights

\nUSAGE: hamstrsearch_local.pl -sequence_file=<> -hmmset=<> -taxon=<>  -refspec=<> [-est|-protein] [-hmm=<>] [-representative] [-h]

OPTIONS:

REQUIRED
-sequence_file 
		path and name of the file containing the sequences hmmer is run against.
-est
		set this flag if you are searching in ESTs. Note, if neither the -est nor the -protein flag is set, HaMStR will
		guess the sequence type.
-protein
		set this flag if you are searching in protein sequences. Note, if neither the -est nor the -protein flag is set, HaMStR will
		guess the sequence type.
-hmmset
		specifies the name of the core-ortholog set.
		The program will look for the files in the default directory 'core-orthologs' unless you specify
		a different path via the option -hmmpath.
-taxon
		You need to specify a default taxon name from which your ESTs or protein sequences are derived.
-refspec
		sets the reference species. Note, it has to be a species that contributed sequences 
		to the hmms you are using. NO DEFAULT IS SET! For a list of possible reference
		taxa you can have a look at the speclist.txt file in the default core-ortholog sets
		that come with this distribution. Please use the abreviations in this list. If you choose
		to use core-orthologs where not every taxon is represented in all core-orthologs, you
		can provide a comma-separated list with the preferred refspec first. The lower-ranking 
		reference species will only be used if a certain gene is not present in the preferred 
		refspecies due to alternative paths in the transitive closure to define the core-orthologs.
		CURRENTLY NO CHECK IS IMPLEMENTED!
		NOTE: A BLAST-DB FOR THE REFERENCE SPECIES IS REQUIRED!

USING NON-DEFAULT PATHS
-blastpath
		Lets you specify the absolute or relative path to the blast databases. DEFAULT: $blastpath
-hmmpath
		Lets you specify the absolute or relative path to the core ortholog set. DEFAULT: $hmmpath
        
ADDITIONAL OPTIONS
-append
		set this flag if the output should be appended to the files *.out and *_cds.out. This becomes relevant when running
		hamstrsearch with individual hmms and you want to combine the results.
-eval_limit=<>
		This options allows to set the e-value cut-off for the HMM search.
		DEFAULT: 1
-filter=<T|F>
		set this flag to F if the re-blast should be performed without low-complexity filtering. Default is T.
-hit_limit=<>
		By default, HaMStR will re-blast all hmmsearch hits against the reference proteome. Reduce the number
		of hits for reblast with this option.
-hmm
		option to provide only a single hmm to be used for the search. 
		Note, this file has to end with .hmm 
-show_hmmsets
		setting this flag will list all available core ortholog sets in the specified path. Can be combined with -hmmpath. 
-rbh
		set this flag if you want to use a reciprocal best hit criterion. Only the highest scoring
		hit from the hmmer search will be used for re-blast.
-relaxed
		set this flag if the reciprocity criterion is fulfilled when the re-blast against
		any of the primer taxa was successfull. Note that setting this flag will substantially decrease the
		stringency of the ortholog assignment with the consequence of an increased number of false positives.
-representative
		From all sequences that fulfill the reciprocity criterion the one showing the highest similarity to the
		core ortholog sequence in the reference species is identified and selected as representative. 
-strict
		set this flag if the reciprocity criterion is only fulfilled when the re-blast against
		all primer taxa was successfull
-longhead
		set this flag in the case your sequence identifier contain whitespaces and you whish to keep
		the entire sequence identifier throughout your analysis. HaMStR will then replace the whitespaces with 
		a '__'. If this flag is not set, HaMStR will truncate the sequence
        Identifier at the first whitespace, however if and only if the sequence identifier then remain unique.
		\n\n";
GetOptions ("h"        => \$help,
            "hmm=s"    => \$hmm,
            "est"    => \$estflag,
            "protein"=> \$proteinflag,
            "sequence_file=s" => \$dbfile,
            "fasta_file=s" => \$fafile,
            "hmmset=s" => \$hmmset,
            "hmmpath=s" => \$hmmpath,
            "taxon_file=s" => \$taxon_file,
            "taxon=s"  => \$taxon_global,
            "eval_limit=s" => \$eval,
            "refspec=s" => \$refspec_string,
            "estfile=s" => \$estfile,
            "representative" => \$rep,
	        "strict" => \$strict,
	        "relaxed" => \$relaxed,
	        "append" => \$append,
	        "filter=s" => \$filter,
	        "longhead" => \$longhead,
	        "rbh" => \$bhh,
	        "hit_limit=s" => \$hitlimit,
            "blastpath=s" => \$blastpath,
            "show_hmmsets" => \$show_coreortholog_sets);

if ($help) {
  print $helpmessage;
  exit;
}

## 1) check if all information is available to run HaMStR
($check, @log) = &checkInput();
if ($check == 0) {
  print "$helpmessage";
  print "#######################################\n";
  print "There was an error running $version\n";
  print join "\n", @log;
  exit;
}
else {
  open (OUT, ">hamstrsearch.log") or die "could not open logfile\n";
  print OUT join "\n", @log;
  close OUT;
}
### read in of the core-ortholog sequences
my $co_seqs = parseSeqfile("$fafile");

## 2) loop through the hmms
## process each hmm file separately
for (my $i = 0; $i < @hmms; $i++) {
  $fileobj = undef;
  my @seqs = qw();
  my @newseqs = qw();## var to contain the sequences to be added to the orthologous cluster
  my @newcds = qw();
  my $hmm = $hmms[$i];
  my $hmmout = $hmm;
  $hmmout =~ s/\.hmm/\.out/;
  ## 3) run the hmm search
  if (!(-e "$hmmsearch_dir/$hmmout")) {
    print "now running $prog using $hmm\n";
#	print "$prog $hmm_dir/$hmm $dbfile >$hmmsearch_dir/$hmmout";
    !`$prog $hmm_dir/$hmm $dbfile >$hmmsearch_dir/$hmmout` or die "Problem running hmmsearch\n";
  }
  else {
    print "an hmmresult $hmmout already exists. Using this one!\n";
  }
  
  ## 4) process the hmm search result
  my $hitcount = 0;
  ## 4a) loop through the individual results
  ## now the modified version for hmmer3 comes
  my $hitlimit_local = $hitlimit;
  my ($query_name, @results) = parseHmmer3($hmmout, $hmmsearch_dir);
  if (! @results) {
    print "no hit found for $query_name\n";
    next;
  }
  chomp $query_name;
  print "Results for $query_name\n";
  my ($check, $refspec_final) = &determineRefspecFinal($query_name, @refspec);
  if ($check == 0) {
    die "error in retrieving refspec data\n";
  }
  if (!defined $hitlimit_local or $hitlimit_local > scalar(@results)) {
	$hitlimit_local = scalar(@results);
}
  for (my $k = 0; $k < $hitlimit_local; $k++) {
    my $hitname = $results[$k];
    print "$hitname\n";
    my $keep = 0;
    my $hitseq = '';
    $refseq = '';
    ## 4b) test for the reciprocity criterion fulfilled
    ($keep, $hitseq)  = &check4reciprocity($query_name, $hitname, $refspec_final, @refspec);
    if ($keep == 1) {
	## blast search with the hmm hit identifies the core-ortholog sequence of the reference species
	## check for the taxon from the taxon_file.txt. IN FACT, THIS ROUTINE IS OUTDATED. I ALWAYS USE THE
	## GLOBAL TAXON NAME FOR THE HAMSTERED SEQUENCES.
      my $taxon = '';
      if ($taxon_check){
	if ($taxon_check == 1) {
	  $taxon = &getTaxon($hitname);
	}
	elsif ($taxon_check == 2) {
	  $taxon = $taxon_global;
	}
      }
      ## put the info about the hits into an object for later post-processing
      ### HERE COMES THE NEW STUFF THAT DEALS WITH THE DIFFERENT POSSIBILITIES: STRICT, RELAXED OR WHATEVER...
      $fileobj = &determineReferences ($fileobj, $taxon, $refspec_final, $hitname, $hitseq, $hitcount);
      $hitcount++;
    }
    else {
      print "Reciprocity not fulfilled!\n";
    }
  }
  ## 5) do the rest only if at least one hit was obtained
  if (defined $fileobj) {
    ## 5a) if the hits are derived from ESTs, get the best ORF
    if ($estflag) {
    $fileobj =  &predictORF();
    }
    ## 5b) if the user has chosen to postprocess the results
    if ($rep) {
      &processHits($fileobj);
    }
    ## 6) prepare the output
    my @taxa = keys(%$fileobj);
    for (my $i = 0; $i< @taxa; $i++) {
      if ($rep) {
	push @newseqs, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
	push @newseqs, $fileobj->{$taxa[$i]}->{refprot};
	if ($estflag) {
	  push @newcds, ">$query_name|$fileobj->{$taxa[$i]}->{refspec_final}|$taxa[$i]|$fileobj->{$taxa[$i]}->{refid}";
	  push @newcds, $fileobj->{$taxa[$i]}->{refcds};
	}
      }
      else {
	my $idobj = $fileobj->{$taxa[$i]}->{ids};
	my $protobj = $fileobj->{$taxa[$i]}->{prot};
	my $cdsobj  = $fileobj->{$taxa[$i]}->{cds};
	my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
	for (my $j = 0; $j < @$idobj; $j++) {
	  push @newseqs, ">$query_name|$refspecobj->[$j]|$taxa[$i]|$idobj->[$j]";
	  push @newseqs, $protobj->[$j];
	  if ($estflag) {
	    push @newcds, ">$query_name|$taxa[$i]|$idobj->[$j]|$refspecobj->[$j]";
	    push @newcds, $cdsobj->[$j];
	  }
	}
      }
      my $refs = $co_seqs->{$query_name};
      for (keys %$refs) {
	my $line = ">$query_name|$_|" . $refs->{$_}->{seqid} . "\n" . $refs->{$_}->{seq};
	push @seqs, $line;
      }
      chomp @seqs;
      print "\n";
      @seqs = (@seqs, @newseqs);
      open (OUT, ">$fa_dir_neu/$query_name.fa");
      print OUT join "\n", @seqs;
      print OUT "\n";
      close OUT;
      if ($estflag) {
	open (OUT, ">$fa_dir_neu/$query_name.cds.fa");
	print OUT join "\n", @newcds;
	close OUT;
      }
      for (my $i = 0; $i < @newseqs; $i+= 2) {
	my $line = $newseqs[$i] . "|" . $newseqs[$i+1];
	$line =~ s/>//;
	push @seqs2store, $line;
	if ($estflag) {
	  my $cdsline = $newcds[$i] . "|" . $newcds[$i+1];
	  $cdsline =~ s/>//;
	  push @cds2store, $cdsline;
	}
      }
    }
  }
}
if (@seqs2store > 0) {
  if ($append) {
    open (OUT, ">>$seqs2store_file") or die "failed to open output file\n";
  }
  else {
    open (OUT, ">$seqs2store_file") or die "failed to open output file\n";
  }
  print OUT join "\n", @seqs2store;
  print OUT "\n";
  close OUT;
  if ($estflag) {
    if ($append) {
      open (OUT, ">>$cds2store_file") or die "failed to open output file\n";
    }
    else {
    open (OUT, ">$cds2store_file") or die "failed to open output file\n";
    }
    print OUT join "\n", @cds2store;
    print OUT "\n";
    close OUT;
  }
}
else {
  print "no hits found\n";
}
exit;
##################### start sub ###############
####### checkInput performs a number of checks whether sufficient information
### and all data are available to run HaMStR
sub checkInput {
	######### check a number of flags that only serve for providing the user with some information
	if (defined $show_coreortholog_sets) {
		## Do nothing but just list all available core ortholog sets in $hmmpath
		my @coresets = (`ls $hmmpath`);
		chomp @coresets;
		print "You can choose from the following core ortholog sets available in $hmmpath:\n\n";
		print join "\n", @coresets;
		print "\n\n";
		exit;
	}
  	my @log;
  	my $check = 1;
  	$dbfile_short = $dbfile;
  	## extract the path from the dbpath if available and prune of trailing '/'
  	if ($dbfile =~ /(.*\/)/) {
		$dbpath = $1;
		$dbpath =~ s/\/$//;
  	}
  	$dbfile =~ s/.*\///;	
  	$dbfile_short =~ s/\..*//;
  	## 
  	## 0) Check for presence of the file with the sequences that should be hamstered
  	if (-e "$dbpath/$dbfile") {
		#the file exists
		push @log, "INFILE PROCESSING\n";
		if (-e "$dbpath/$dbfile.mod") {
			push  @log, "\tA modified infile $dbfile.mod already exists. Using this one\n\n";
		}
		else {
			print "removing newlines from the infile $dbfile such that a sequence forms a consecutive string\n";
			`$path/bin/nentferner.pl -in=$dbpath/$dbfile -out=$dbpath/$dbfile.mod`;
			#`nentferner.pl -in=$dbpath/$dbfile -out=$dbpath/$dbfile.mod`;
		}
		if 	(! -e "$dbpath/$dbfile.mod") { 
			push @log, "Problems running the script nentferner.pl\n";
			$check = 0;
		}
		else {
			print "nentferner.pl succeeded.\n";
			push @log, "\tNewlines from the infile have been removed\n";
			$dbfile = $dbfile . '.mod';
				if (defined $longhead) {
				`$sedprog -i -e "s/[[:space:]]\\+/$idsep/g" -e 's/\\(>.\\{20\\}\\).*/\\1/' $dbfile`;
				push @log, "\tOption -longhead was chosen. Replaced whitespaces in the sequence identifier with '$idsep'\n";
			}
		}
	  }
	  else {
		#the provided infile does not exist:
		push @log, "The specified infile $dbpath/$dbfile does not exist. PLEASE PROVIDE A VALID INFILE!\n";
		$check = 0;
		return ($check, @log);
 	 }
	  ## 1) check for filetype
	  print "Checking for filetype:\t";
	  if (!defined $estflag and !defined $proteinflag) {
	    push @log, "No file sequence type was determined. HaMStR will guess whether EST or protein sequences are analyzed";
		my $seq = `head -n 2 $dbpath/$dbfile |tail -n 1`;
		my $orilength = length($seq);
		$seq =~ s/[AGCTN]//ig;
		if (length($seq) / $orilength >0.1) {
			$proteinflag = 1;
			print "Guessing sequence type: Protein\n";
			push @log, "More than 10% of the first sequence in the file are non-AGCTN. Guessing sequence type: Protein";
		}
		else {
			$estflag = 1;
			print "Guessing sequence type: DNA\n";
		}
    	$check = 1;
	  }
	  if ($estflag) {
      	$estfile = $dbfile;
      	$dbfile = "$dbfile.tc";
      	push @log, "HaMStR will run on the ESTs in $estfile";
      	push @log, "Translating ESTs";
      	if (!(-e "$dbpath/$dbfile")) {
      		print "translating $estfile, this may take a while\n";
  			`$path/bin/translate.pl -in=$dbpath/$estfile -out=$dbpath/$dbfile`;
			open (LOG, "$dbpath/hamstrsearch.log");
			my @info = <LOG>;
			@log = (@log, @info);
			close LOG;
  	    }
  	    else {
      		push @log, "Translated file already exists, using this one\n";
    	}	
    	if (! -e "$dbpath/$dbfile") {
    	 	push @log, "The translation of $estfile failed. Check the script translate.pl";
				print "failed\n";
				$check = 0;
    	}
    	else {
     	 ## file type is protein
     	 print "succeeded\n";
    	}
	  }
	  ## 2) Check for presence of blastall
	  print "Checking for the blast program\t";
	  if (`which $blast_prog` =~ / no /) {
    	push @log, "could not execute $blast_prog. Please check if this program is installed and executable";
    	print "failed\n";
    	$check = 0;
	  }
	  else {
    	push @log, "check for $blast_prog succeeded";
    	print "succeeded\n";
	  }
	  ## 3) Check for presence of hmmsearch
	  print "Checking for hmmsearch\t";
	  my $hmmcheck = `$prog -h |grep -c 'HMMER 3'`;
	  if (! `$prog -h`) {
	  	push @log, "could not execute $prog. Please check if this program is installed and executable";
    	print "failed: $prog is not installed or not executable\n";
    	$check = 0;
	  }
	  elsif ($hmmcheck != 1) {
		push @log, "It seems that $prog is not from the HMMER 3 package. Please check!";
		print "failed: $prog is not from the HMMER 3 package\n";
		$check = 0;
	  }
	  else {
      	push @log, "check for $prog succeeded\n";
      	print "succeeded\n";
	  }
	  ## 3b) Check for genewise
	  print "Checking for genewise\t";
	  if (! `genewise -help`) {
    	push @log, "Could not execute genewise. Please check if this program is installed and executable";
    	print "failed: genewise is not executable\n";
    	$check = 0;
	  }
	  else {
	  	my $gwcheck = `echo \$WISECONFIGDIR`;
	  	if (length($gwcheck) < 1) {
	  		push @log, "The environmental variable WISECONFIGDIR has not been set. I am expecting troubles when invoking genewise. 
Please consult the installation manual for genewise and set this variable";
		print "failed: the environmental variable WISECONFIGDIR has not been set.\n";
		$check = 0;
	  	}
	  	else {
	  		print "succeeded\n";
	  	}
	  }
	  ## 4) Check for presence of the directory structure
	  print "checking for presence of the hmm files:\t";
	  if ( ! defined $hmmset or ! -e "$hmmpath/$hmmset") {
	  	push @log, "You need to specify a valid core ortholog set. Make also sure that you provide the path to this set if it is not in the default location $hmmpath.";
  		print "failed\n";
  		$check = 0;
	  }
	  else {
  		$hmmpath = "$hmmpath/$hmmset";
  		$fafile = "$hmmpath/$hmmset" . '.fa';
  		$hmm_dir = "$hmmpath/$hmm_dir";
  		$hmmsearch_dir = $dbpath .'/hmm_search_' . $dbfile_short . '_' . $hmmset;
  		## 4b) check for the presence of the hmm-files and the fasta-file
  		if (!(-e "$hmm_dir")) {
  			push @log, "Could not find $hmm_dir";
    		print "failed\n";
    		$check = 0;
  		}
  		else {
  			if (defined $hmm) {
    			@hmms = split ',', $hmm;
    			chomp @hmms;
				### check for the presence of all hmms
				for (my $k = 0; $k < @hmms; $k++) {
					if (! -e "$hmm_dir/$hmms[$k]") {
						push @log, "$hmms[$k] has been defined but could not be found in $hmm_dir/$hmms[$k]";
						$check = 0;
						last;
      					}
      				else {
      					push @log, "$hmms[$k] has been found\n";
     	 			}	
				}
			}	
    		else {
      			push @log, "running HaMStR with all hmms in $hmm_dir";
     	 		@hmms = `ls $hmm_dir`;
    		}
    		chomp @hmms;
    		print "succeeded\n";
 	 	}
  	}
  	## 6) Test for presence of the fasta file containing the sequences of the core-ortholog cluster
  	print "checking for presence of the core-ortholog file:\t";
  	if (defined $fafile) {
    	if (! -e "$fafile") {
      		push @log, "Could not find the file $fafile";
      		print "failed\n";
      		$check = 0;
    	}
    	else {
    		push @log, "check for $fafile succeeded";
      		print "succeeded\n";
    	}
  	}
  	else {
    	push @log, "Please provide path and name of fasta file containing the core-ortholog sequences";
    	$check = 0;
    	print "failed\n";
  	}
  	## 7) Checks for the taxon_file
  	print "testing whether the taxon has been determined:\t";  
  	if (!(defined $taxon_file) or (!(-e "$taxon_file"))) {
    	if (defined $taxon_global) {
      		push @log, "using default taxon $taxon_global for all sequences";
      		print "succeeded\n";
      		$taxon_check = 2;
    	}
    	else {
      		push @log, "No taxon_file found. Please provide a global taxon name using the option -taxon";
     		 print "failed\n";
      		$check = 0;
    	}
  	}
  	else {
    	push @log, "using the file $taxon_file as taxon_file";
    	print "succeeded\n";
    	$taxon_check = 1;
  	}
  	## 8) Check for reference taxon
  	print "Checking for reference species and blast-dbs\t";
  	if (!(defined $refspec_string) and (! defined $strict and ! defined $relaxed)) {
      push @log, "Please provide a reference species for the reblast!";
      print "failed\n";
      $check = 0;
  	}
  	elsif (defined $strict or defined $relaxed) {
      if (! defined $refspec_string) {
      	## The user has not provided a string of reference taxa. Chose all from the fasta file containing
	  	## the core orthologs. 
	  	@refspec = `grep '>'  $fafile |cut -d '|' -f 2 |sort |uniq`;
	 	chomp @refspec;
		$refspec_string = join ',', @refspec;
      }
      else {
      	@refspec = split (/,/, $refspec_string);
      }
      if ($strict) {
      	push @log, "Strict flag has been set. Reference species for the reblast: $refspec_string";
      }
      else {
      	push @log, "Relaxed flag has been set. Reference species for the reblast: $refspec_string";
      }
      if (@refspec == 0) {
	  	print "failed\n";
	 	 $check = 0;
      }
      else {
	  print "succeeded\n";
      }
  	}
  	else {
  		push @log, "Reference species for the re-blast: $refspec_string";
    	@refspec = split(/,/, $refspec_string);
    	$refspec_name = $refspec[0];
    	print "succeeded\n";
  	}
  	## 9) Check for presence of the required blast dbs
  	print "checking for blast-dbs:\t";
  	for (my $i = 0; $i < @refspec; $i++) {
    	my $blastpathtmp = "$blastpath/$refspec[$i]/$refspec[$i]".'_prot';
    	if (! (-e "$blastpathtmp.pin")) {
      		push @log, "please edit the blastpath. Could not find $blastpathtmp";
      		print "$blastpathtmp failed\n";
      		$check = 0;
    	}
    	else {
      		push @log, "check for $blastpathtmp succeeded";
      		print "succeeded\n";
    	}
  	}
  	## 10) Set the file where the matched seqs are found
  	my $strictstring = '';
  	if (defined $strict) {
        $strictstring = '.strict';
	}
  	$seqs2store_file = $dbpath . '/hamstrsearch_' . $dbfile_short . '_' . $hmmset . $strictstring . '.out';
  	$cds2store_file = $dbpath . '/hamstrsearch_' . $dbfile_short . '_' . $hmmset . '_cds' . $strictstring . '.out';

  	## 11) check for filter setting for BLAST
  	print "checking for low complexity filter setting:\t";
	$filter =~ tr/ft/FT/;
	if ($filter ne 'T' and $filter ne 'F') {
		push @log, "Filter is set to $filter. Please set the low complexity filter either to F or T.";
		print "low complexity filter check failed\n";
		$check = 0;
   	}
	else {
	push @log, "check for low complexity filter setting succeeded. Chosen value is $filter";
	print "succeeded\n";
	}

  	## 12) apply the evalue-cut-off to the hmmsearch program
  	$prog = $prog . " -E $eval";
  	push @log, "hmmsearch: $prog";

  	## 12b) hit limit for the re-blast
  	if ($hitlimit) {
  		push @log, "re-blast hit_limit: $hitlimit";
  	}
  	else {
		push @log, "re-blast hit_limit: none applied";
  	}
 	## 13) setting up the directories where the output files will be put into.
 	$fa_dir_neu = $dbpath . '/fa_dir_' . $dbfile_short . '_' . $hmmset . '_' . $refspec[0];
  	$tmpdir = $dbpath . '/' . $tmpdir;

  	if ($strict) {
      $fa_dir_neu .= '_strict';
  	}
  	if ($relaxed) {
      $fa_dir_neu .= '_relaxed';
  	}
  	if ($check == 1) {
    	if (!(-e "$hmmsearch_dir")) {
      		`mkdir "$hmmsearch_dir"`;
    	}
    	if (!(-e "$fa_dir_neu")) {
      		`mkdir "$fa_dir_neu"`;
    	}
    	if (!(-e "$tmpdir")) {
      		`mkdir "$tmpdir"`;
    	}
  	}
  	## 14) determin whether or not the -representative flag has been set
  	if (defined $rep) {
		push @log, "HaMStR will run with the -representative option";
	}
  	else {
		push @log, "HaMStR was called without the -representative option. More than one ortholog may be identified per core-ortholog group!";
	} 
  	return ($check, @log);
}
#################
## check4reciprocity is the second major part of the program. It checks
## whether the protein sequence that has been identified by the hmmsearch
## identifies in turn the protein from the reference taxon that was used to
## build the hmm.
sub check4reciprocity {
  my ($query_name, $hitname, $refspec_final, @refspec) = @_;
  my $searchdb;
  my $strict_suc = -1; # keeps track of success for all taxa
  my $relaxed_suc = 0; # keeps track of success for at least one taxon
  ## get the sequence that was identified as hit in the pHMM search from the db_file
  my $hitseq = `grep -m 1 -A 1 ">$hitname\$" $dbfile | tail -n 1`;
  if (!defined $hitseq) {
    print "could not retrieve a sequence for $hitname. Skipping...\n";
    return(0, '', '', '');
  }
  
  ## continue with the blast
  chomp $hitseq;
  ## now run the blast
  open (OUT, ">$tmpdir/$$.fa") or die "could not open out for writing\n";
  print OUT ">$hitname\n$hitseq";
  close OUT;
  
  ## now comes the new part that does one to many blast searches. We need to iterate through all
  ## entries in the file $refspec_final and perform the Blast against each reftaxon. Note, unless
  ## $strict or $relaxed flags are set, there will be only a single reftaxon. If $relaxed is chosen
  ## then we can stop the blast searches as soon as the reciprocity is fulfilled.
  for (my $k = 0; $k < @$refspec_final; $k++) {
    my $orthocount = $refspec_final->[$k]->{orthocount};
    ## 1) Perform the blast search with the k-th reftaxon
    print "Reftaxon: $refspec_final->[$k]->{refspec}\n";
    if ($blast_prog =~ /blastp/) {
      !`$blast_prog -db $refspec_final->[$k]->{searchdb} -num_descriptions 10 -num_alignments 10 -query $tmpdir/$$.fa  -out $tmpdir/$$.blast` or die "Problem running blastp\n";
    ### postprocess the outfile
	`$sedprog -i -e 's/Length=\\([0-9]*\\)/  (\\1 letters)/' -e 's/^\\(>*\\)lcl|/\\1/' $tmpdir/$$.blast`;	
    }
    else {
      !`blastall -p blastp -d $refspec_final->[$k]->{searchdb} -F $filter -i $tmpdir/$$.fa -o $tmpdir/$$.blast` or die "Problem running blast\n";
    }
    ## 2) now parse the best blast hit

    my @hits = &getBestBlasthit("$tmpdir/$$.blast");
    if (@hits > 0) {
      my $idsref = $refspec_final->[$k]->{refid};
      my @original_ids = @$idsref;
      print "core_orthologs: ", join "\t", @original_ids , "\n";
      ## now loop through the best hits with the same evalue and check whether
      ## among these I find the same seq as in $original
      my $i = 0;
      my $suc = 0; # keeps track of success for a single taxon
      while ($suc == 0 and $i <@hits) {
	print "blast-hit: $hits[$i]";
	## now loop through all the refspec-sequences in the hmm file
	my $j = 0;
	while ($suc == 0 and $j < @original_ids) {
	  if ($original_ids[$j] eq $hits[$i]) {
	    print "\thitting\n";
	    $refspec_final->[$k]->{hit} = $j;
	    $suc = 1;
	    $relaxed_suc = 1;
	  }
	  else {
	    print "\nnot hitting $original_ids[$j]\n";
	    $j ++;
	  }
	  if ($suc == 1) {
	    $relaxed_suc = 1;
	    if ($strict_suc == -1) {
	      $strict_suc = 1;
	    }
	  }
	}
	$i++;
      }
      if ($suc == 0) {
	$strict_suc = 0; # none of the blast hits matched against the
	# the reftaxon seq
      }
    }
    else {
      print "no hit obtained\n";
      $strict_suc = 0;
    }	
    ## when the user has chosen the strict flag, there is no reason to continue when $suc
    ## has remained 0 (reciprocity criterion not fulfilled). Thus, return to main.
    if ($strict and $strict_suc == 0) {
      return (0, $hitseq);
    }
  }

  if ($relaxed_suc == 1) {
    return (1, $hitseq);
  }
  else {
    return (0, $hitseq);
  }
}
#############
sub getBestBlasthit {
    my @hits;
    my ($file) = @_;
    my $searchio = Bio::SearchIO->new(-file        => $file,
				      -format      => 'blast',
				      -report_type => 'blastp') or die "parse failed";
    while(my $result = $searchio->next_result){
	my $count = 0;
	my $sig;
	my $sig_old;
	while( my $hit = $result->next_hit){
	    ## now I enter all top hits having the same evalue into the result
	    $sig = $hit->score;
	    if (!defined $sig_old) {
		$sig_old = $sig;
	    }
	    if ($sig == $sig_old) {
		push @hits, $hit->accession;
	    }
	    else {
		last;
	    }
	}
    }
    return(@hits);
}
##################
sub getTaxon {
    my ($hitname) = @_;
#    my $q = "select name from taxon t, est_project e, est_info i, annotation_neu a where a.id = $hitname and a.contig_id = i.contig_id and i.project_id = e.project_id and e.taxon_id = t.taxon_id";
    if ($hitname =~ /\D/) {
	$hitname =~ s/_.*//;
    }
    my $taxon = `grep -m 1 "^$hitname," $taxon_file | $sedprog -e 's/^.*,//'`;
    chomp $taxon;
    $taxon =~ s/^[0-9]+,//;
    $taxon =~ s/\s*$//;
    $taxon =~ s/\s/_/g;
    if ($taxon) {
	return ($taxon);
    }
    else {
	return();
    }
}
###############
sub determineReferences {
    my ($fileobj, $taxon, $refspec_final, $hitname, $hitseq, $hitcount) = @_;
    my $refseq = '';
    my $refspec;
    ## now we have to distinguish between three cases:
    ## 1) hamstr is running in normal mode and one refspec has been determined. In this case, $refspec_final
    ## contains data only from a single species.
    ## 2) hamstr is running in normal mode and alternative refspecs have been determined by the user.
    ## $refspec_final may contain results from more than one species, but we need to consider only the first
    ## entry.
    ## 3) hamstr is running in the strict mode. In this case $refspec_final contains data from several taxa and we need
    ## to select the taxon and sequence that is most similar to the hamstered sequence.
    ## 4) hamstr is running in the relaxed mode. In this case $refspec_final may contain data from several taxa and
    ## we need to select the taxon and the sequence that is most similar to the hamstered sequence.
    if (defined $strict or defined $relaxed) {
	## more than one refspec. Now find the one that fits best
	my $max_score = 0;
	for (my $i = 0; $i < @$refspec_final; $i++) {
	    ## first, check whether the reciprocity criterion has been fulfilled
	    if (defined $refspec_final->[$i]->{hit}) {
		my $rcn = $refspec_final->[$i]->{hit};
		my $refseq_cand = $refspec_final->[$i]->{sequence}->[$rcn];
		my $refspec_cand_id = $refspec_final->[$i]->{refid}->[$rcn];
		my $refspec_cand = $refspec_final->[$i]->{refspec};
		my $score = &getAlignmentScore($refseq_cand, $hitseq);
		if ($score > $max_score) {
		    $refspec = $refspec_cand;
		    $refseq = $refseq_cand;
		    $max_score = $score;
		}
	    }
	}
    }
    else { ## no choice, just one refspec
	my $rcn = $refspec_final->[0]->{hit};
	$refseq = $refspec_final->[0]->{sequence}->[$rcn];
	$refspec = $refspec_final->[0]->{refspec};
    }
    $fileobj->{$taxon}->{prot}->[$hitcount] = $hitseq;
    $fileobj->{$taxon}->{ids}->[$hitcount] = $hitname;
    $fileobj->{$taxon}->{refseq}->[$hitcount]= $refseq;
    $fileobj->{$taxon}->{refspec}->[$hitcount] = $refspec;
    return($fileobj);
}
###############
sub processHits {
  my ($fileobj) = @_; 
  ## 1) align all hit sequences for a taxon against the reference species
  my @taxa = keys(%$fileobj);
  for (my $i = 0; $i < @taxa; $i++) {
    &orfRanking($taxa[$i]);
  }
}  
  

################
sub predictORF {
  my $fileobj_new;
  my @taxa = keys(%$fileobj);
  for (my $i = 0; $i < @taxa; $i++) {
    my $protobj = $fileobj->{$taxa[$i]}->{prot};
    my $idobj = $fileobj->{$taxa[$i]}->{ids};
    my $refseqobj = $fileobj->{$taxa[$i]}->{refseq};
    my $refspecobj = $fileobj->{$taxa[$i]}->{refspec};
    my @ids = @$idobj;
    for (my $j = 0; $j < @ids; $j++) {
	my $refseq = $refseqobj->[$j];
	my $refspec = $refspecobj->[$j];
      ## determine the reading frame
      my ($rf) = $ids[$j] =~ /.*_RF([0-9]+)/;
	print "rf is $rf\n";
      $ids[$j] =~ s/_RF.*//;
      my $est = `grep -A 1 ">$ids[$j]\\b" $estfile |tail -n 1`;
      if (! $est) {
	die "error in retrieval of est sequence for $ids[$j] in subroutine processHits\n";
      }
      ## the EST is translated in rev complement
      if ($rf > 3) {
	$est = revComp($est);
      }
      ### debuggin IUB code
      if ($est =~ /[^AGCT]/i) {
	$est =~ s/[^AGCTagct]/n/g;
      }
	print "running genewise\n";
      my $gw = run_genewise_hamstr->new($est, $refseq, "$tmpdir");
      my $translation = $gw->translation;
      my $cds = $gw->codons;
      $translation =~ s/[-!]//g;
      $fileobj_new->{$taxa[$i]}->{ids}->[$j] = $ids[$j];
      $fileobj_new->{$taxa[$i]}->{prot}->[$j] = $translation;
      $fileobj_new->{$taxa[$i]}->{cds}->[$j] = $cds;
      $fileobj_new->{$taxa[$i]}->{refseq}->[$j] = $refseq;
      $fileobj_new->{$taxa[$i]}->{refspec}->[$j] = $refspec;
    }
  }
  return($fileobj_new);
}
############################
sub orfRanking {
  my ($spec) = @_;
  my $result;
  my $refprot;
  my $refcds;
  my @toalign;
  my $protobj = $fileobj->{$spec}->{prot};
  my $idobj = $fileobj->{$spec}->{ids};
  my $refcluster; ## variables to take the cluster and its id for later analysis
  my $refid;
  if (@$protobj == 1) {
    ## nothing to chose from
    $refprot = $protobj->[0];
    $refcds = $fileobj->{$spec}->{cds}->[0];
    my $length = length($refprot);
    $refid = $idobj->[0] . "-" . $length;
  }
  else {
    ## more than one cluster
      ## note, I set the refseq fix to the first entry. This is to avoid that in this routine 
      ## sequences from different taxa are used.  
      push @toalign, ">$fileobj->{$spec}->{refspec}->[0]";
      push @toalign, $fileobj->{$spec}->{refseq}->[0];
      ## now walk through all the contigs
      for (my $i = 0; $i < @$protobj; $i++) {
	  my @testseq = (">$idobj->[$i]", $protobj->[$i]);
	  @testseq = (@testseq, @toalign);
	  open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
	  print OUT join "\n", @testseq;
	  close OUT;
	  ## run clustalw
	  !(`$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running clustalw\n";
	  ## get the alignment score
	  $result->[$i]->{score} =  `grep "Alignment Score" $tmpdir/$pid.ref.log |$sedprog -e 's/[^0-9]//g'`;
	  if (!$result->[$i]->{score}) {
	      die "error in determining alignment score\n";
	  }
	  chomp $result->[$i]->{score};
	  ## get the aligned sequence
	  open (ALN, "$tmpdir/$pid.ref.aln") or die "failed to open alignment file\n";
	  my @aln = <ALN>;
	  close ALN;
	  my $aseq = extractSeq($idobj->[$i], @aln);
	  ## remove the terminal gaps
	  $aseq =~ s/-*$//;
	  $result->[$i]->{aend} = length $aseq;
	  my ($head) = $aseq =~ /^(-*).*/;
	  ($result->[$i]->{astart}) = length($head)+1;
      }
      ### the results for all seqs has been gathered, now order them
      $result = &sortRef($result);
      ($refprot, $refcds, $refid) = &determineRef($result,$spec);
  }
  $fileobj->{$spec}->{refprot} = $refprot;
  $fileobj->{$spec}->{refcds}  = $refcds;
  $fileobj->{$spec}->{refid}   = $refid;
  $fileobj->{$spec}->{refspec_final} = $fileobj->{$spec}->{refspec}->[0];
  return();
}
###########################
sub sortRef {
    my $result = shift;
    my @sort;
    for (my $i = 0; $i < @$result; $i++) {
	push @sort, "$i,$result->[$i]->{astart},$result->[$i]->{aend},$result->[$i]->{score}";
    }
    open (OUT, ">$tmpdir/$pid.sort") or die "failed to write for sorting\n";
    print OUT join "\n", @sort;
    close OUT;
    `sort -n -t ',' -k 2 $tmpdir/$pid.sort >$tmpdir/$pid.sort.out`;
    @sort = `less $tmpdir/$pid.sort`;
    chomp @sort;
    $result = undef;
    for (my $i = 0; $i < @sort; $i++) {
	($result->[$i]->{id}, $result->[$i]->{start}, $result->[$i]->{end}, $result->[$i]->{score}) = split ',', $sort[$i];
    }
    return($result);
}
########################
sub determineRef {
  my ($result, $spec) = @_;
  my $lastend = 0;
  my $lastscore = 0;
  my $final;
  my $count = 0;
  my $id = '';
  for (my $i = 0; $i < @$result; $i++) {
    if ($result->[$i]->{start} < $lastend or $lastend == 0) {
      if ($result->[$i]->{score} > $lastscore) {
	$lastend = $result->[$i]->{end};
	$lastscore = $result->[$i]->{score};
	$id = $result->[$i]->{id};
      }
    }
    elsif ($result->[$i]->{start} > $lastend) {
      ## a new part of the alignment is covered. Fix the results obtained so far
      $final->[$count]->{id} = $id;
      $lastend = $result->[$i]->{end};
      $id = $result->[$i]->{id};
      $count++;
    }
  }
  $final->[$count]->{id} = $id;
  ## now concatenate the results
  my $refprot = '';
  my $refid = '';
  my $refcds = '';
  for (my $i = 0; $i < @$final; $i++) {
    my $seq = $fileobj->{$spec}->{prot}->[$final->[$i]->{id}];
    my $cdsseq = $fileobj->{$spec}->{cds}->[$final->[$i]->{id}];
    my $length = length($seq);
    $refid .= "$fileobj->{$spec}->{ids}->[$final->[$i]->{id}]-$length" . "PP";
    $refprot .= $seq;
    if ($estflag) {
      $refcds .= $cdsseq;
    }
  }
  $refid =~ s/PP$//;
  return($refprot, $refcds, $refid);
}
#############################
sub extractSeq {
  my ($id, @aln) = @_;
  my $seq = '';
  my $start = 0;
  for (my $i = 0; $i < @aln; $i++) {
    if ($aln[$i] =~ $id) {
      $start = 1;
    }
    elsif ($aln[$i] =~ />/ and $start == 1) {
      last;
    }
    elsif ($start == 1) {
      $seq .= $aln[$i];
    }
  }
  $seq =~ s/\s//g;
  return ($seq);
}
##############################
sub revComp {
    my ($seq) = @_;
    chomp($seq);
    $seq =~ tr/AGCTYRKMWSagct/TCGARYMKWSTCGA/;
    $seq = reverse($seq);
    return($seq);
}
##############################
sub parseHmmer3 {
  my ($file, $path) = @_;
  if (!defined $path) {
    $path = '.';
  }
  open (IN, "$path/$file") or die "failed to open $file\n";
  my @data = <IN>;
  close IN;
  ### extract the hits
  my @hit;
  my $start = 0;
  my $stop = 0;
  my $i = 0;
  for (my $i = 0; $i < @data; $i++) {
    if (!($data[$i] =~ /\S/)) {
      next;
    }
    else {
      if ($data[$i] =~ /Scores for complete sequence/) {
	$start = 1;
	$i += 4;
      }
      elsif (($data[$i] =~ /inclusion threshold/) or ($data[$i] =~ /Domain/i)) {
	last;
      }
      if ($start == 1 and $stop == 0) {
	$data[$i] =~ s/^\s+//;
	my @list = split /\s+/, $data[$i];
	push @hit, $list[8];
	if (@hit == 1 and $bhh) {
		last;
	}
      }
    }
  }
  ### get the query_id
  my ($query) = grep /^Query:/, @data;
  $query =~ s/^Query:\s+//;
  $query =~ s/\s.*//;
  if (defined $hit[0]) {
    chomp @hit;
    return ($query, @hit);
  }
  else {
    return ($query);
  }
}
#####################
sub parseSeqfile {
  my $seqref;
  my $id;
  my $spec;
  my $seqid;
  my $seq;
  my $file = shift;
  open (IN, "$file") or die "failed to open $file\n";
  my @seqs = <IN>;
  close IN;
  chomp @seqs;
  for (my $i = 0; $i < @seqs; $i++) {
    if ($seqs[$i] =~ />/) {
	$seqs[$i] =~ s/>//;
      if (defined $id and defined $seq) {
	$seqref->{$id}->{$spec}->{seqid} = $seqid;
	$seqref->{$id}->{$spec}->{seq} = $seq;
	$seq = undef;
      }
      ($id, $spec, $seqid) = split (/\|/, $seqs[$i]);
    }
    else {
      $seq .= $seqs[$i];
    }
  }
  if (defined  $id and defined $seq) {
	$seqref->{$id}->{$spec}->{seqid} = $seqid;
	$seqref->{$id}->{$spec}->{seq} = $seq;
	$seq = undef;
      }
  return ($seqref);
}
##################
sub getAlignmentScore{ 
    my ($refseq_cand, $hitseq) = @_;
    my @testseq = ('>hitseq', $hitseq, '>refseq', $refseq_cand);
    open (OUT, ">$tmpdir/$pid.ref.fa") or die "could not open file for writing refseqs\n";
    print OUT join "\n", @testseq;
    close OUT;
    ## run clustalw
    !(`$alignmentprog -infile=$tmpdir/$pid.ref.fa -output=fasta -outfile=$tmpdir/$pid.ref.aln 2>&1 >$tmpdir/$pid.ref.log`) or die "error running clustalw\n";
    ## get the alignment score
    my $score =  `grep "Alignment Score" $tmpdir/$pid.ref.log |$sedprog -e 's/[^0-9]//g'`;
    if (!$score) {
	die "error in determining alignment score! Problem with ClustalW\n";
    }
    chomp $score;
    return ($score);
}
######################3
sub determineRefspecFinal {
  my ($query_name, @refspec) = @_;
  my $refspec_final;
  ## now get the id and the sequence used for building the hmm. Note, the latter will be
  ## needed at a later step to determine the best hit
  my @original;
  my $ac = 0;
  for (my $i = 0; $i < @refspec; $i++) {
    @original = `grep -A 1 "^>$query_name|$refspec[$i]" $fafile |$sedprog -e "s/.*$refspec[$i]\|//"`;
    chomp @original;
    
    if (@original > 0) {
      $refspec_final->[$ac]->{refspec} = $refspec[$i];
      $refspec_final->[$ac]->{searchdb} = "$blastpath/$refspec[$i]/$refspec[$i]" . "_prot";
      ## now allow for more than one sequence per core-ortholog cluster and species
      $refspec_final->[$ac]->{orthocount} = 0;
      for (my $j = 0; $j < @original; $j+= 2) {
	$refspec_final->[$ac]->{refid}->[$refspec_final->[$ac]->{orthocount}] = $original[$j];
	$refspec_final->[$ac]->{sequence}->[$refspec_final->[$ac]->{orthocount}] = $original[$j+1];
	$refspec_final->[$ac]->{orthocount} += 1;
      }
      $ac++;
      @original = qw();
      if (!defined $strict and !defined $relaxed) {
	## one reftaxon is enough
	last;
      }
    }
    else {
      print "original sequence not be found with grepping for ^>$query_name|$refspec[$i]. Proceeding with next refspec\n";
    }
  }
  if (! defined $refspec_final->[0]->{refid}) {
    print "original sequence not found\n";
    return (0, $refspec_final);
  } 
  ## now print some wordy information...
  if (!defined $strict and !defined $relaxed) {
    print "REFSPEC is $refspec_final->[0]->{refspec}\n";
  }
  return(1, $refspec_final);
}

############## co-ortholog prediction using a alignment score criterion as in InParanoid
sub identifyCoorthologs{
	
}
