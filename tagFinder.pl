#!/usr/bin/perl -w

use strict;
use warnings;

sub usage()
{
print STDERR << "EOF";

usage: perl $0 [-lsSiNODIXLEcewWrRh] -f fastqFile
       [-t tagsFile] [-h headpieces] [-o overhangs] [-p closingprimers]
       [-q basequality] [-a anchorsize]
       [-v valid patterns] [-V invalid patterns]
       [-d 1.ddd-www_2.ddd-www_3.ddd-www]
       [-T ReadsToTest] [-x threads]

what it does:
- counts all tag combinations found in sequence reads
- detects over-represented tag combinations' tendencies
- removes duplicate counts if degenerated bases are found

what it needs:
- [mandatory] raw sequencing results in fastq standard format
- [mandatory] tag inventary in a 2 column-tabulated text file
  format: CycleId(dot)NumericId(tab)sequence (e.g. 1.001 ACGTACGT)
- [optional] configuration file defining experiments' information
  format: FastqFile ReverseCycles(0/1) TagFilePath HeadPieces Overhangs ClosingPrimers ValidTags InvalidTags
  note1: multiple tagfiles, headpieces, overhangs and/or closing primers must be splitted by comma (",")
  note2: valid and invalid tags accepts comma (",") separated lists and also regex patterns

what to know:
- options -slr are recommended for low quality data (e.g. nanopore)
- max errors allowed are 1 indel/tagstring and 1bp variant/tagcycle
- valid and invalid tag info is recommended for similar reads options

by Jorge Amigo Lechuga

OPTIONS:
  -f  fastq file (plain text or gzipped)
  -t  tags file
  -h  headpieces, splitted by ","
  -o  overhangs, splitted by ","
  -p  closing primers, splitted by ","
  -a  anchor size (default: 7)
  -q  minimum base quality (default: 0)
      10 means 1 in 10   (90%)
      20 means 1 in 100  (99%)
      30 means 1 in 1000 (99.9%)
  -l  allow left-anchored-only reads
  -s  detect similar tags (1 error/tagcycle allowed)
  -S  detect very similar tags (1 error/tagstring allowed)
  -i  reverse paired cycles
  -N  DO NOT deal with degen region at all
  -O  DO NOT detect overrepresented tendencies
  -D  DO NOT clean degen region for possible sequencing errors (limit: 10000)
  -I  print invalid reads
  -X  print existing tags
  -L  print tag lengths
  -E  print error types (for similar search)
  -c  print chimeric counts
  -e  print expected library tags
  -w  print raw counts
  -v  valid tags
  -V  invalid tags
  -W  leave invalid tags out from analysis
  -d  degen count distribution for a particular tag combination
  -r  recovery mode (try to find new matches in processed reads)
  -R  print recovered cycles
  -T  number of reads to test analysis (default: 0)
  -x  number of threads (default: 1)
  -h  this (help) message

EOF
exit;
}

my $start = time;
my $elapsed = $start;
my $command = my $commandOptions = my $commandString = "";
$commandOptions .= "$_ " foreach (@ARGV);
if ($commandOptions) {
	$command = "$0 $commandOptions";
	($commandString = $command) =~ s/ -/ \\\n-/g;
}

# get options
use Getopt::Std;
my %opt = ();
my $opt_string = 'f:t:h:o:p:a:q:lsSiNODIXLEcewv:V:Wd:rRT:x:h';
getopts("$opt_string", \%opt) or usage();
usage() if $opt{h};
my $fastqFile = "";
if ( $opt{f} ) { $fastqFile = $opt{f} } else { usage() }
my $tagFiles = $opt{t} ? $opt{t} : "";
my $headpieces = $opt{h} ? $opt{h} : "";
my $overhangs = $opt{o} ? $opt{o} : "";
my $closingprimers = $opt{p} ? $opt{p} : "";
my $mbq = defined $opt{q} ? $opt{q} : 0;
my $anchorSize = $opt{a} ? $opt{a} : 7;
my $leftAnchored = $opt{l} ? 1 : 0;
my $similarTags = $opt{s} ? 1 : 0;
my $similarTagsStrict = $opt{S} ? 1 : 0; if ($similarTagsStrict) { $similarTags = 1 }
my $reverseCycles = $opt{i} ? 1 : 0;
my $detectDegen = $opt{N} ? 0 : 1;
my $printOver = $opt{O} ? 0 : 1;
my $cleanDegen = $opt{D} ? 0 : 1; $cleanDegen = 0 unless $detectDegen;
my $printInvalid = $opt{I} ? 1 : 0;
my $printExistingTags = $opt{X} ? 1 : 0;
my $printLengths = $opt{L} ? 1 : 0;
my $printErrorTypes = $opt{E} ? 1 : 0;
my $printChimeras = $opt{c} ? 1 : 0;
my $printExpected = $opt{e} ? 1 : 0;
my $printRawCounts = $opt{w} ? 1 : 0;
my $validTags = $opt{v} ? $opt{v} : "";
my $invalidTags = $opt{V} ? $opt{V} : "";
my $invalidTagsOut = $opt{W} ? 1 : 0;
my $tags2focus = $opt{d} ? $opt{d} : "";
my $recoveryMode = $opt{r} ? 1 : 0;
my $printRecovery = $opt{R} ? 1 : 0;
my $testReads = $opt{T} ? $opt{T} : 0;
my $threads = $opt{x} ? $opt{x} : 1;

# print degen regions counts
my $printDegen = 0;

# max number of compounds to sort at output
my $outSortLimit = 100000;

# non recursive tag library search if multiple tag files present
# 1 - saves a little time not reanalyzing matched reads (recommended)
# 0 - allows detecting different libraries' joined reads (very rare)
#     (note that short and invalid reads count will not be accurate)
my $nonRecursiveLibSearch = 1;

# initialize the label for empty anchor tags
my $noFwAnchorTag = "";

# get prefix from file name
my $prefix = "";
if ($fastqFile =~ /([^\/]+)\.[fastq]+(\.gz)?(\.(\d+))?/i) {
	$prefix = $1;
	$prefix .= ".$4" if $4;
	$prefix =~ s/\.unaligned\.\d+//;
	$prefix =~ s/\.aligned\.\d+//;
	$prefix =~ s/\.unaligned//;
	$prefix =~ s/\.aligned//;
	$prefix =~ s/.+user_(XEN|MX)-+\d+-//;
	print STDERR "$prefix\n";
} else {
	die "could not find a valid prefix in $fastqFile";
}

# deal with threading
my $fastqFiles = 0;
foreach (<$fastqFile.*>) { $fastqFiles++ }
my $threaded = 0;
if ($threads > 1) {
	if ($fastqFile =~ /\.\d+$/) {
		$threaded = 1;
	} else {
		unless ($fastqFiles) {
			split_and_fork($fastqFile, $prefix, $threads);
			foreach (<$fastqFile.*>) { $fastqFiles++ }
		}
	}
}

# fastq qualities
my $lowQualRegex = "";
if ($mbq) {
	# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
	my $fastqQuals = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~';
	my $lowQualChars = substr($fastqQuals, 0, $mbq + 1);
	$lowQualRegex = qr/[$lowQualChars]/
}

# prepare strings and file headers for overrepresented sequences
my $tagsOverHeader = "";
my $tagsOverTabBlanks = "";
my @tagsOverTypes = ();
my @overStrings = ();
if ($printOver) {
	$tagsOverHeader = "\tSDCOUNT_RAW\tSDCOUNT_DEDUP";
	$tagsOverTabBlanks = "\t0\t0";
	@tagsOverTypes = ("raw", "dedup", "unique");
	@overStrings = ("LINES", "PLANES");
	foreach my $overString (@overStrings) {
		foreach my $tagsOverType (@tagsOverTypes) {
			$tagsOverHeader .= "\tOVER_".uc($tagsOverType)."_$overString";
			$tagsOverTabBlanks .= "\t0";
		}
	}
} else {
	$tagsOverHeader = "";
	$tagsOverTabBlanks = "";
}

# find relative paths
my $configFile = "";
($configFile = $0) =~ s/([^\/]+)\.pl$/$1.config.ini/;
my $tagsPath = "";
($tagsPath = $0) =~ s/bins\/[^\/]+$/tags\//;

# read config file
if (open CONFIG, $configFile) {
	<CONFIG>;
	while (<CONFIG>) {
		if (/^[^#;]\S/) {
			s/[\r\n]//g;
			my ($f, $i, $t, $h, $o, $p, $v, $V) = split "\t";
			if ($fastqFile =~ /$f/i or $f =~ /$fastqFile/i) {
				$reverseCycles = $reverseCycles ? $reverseCycles : $i;
				$tagFiles = $tagFiles ? $tagFiles : $t;
				$headpieces = $headpieces ? $headpieces : $h if $h;
				$overhangs = $overhangs ? $overhangs : $o if $o;
				$closingprimers = $closingprimers ? $closingprimers : $p if $p;
				$validTags = $validTags ? $validTags : $v if $v;
				$invalidTags = $invalidTags ? $invalidTags : $V if $V;
				last;
			}
		}
	}
	close CONFIG;
}
usage() unless $tagFiles;

# convert dots to pattern
$validTags =~ s/\./\\./g;
$invalidTags =~ s/\./\\./g;
my $checkValidTags = $validTags or $invalidTags ? 1 : 0;

# check for everything needed
my @tagFilesArray = split ",", $tagFiles;
my @headpiecesArray = split ",", $headpieces;
my @overhangsArray = split ",", $overhangs;
my @closingprimersArray = split ",", $closingprimers;
my @validTagsArray = split ",", $validTags;
my @invalidTagsArray = split ",", $invalidTags;

# loop through tag files
my $tagFileCounter = 0;
my %readsProcessed = ();
my $totalReads = 0;
my $validReads = 0;
my $maxLengthSurround = 0;
my $maxTagLength = 0;
my $maxLengthCount = 0;
my $openedReads = 0;
my $forwardReads = 0;
my $reverseReads = 0;
my $similarReads = 0;
my $shorterReads = 0;
my $longerReads = 0;
my $reducedReads = 0;
my $lowQualReads = 0;
my $invalidReads = 0;
my $recoverReads = 0;
my %similarCount = ();
my $chimeraReads = 0;
my $unfoundReads = 0;
my $undedupReads = 0;
my $matchedReads = 0;
my $dedupedReads = 0;
my $matchedReadsRecovered = 0;
my $logFile = "tags_$prefix.log";
my %expectedMatchUniq = ();
my %expectedMatchCounts = ();
my %expectedDedupMatchCounts = ();
my %nonexpectedMatchUniq = ();
my %nonexpectedMatchCounts = ();
my %nonexpectedDedupMatchCounts = ();
my %librarySizeCP = ();
my %libraryTagsCPstats = ();
foreach my $tagFile (@tagFilesArray) {

# count tag files
$tagFileCounter++;

# reset primers' arrays
@headpiecesArray = split ",", $headpieces;
@closingprimersArray = split ",", $closingprimers;

# look for expected libraries inside the tag file name (all other libraries will be considered not expected)
my %selectedLibs = ();
my $tagFileLabel = "";
my $prefixLabel = "";
my @libs2select = split ";", $tagFile;
$tagFile = shift @libs2select;
($tagFile, $tagFileLabel) = split ":", $tagFile;
foreach (@libs2select) { $selectedLibs{$_} = 1 }

if ($tagFileLabel) {
	$prefixLabel = "$prefix.$tagFileLabel";
} else {
	$prefixLabel = $prefix;
}

# define file names
my $txtFileDWallTags = "tags_$prefixLabel.allTags.dwar";
my $txtFileDWfilter = "tags_$prefixLabel.filtered.dwar";
my $invalidFile = "tags_$prefixLabel.invalid.txt";
my $chimFile = "tags_$prefixLabel.chimeras.txt";
my $distFile = "tags_$prefixLabel.$tags2focus.txt";
my $overTagsFile = "tags_$prefixLabel.over.txt";
my $existingTagsFile = "tags_$prefixLabel.existingtags.txt";
my $lengthsFile = "tags_$prefixLabel.lengths.txt";
my $degenFile = "tags_$prefixLabel.degen.txt";
my $errorTypeFile = "tags_$prefixLabel.errors.txt";
my $recoveryFile = "tags_$prefixLabel.recovery.txt";
my $expectedFile = "tags_$prefixLabel.expected.txt";
my $tagCountsFile = "tags_$prefixLabel.tagcounts.txt";

# read tag files
my %headerIndex = ();
my %libraryHPs = ();
my %libraryCPs = ();
my %libraryCPId = ();
my %libraryTagsCP = ();
my %validTagCodesCP = ();
my %invalidTagCodesCP = ();
my %tags = ();
my %cycleTags = ();
my %cycleLengths = ();
my $existingTags = "";
my %libraryTags = ();
my %validTagCodes = ();
if ( ! -f $tagFile and -f "$tagsPath$tagFile" ) { $tagFile = "$tagsPath$tagFile" }
open TAGS, $tagFile or die "could not open $tagFile: $!";
while (<TAGS>) {
	if (/^#ID\tSEQUENCE\t/) {
		s/[\r\n]//g;
		my @cols = split "\t";
		$checkValidTags = 1;
		for (my $i = 2; $i <= $#cols; $i++) {
			$headerIndex{$i} = $cols[$i];
		}
	} elsif (/^HPL/) {
		s/[\r\n]//g;
		my @cols = split "\t";
		my $headPiece = $cols[1];
		for (my $i = 2; $i <= $#cols; $i++) {
			if ( $cols[$i] ) {
				$libraryHPs{$headPiece} = 1;
			}
		}
	} elsif (/^CPL/) {
		s/[\r\n]//g;
		my @cols = split "\t";
		my $closingPrimer = "$cols[0]-$cols[1]";
		(my $closingPrimerId = $closingPrimer) =~ s/N[^-]+$//;
		for (my $i = 2; $i <= $#cols; $i++) {
			if ( $cols[$i] ) {
				$libraryCPs{$closingPrimer} = 1;
				$libraryCPId{ $headerIndex{$i} }{$closingPrimerId} = 1;
			}
		}
	} elsif (/(^\S*?([0-9]+)[.-][0-9]+)\t([ACGT]+)/) {
		my $tagCode = $1;
		my $cycle = $2;
		my $tag = $3;
		if ($reverseCycles && $cycle % 2 == 0) {
			($tag = reverse $tag) =~ tr/ACGT/TGCA/;
		}
		if ( exists $tags{$cycle}{$tag} ) {
			$existingTags .= "$tag already exists: $tagCode vs $tags{$cycle}{$tag}\n";
		} else {
			$cycleTags{$cycle} = ";" unless exists $cycleTags{$cycle};
			$cycleTags{$cycle} .= "$tag;";
			$tags{$cycle}{$tag} = $tagCode;
			if (%libraryCPId) {
				s/[\r\n]//g;
				my @cols = split "\t";
				for (my $i = 2; $i <= $#cols; $i++) {
					if ( $cols[$i] ) {
						if ( not %selectedLibs or exists $selectedLibs{ $headerIndex{$i} } ) {
							foreach my $closingPrimerId ( keys %{ $libraryCPId{ $headerIndex{$i} } } ) {
								$libraryTagsCP{$closingPrimerId}{$cycle}++;
								$validTagCodesCP{$closingPrimerId}{$tagCode} = 1;
							}
						} else {
							foreach my $closingPrimerId ( keys %{ $libraryCPId{ $headerIndex{$i} } } ) {
								$libraryTagsCP{$closingPrimerId}{$cycle}++;
								$invalidTagCodesCP{$closingPrimerId}{$tagCode} = 1;
							}
						}
					}
				}
			}
			if ($checkValidTags) {
				my $validFlag = 0;
				foreach (@validTagsArray) {
					my @cols = split ";";
					my $validTag = $cols[$#cols];
					if ($tagCode =~ /$validTag/) {
						if ($#cols) {
							for (my $i = 0; $i < $#cols; $i++) {
								my $closingPrimerId = $cols[$i];
								$closingPrimerId =~ s/N[^-]+$//;
								unless ( exists $validTagCodesCP{$closingPrimerId}{$tagCode} ) {
									$libraryTagsCP{$closingPrimerId}{$cycle}++;
									$validTagCodesCP{$closingPrimerId}{$tagCode} = 1;
								}
							}
						} else {
							$validFlag = 1;
						}
					}
				}
				foreach (@invalidTagsArray) {
					my @cols = split ";";
					my $invalidTag = $cols[$#cols];
					if ($tagCode =~ /$invalidTag/) {
						if ($#cols) {
							for (my $i = 0; $i < $#cols; $i++) {
								my $closingPrimerId = $cols[$i];
								$closingPrimerId =~ s/N[^-]+$//;
								if ( exists $validTagCodesCP{$closingPrimerId}{$tagCode} ) {
									$libraryTagsCP{$closingPrimerId}{$cycle}--;
									delete $validTagCodesCP{$closingPrimerId}{$tagCode};
									$invalidTagCodesCP{$closingPrimerId}{$tagCode} = 1;
								} else {
									foreach my $validCpId (keys %validTagCodesCP) {
										if ($validCpId =~ /$closingPrimerId/) {
											if ( exists $validTagCodesCP{$validCpId}{$tagCode} ) {
												$libraryTagsCP{$validCpId}{$cycle}--;
												delete $validTagCodesCP{$validCpId}{$tagCode};
												$invalidTagCodesCP{$validCpId}{$tagCode} = 1;
											}
										}
									}
								}
							}
						} else {
							$validFlag = 0;
						}
					}
				}
				if ($validFlag) {
					$libraryTags{$cycle}++;
					$validTagCodes{$tagCode} = 1;
				}
			}
			unless (%libraryCPId or $checkValidTags) {
				$libraryTags{$cycle}++;
			}
			if ( not exists $cycleLengths{$cycle} ) {
				$cycleLengths{$cycle} = length($tag);
			} elsif (length($tag) != $cycleLengths{$cycle}) {
				die "diferent tag lengths in cycle $cycle";
			}
		}
	} else {
		die "unknown tag format $_";
	}
}
close TAGS;

# review primers
if (%libraryHPs) {
	foreach my $headPiece (@headpiecesArray) {
		die "headpiece $headPiece not found in $tagFile" if not exists $libraryHPs{$headPiece};
	}
	@headpiecesArray = sort keys %libraryHPs;
}
if (%libraryCPs) {
	foreach my $closingPrimer (@closingprimersArray) {
		die "closing primer $closingPrimer not found in $tagFile" if not exists $libraryCPs{$closingPrimer};
	}
	@closingprimersArray = sort keys %libraryCPs;
}

# calculate the expected size for a tag combination
my $cycles = 0;
my $cycleLength = 0;
my $tagsLength = 0;
my @sortedCycles = ();
my @sortedCyclesLengths = ();
foreach my $cycle (sort {$a<=>$b} keys %cycleLengths) {
	$cycles++;
	push @sortedCycles, $cycle;
	$cycleLength = $cycleLengths{$cycle};
	push @sortedCyclesLengths, $cycleLength;
	$tagsLength += $cycleLength;
}
my $overhangLength = 0;
my @overhangsLengths = ();
foreach my $overhang (@overhangsArray) {
	$overhangLength = length($overhang);
	push @overhangsLengths, $overhangLength;
	$tagsLength += $overhangLength;
}
my $tagsLengthTop = $tagsLength + 1;
my $tagsLengthLow = $tagsLength - 1;
my $anchoredTagsLength = $tagsLength + $anchorSize;
if ($similarTags) { $anchoredTagsLength++ }

# deal with head pieces
unless ( @headpiecesArray ) { die "headpieces needed" }
my $primerNumber = 0;
my $headPieceNumber = 0;
my $minPrimerLength = 0;
my @anchor5regexes = ();
my @anchor3regexes = ();
my @staticSeqs = ();
my @staticSeqLengths = ();
my @anchor5regexesSimilar = ();
my @anchor3regexesSimilar = ();
my %regexNcp = ();
my $maxDegenErrors = 0;
my @closingPrimerIds = ();
foreach my $headPiece (@headpiecesArray) {
	$primerNumber++;
	$headPieceNumber++;
	if ($anchorSize > length($headPiece)) {
		die "anchor size must be smaller than " . length($headPiece) . ", cause head piece";
	}
	# store minimum primer size
	$minPrimerLength = length($headPiece) unless $minPrimerLength;
	$minPrimerLength = length($headPiece) if length($headPiece) < $minPrimerLength;
	# extract anchors from headPieces
	(my $headPieceRev = reverse $headPiece) =~ tr/ACGT/TGCA/;
	my $fwAnchor5 = substr($headPiece, -$anchorSize);
	my $rvAnchor3 = substr($headPieceRev, 0, $anchorSize);
	# store anchor patterns
	push @anchor5regexes, qr/$fwAnchor5(.+)/;
	push @anchor3regexes, qr/(.+?)$rvAnchor3/;
	# store static sequences
	push @staticSeqs, $headPiece;
	push @staticSeqLengths, length($headPiece);
	# deal with similar patterns
	if ($similarTags) {
		my $fwAnchor5pattern = "";
		for (my $i = length($fwAnchor5) - 1; $i >= 0; $i--) {
			my $anchorPattern = $fwAnchor5;
			substr($anchorPattern, $i, 1) = substr($anchorPattern, $i, 1)."?.";
			$fwAnchor5pattern .= $anchorPattern . "|";
		}
		chop $fwAnchor5pattern;
		push @anchor5regexesSimilar, qr/($fwAnchor5pattern)(.+)/;
		push @anchor3regexesSimilar, qr/(.{$tagsLengthLow,$tagsLengthTop})$rvAnchor3/;
	}
}

# deal with closing primers
unless ( @closingprimersArray ) { die "closingprimers needed" }
foreach my $closingPrimer (@closingprimersArray) {
	$primerNumber++;
	# extract id
	my $closingPrimerId = "";
	if ($closingPrimer =~ /^(\S*-)(\S+)/) {
		$closingPrimerId = $1;
		$closingPrimer = $2;
	} else {
		$closingPrimerId = "";
	}
	# check for degenerated bases on closing primer
	my $staticSeq = "";
	if ($closingPrimer =~ /^([ACGT]+)(N+)(.+)/) {
		my $anchorN = $1;
		$staticSeq = $3;
		if ($anchorSize > length($anchorN)) {
			die "anchor size must be smaller than " . length($anchorN) . ", cause non degen region";
		}
		my $lengthN = length($2);
		$closingPrimerId .= $anchorN;
		$regexNcp{$closingPrimerId} = qr/$anchorN(.{$lengthN})/;
		if ($maxDegenErrors) {
			if ( $lengthN > $maxDegenErrors + 1 ) { $maxDegenErrors = $lengthN - 1 }
		} else {
			$maxDegenErrors = $lengthN - 1;
		}
	} else {
		$staticSeq = $closingPrimer;
		if ($anchorSize > length($closingPrimer)) {
			die "anchor size must be smaller than " . length($closingPrimer) . ", cause closing primer";
		}
		$closingPrimerId .= $noFwAnchorTag;
		$regexNcp{$closingPrimerId} = "";
	}
	# store closing primer id
	push @closingPrimerIds, $closingPrimerId;
	# store minimum primer size
	$minPrimerLength = length($closingPrimer) if length($closingPrimer) < $minPrimerLength;
	# extract anchors from closingPrimers
	(my $closingPrimerRev = reverse $closingPrimer) =~ tr/ACGT/TGCA/;
	my $fwAnchor3 = substr($closingPrimer, 0, $anchorSize);
	my $rvAnchor5 = substr($closingPrimerRev, -$anchorSize);
	# store anchor patterns
	push @anchor3regexes, qr/(.+?)$fwAnchor3/;
	push @anchor5regexes, qr/$rvAnchor5(.+)/;
	# store static sequences
	($staticSeq = reverse $staticSeq) =~ tr/ACGT/TGCA/;
	push @staticSeqs, $staticSeq;
	push @staticSeqLengths, length($staticSeq);
	# deal with similar patterns
	if ($similarTags) {
		my $rvAnchor5pattern = "";
		for (my $i = length($rvAnchor5) - 1; $i >= 0; $i--) {
			my $anchorPattern = $rvAnchor5;
			substr($anchorPattern, $i, 1) = substr($anchorPattern, $i, 1)."?.";
			$rvAnchor5pattern .= $anchorPattern . "|";
		}
		chop $rvAnchor5pattern;
		push @anchor3regexesSimilar, qr/(.{$tagsLengthLow,$tagsLengthTop})$fwAnchor3/;
		push @anchor5regexesSimilar, qr/($rvAnchor5pattern)(.+)/;
	}
}
push @closingPrimerIds, $noFwAnchorTag;
my $minSeqLength = $minPrimerLength + $anchoredTagsLength;

# make sure tag info depends on closing primer
foreach my $closingPrimerId (@closingPrimerIds) {
	foreach my $cycle (keys %libraryTags) { $libraryTagsCP{$closingPrimerId}{$cycle} += $libraryTags{$cycle} }
	foreach my $tagCode (keys %validTagCodes) { $validTagCodesCP{$closingPrimerId}{$tagCode} = 1 }
}

# get library size
foreach my $closingPrimerId (keys %libraryTagsCP) {
	$librarySizeCP{$closingPrimerId} = 1;
	foreach my $cycle ( keys %{ $libraryTagsCP{$closingPrimerId} } ) {
		$libraryTagsCPstats{$closingPrimerId}{$cycle} = $libraryTagsCP{$closingPrimerId}{$cycle};
		$librarySizeCP{$closingPrimerId} *= $libraryTagsCP{$closingPrimerId}{$cycle};
	}
}

# check for valid overhangs
if (@overhangsArray) {
	if ($cycles > $#overhangsArray + 2 || $cycles < $#overhangsArray + 1) {
		die "incorrect number of overhangs";
	}
	push @overhangsArray, "";
}

my %cpMatchesCounts = ();
my %matchedCPreads = ();
my %cpMatchesStrand = ();
my %cpMatchesDedupStrings = ();
my %degenCounts = ();
$totalReads = 0;
# $validReads = 0;
# $maxLengthSurround = 0;
$maxTagLength = 0;
$maxLengthCount = 0;
# $openedReads = 0;
# $forwardReads = 0;
# $reverseReads = 0;
# $similarReads = 0;
$shorterReads = 0;
# $longerReads = 0;
# $reducedReads = 0;
# $lowQualReads = 0;
$invalidReads = 0;
# $recoverReads = 0;
# %similarCount = ();
# $chimeraReads = 0;
# $unfoundReads = 0;
# $undedupReads = 0;
# $matchedReads = 0;
# $dedupedReads = 0;
# $matchedReadsRecovered = 0;
my %beginningSeqs = ();
my %recoverReadCount = ();
my %tagLengths = ();
my %matchedTags = ();
my %unmatchedTags = ();
my %chimeras = ();
my %similarTypes = ();
my %beginningSeqsBaseErrors = ();
my %cpMatchesDedups = ();
my %countsUniq = ();
my %countsAverage = ();
my %tagsFound = ();
my %structuresFound = ();
my %countsStdDev = ();
my %tagsOver = ();
my %cpMatchesSD = ();
my %idTagsMissing = ();
my $filterMatches = 0;
if ($fastqFiles) {
	foreach my $splitResult (<tags_$prefix.*.allTags.dwar>) {
		open SPLIT, "gunzip -fc $splitResult |" or die "could not read $splitResult: $!";
		# open SPLIT, $splitResult or die "could not read $splitResult: $!";
		while (<SPLIT>) {
			my ($match, $closingPrimerId, $matchCount, $strandCount, $degenSeqs) = split ",";
			my $cpMatch = "$closingPrimerId;$match";
			$cpMatchesCounts{$cpMatch} += $matchCount;
			$matchedCPreads{$closingPrimerId} += $matchCount;
			$cpMatchesStrand{$cpMatch} += $strandCount;
			if ($degenSeqs) {
				$cpMatchesDedupStrings{$cpMatch} .= "$degenSeqs;";
				if ($printDegen) {
					foreach my $degenSeq (split ";", $degenSeqs) {
						$degenCounts{$closingPrimerId}{$degenSeq}++;
					}
				}
			}
		}
		close SPLIT;
		unlink $splitResult;
	}
	foreach my $splitResult (<tags_$prefix.*.log>) {
		open SPLIT, $splitResult or die "could not read $splitResult: $!";
		while (<SPLIT>) {
			if (/^Total:\s+(\d+)/) { $totalReads += $1 }
			elsif (/^Valid:\s+(\d+)/) { $validReads += $1 }
			elsif (/^Almost:\s+(\d+)/) { $maxLengthSurround += $1 }
			elsif (/^MaxTagLength:\s+(\d+)/) { $maxTagLength += $1 / $fastqFiles; }
			elsif (/^MaxLengthCount:\s+(\d+)/) { $maxLengthCount += $1 }
			elsif (/^Opened:\s+(\d+)/) { $openedReads += $1 }
			elsif (/^Forward:\s+(\d+)/) { $forwardReads += $1 }
			elsif (/^Reverse:\s+(\d+)/) { $reverseReads += $1 }
			elsif (/^Similar:\s+(\d+)/) { $similarReads += $1 }
			elsif (/^Shorter:\s+(\d+)/) { $shorterReads += $1 }
			elsif (/^Longer:\s+(\d+)/) { $longerReads += $1 }
			elsif (/^Reduced:\s+(\d+)/) { $reducedReads += $1 }
			elsif (/^LowQual:\s+(\d+)/) { $lowQualReads += $1 }
			elsif (/^Invalid:\s+(\d+)/) { $invalidReads += $1 }
			elsif (/^Recover:\s+(\d+)/) { $recoverReads += $1 }
			elsif (/^Deletions:\s+(\d+)/) { $similarCount{"del"} += $1 }
			elsif (/^Insertions:\s+(\d+)/) { $similarCount{"ins"} += $1 }
			elsif (/^1bVariants:\s+(\d+)/) { $similarCount{"var"} += $1 }
			elsif (/^Chimera:\s+(\d+)/) { $chimeraReads += $1 }
			elsif (/^Unfound:\s+(\d+)/) { $unfoundReads += $1 }
			elsif (/^Undedup:\s+(\d+)/) { $undedupReads += $1 }
			elsif (/^Matched:\s+(\d+)/) { $matchedReads += $1 }
			elsif (/^Matched recovered:\s+(\d+)/) { $matchedReadsRecovered += $1 }
		}
		close SPLIT;
		unlink $splitResult;
	}
	foreach (<$fastqFile.*>) { unlink $_; }
} else {
	# get reads
	if ($printInvalid) { open INVALID, ">$invalidFile" or die "could not write $invalidFile: $!" }
	if (-f $fastqFile) {
		open FASTQ, "seqtk seq $fastqFile |"
		or open FASTQ, "gunzip -fc $fastqFile |"
		or die "could not open $fastqFile";
	} else {
		die "could not find $fastqFile";
	}
	while (<FASTQ>) {
		my $seq = <FASTQ>;
		<FASTQ>;
		my $qual = <FASTQ>;
		$totalReads++;
		if ($tagFileCounter > 1) { if ($nonRecursiveLibSearch) { next if exists $readsProcessed{$totalReads} } }
		my $recoverRead = 0;
		if ($testReads) {
			if ($totalReads > $testReads) {
				print STDERR "\n* TEST RUN: $testReads reads analyzed *\n\n";
				$totalReads--;
				last;
			}
		}
	RECOVERY:
		my $forward = "";
		my $tagString = "";
		my $staticSeq = "";
		my $shortRead = 0;
		my $openedRead = 0;
		my $tagPosMatch = 0;
		my $closingPrimerId = "";
		my $jMin = 0;
		my $jMax = 0;
		my $staticSeqLength = 0;
		if (length($seq) < $minSeqLength) {
			$shorterReads++ unless $recoverRead;
			next;
		} else {
			ANCHORLOOP: for (my $i = 0; $i < $primerNumber; $i++) {
				# qr/$fwAnchor5(.+)/ qr/$rvAnchor5(.+)/
				if ($seq =~ $anchor5regexes[$i]) {
					$openedRead = 1;
					my $anchoredSeq = $1;
					if (length($anchoredSeq) < $anchoredTagsLength) {
						$shortRead = 1 unless $tagPosMatch;
					} else {
						if ($i < $headPieceNumber) {
							$forward = 1; $jMin = $headPieceNumber; $jMax = $primerNumber;
						} else {
							$forward = 0; $jMax = $headPieceNumber;
						}
						$shortRead = 0;
						$tagPosMatch = $-[0] + $anchorSize;
						for (my $j = $jMin; $j < $jMax; $j++) {
							# qr/(.+?)$rvAnchor3/ qr/(.+?)$fwAnchor3/
							if ($anchoredSeq =~ $anchor3regexes[$j]) {
								$openedRead = 0;
								$tagString = $1;
								$staticSeq = $staticSeqs[$i];
								$staticSeqLength = $staticSeqLengths[$i];
								$closingPrimerId = $forward ? $closingPrimerIds[$j - $headPieceNumber] : $closingPrimerIds[$i - $headPieceNumber];
								last ANCHORLOOP;
							} elsif ($leftAnchored) {
								$tagString = $anchoredSeq;
							}
						}
					}
				}
			}
			if ($similarTags and not $tagString) {
				SIMILARLOOP: for (my $i = 0; $i < $primerNumber; $i++) {
					# qr/($fwAnchor5pattern)(.+)/ qr/($rvAnchor5pattern)(.+)/
					if ($seq =~ $anchor5regexesSimilar[$i]) {
						$openedRead = 1;
						my $anchoredSeq = $2;
						if (length($anchoredSeq) < $anchoredTagsLength) {
							$shortRead = 1 unless $tagPosMatch;
						} else {
							if ($i < $headPieceNumber) {
								$forward = 1; $jMin = $headPieceNumber; $jMax = $primerNumber;
							} else {
								$forward = 0; $jMin = 0; $jMax = $headPieceNumber;
							}
							$shortRead = 0;
							$tagPosMatch = $-[0] + $anchorSize;
							for (my $j = $jMin; $j < $jMax; $j++) {
								# qr/(.{$tagsLengthLow,$tagsLengthTop})$rvAnchor3/ qr/(.{$tagsLengthLow,$tagsLengthTop})$fwAnchor3/
								if ($anchoredSeq =~ $anchor3regexesSimilar[$j]) {
									$openedRead = 0;
									$tagString = $1;
									$staticSeq = $staticSeqs[$i];
									$staticSeqLength = $staticSeqLengths[$i];
									$closingPrimerId = $forward ? $closingPrimerIds[$j - $headPieceNumber] : $closingPrimerIds[$i - $headPieceNumber];
									last SIMILARLOOP;
								} elsif ($leftAnchored) {
									$tagString = $anchoredSeq;
								}
							}
						}
					}
				}
			}
		}
		$openedReads++ if $openedRead and not $recoverRead;
		if ($shortRead) {
			$shorterReads++ unless $recoverRead;
			next;
		} elsif ($tagString) {
			if ($nonRecursiveLibSearch) { $readsProcessed{$totalReads} = 1 }
			unless ($forward) {
				($tagString = reverse $tagString) =~ tr/ACGT/TGCA/;
			}
			if ($cleanDegen) {
				$beginningSeqs{$staticSeq}{ substr($seq, 0, $staticSeqLength) }++ if $staticSeq;
			}
			if ($recoverRead) {
				$recoverReads++;
				$recoverReadCount{$recoverRead}++;
				$recoverReadCount{$recoverRead-1}-- if $recoverRead > 1;
			}
			if ($mbq) {
				if (substr($qual, $tagPosMatch, $tagsLength) =~ $lowQualRegex) {
					$lowQualReads++ unless $recoverRead;
					next;
				}
			}
			my $longer = 0;
			my $shorter = 0;
			my $similar = "";
			my $searchChimeras = 0;
			my @tagStrings = ();
			my @tagStringsErrorPos = ();
			my $tagLength = length($tagString);
			$tagLengths{$tagLength}++;
			if ($tagLength < $tagsLength) {
				if ($similarTags and $tagLength == $tagsLength - 1) {
					for (my $baseno = $tagLength - 1; $baseno >= 0; $baseno--) {
						foreach ("A", "C", "G", "T") {
							my $tempTagString = $tagString;
							substr($tempTagString, $baseno, 1) .= $_;
							push @tagStrings, $tempTagString;
							push @tagStringsErrorPos, $baseno + 2;
						}
					}
					foreach ("A", "C", "G", "T") {
						push @tagStrings, $_.$tagString;
						push @tagStringsErrorPos, 1;
					}
					$similar = "del";
				} else {
					$shorter = 1;
				}
			} elsif ($tagLength > $tagsLength) {
				if ($similarTags and $tagLength == $tagsLength + 1) {
					for (my $baseno = $tagLength - 1; $baseno >= 0; $baseno--) {
						my $tempTagString = $tagString;
						substr($tempTagString, $baseno, 1) = "";
						push @tagStrings, $tempTagString;
						push @tagStringsErrorPos, $baseno + 1;
					}
					$similar = "ins";
				} else {
					$longer = 1;
					$searchChimeras = 1;
					push @tagStrings, $tagString;
				}
			} else {
				push @tagStrings, $tagString;
			}
			my $match = "";
			my $matched = 0;
			my $chimera = 0;
			my $tagStart = 0;
			my $tagStringIndex = 0;
			TAGSTRINGS: foreach my $tagString (@tagStrings) {
				for (my $cycleIndex = 0; $cycleIndex < $cycles; $cycleIndex++) {
					my $cycle = $sortedCycles[$cycleIndex];
					my $cycleLength = $sortedCyclesLengths[$cycleIndex];
					my $tag = substr($tagString, $tagStart, $cycleLength);
					my $overhang = "";
					my $overhangLength = 0;
					my $postTag = "";
					if ( $overhang = $overhangsArray[$cycleIndex] ) {
						$overhangLength = $overhangsLengths[$cycleIndex];
						$postTag = substr($tagString, $tagStart + $cycleLength, $overhangLength);
					} else {
						$overhang = "";
					}
					if ($searchChimeras) {
						if (my $tagCode = $tags{$cycle}{$tag}) {
							my $tagRepeats = $tagString =~ s/$tag$overhang//g;
							if ($tagRepeats > 1) {
								$chimera = 1;
								$match .= "$tagCode-$tag*\t";
							} else {
								$matched++;
								$match .= "$tagCode-$tag\t";
							}
						} else {
							$match .= ".\t";
							$tagStart += $cycleLength + $overhangLength;
						}
					} else {
						unless ($similarTags) {
							next TAGSTRINGS if $overhang ne $postTag;
						}
						if (my $tagCode = $tags{$cycle}{$tag}) {
							if ($checkValidTags and $similar) { next TAGSTRINGS unless exists $validTagCodesCP{$closingPrimerId}{$tagCode} or exists $invalidTagCodesCP{$closingPrimerId}{$tagCode} }
							$matched++;
							$match .= "$tagCode-$tag\t";
							$tagStart += $cycleLength + $overhangLength;
							$matchedTags{$cycle}{$tag}++;
						} elsif ($similarTags) {
							my $similarTagUnmatch = 1;
							unless ($similarTagsStrict and $similar) {
								for (my $tagPos = $cycleLength - 1; $tagPos >= 0; $tagPos--) {
									my $tagPattern = $tag;
									substr($tagPattern, $tagPos, 1) = "[ACGT]";
									if ($cycleTags{$cycle} =~ /;($tagPattern);/) {
										$tagCode = $tags{$cycle}{$1};
										if ($checkValidTags) { next unless exists $validTagCodesCP{$closingPrimerId}{$tagCode} or exists $invalidTagCodesCP{$closingPrimerId}{$tagCode} }
										my $errorPos = $tagStart + $tagPos + 1;
										$similar .= ";var,$errorPos";
										$matched++;
										$match .= "$tagCode-$1\t";
										$tagStart += $cycleLength + $overhangLength;
										$matchedTags{$cycle}{$1}++;
										$similarTagUnmatch = 0;
										last;
									}
								}
							}
							if ($similarTagUnmatch) {
								$unmatchedTags{$cycle}{$tag}++;
								next TAGSTRINGS;
							}
						} else {
							$unmatchedTags{$cycle}{$tag}++;
							next TAGSTRINGS;
						}
					}
				}
				last TAGSTRINGS if $matched == $cycles or $chimera;
				$tagStringIndex++;
			}
			if ($chimera) {
				$chimeraReads++ unless $recoverRead;
				$chimeras{"$closingPrimerId;$match"}++;
			} else {
				if ($shorter) { $reducedReads++ unless $recoverRead; next }
				if ($longer) { $longerReads++ unless $recoverRead; next unless $leftAnchored }
				if ($similar) { $similarReads++ unless $recoverRead }
				if ($forward) { $forwardReads++ unless $recoverRead } else { $reverseReads++ unless $recoverRead }
				$validReads++ unless $recoverRead;
				if ($matched == $cycles) {
					if ($recoverRead) {
						$matchedReadsRecovered++;
					} else {
						$matchedReads++;
					}
					if ($similar) {
						foreach my $errorTypePos (split ";", $similar) {
							if ($errorTypePos) {
								my ($errorType, $errorPos) = split ",", $errorTypePos;
								$errorPos = $tagStringsErrorPos[$tagStringIndex] unless $errorPos;
								$similarTypes{$errorType}{$errorPos}++;
								$similarCount{$errorType}++;
							}
						}
					}
					$matchedCPreads{$closingPrimerId}++;
					my $cpMatch = "";
					if ( my $regexN = $regexNcp{$closingPrimerId} ) {
						$cpMatch = "$closingPrimerId;$match";
						if ($detectDegen) {
							my $seq4degen = "";
							if ($forward) {
								$seq4degen = substr($seq, $tagPosMatch + $tagsLengthLow);
							} else {
								($seq4degen = reverse substr($seq, 0, $tagPosMatch)) =~ tr/ACGT/TGCA/;
							}
							if ($seq4degen =~ $regexN) {
								$cpMatchesDedupStrings{$cpMatch} .= "$1;";
								$degenCounts{$closingPrimerId}{$1}++ if $printDegen;
							} else {
								$undedupReads++ unless $recoverRead;
							}
						} else {
							$undedupReads++ unless $recoverRead;
						}
					} else {
						$cpMatch = "$noFwAnchorTag;$match";
					}
					$cpMatchesCounts{$cpMatch}++;
					if ($forward) {
						$cpMatchesStrand{$cpMatch}++;
					} else {
						$cpMatchesStrand{$cpMatch}--;
					}
				} else {
					$unfoundReads++ unless $recoverRead;
				}
			}
		} else {
			unless ($recoverRead) {
				$invalidReads++;
				if ($printInvalid) { print INVALID "invalid\t$seq" }
			}
		}
		if ($recoveryMode and $tagPosMatch) {
			$seq = substr($seq, $tagPosMatch + $anchoredTagsLength);
			$qual = substr($qual, $tagPosMatch + $anchoredTagsLength);
			$recoverRead++;
			goto RECOVERY;
		}
	}
	close FASTQ;
	if ($printInvalid) { close INVALID }
}

$elapsed = time - $start;
print STDERR "Elapsed time post-patterns: $elapsed seconds\n";

# get the number of compounds
my $cpMatchesCount = keys %cpMatchesCounts;

# deal with tag lengths
my $tagLengthString = "";
if (%tagLengths) {
	foreach my $tagLength (sort {$a<=>$b} keys %tagLengths) {
		my $tagLengthCount = $tagLengths{$tagLength};
		$tagLengthString .= "$tagLength\t$tagLengthCount\n";
		if ($tagLengthCount > $maxLengthCount) {
			$maxTagLength = $tagLength;
			$maxLengthCount = $tagLengthCount;
		}
	}
	if ( exists $tagLengths{$maxTagLength - 1} ) { $maxLengthSurround += $tagLengths{$maxTagLength - 1} }
	if ( exists $tagLengths{$maxTagLength + 1} ) { $maxLengthSurround += $tagLengths{$maxTagLength + 1} }
}

unless ($threaded) {
	
	# deal with tag lengths
	if ($printLengths and $tagLengthString) {
		open LENGTHS, ">$lengthsFile" or die "could not write $lengthsFile: $!";
		print LENGTHS "#LENGTHS\tCOUNTS\n$tagLengthString";
		close LENGTHS;
		print STDERR "$lengthsFile\n";
	}
	
	# print degen counts
	if ($printDegen) {
		open DEGEN, "| sort -k1,1 -k3nr >$degenFile" or die "could not write $degenFile: $!";
		print DEGEN "#CP\tDEGENSEQ\tCOUNT\n";
		foreach my $closingPrimerId (keys %degenCounts) {
			foreach my $degenSeq ( keys %{ $degenCounts{$closingPrimerId} } ) {
				print DEGEN "$closingPrimerId\t$degenSeq\t$degenCounts{$closingPrimerId}{$degenSeq}\n";
			}
		}
		close DEGEN;
		print STDERR "$degenFile\n";
	}

	# print error types
	if ($printErrorTypes) {
		open ERRORTYPE, ">$errorTypeFile" or die "could not write $errorTypeFile: $!";
		foreach my $type (sort keys %similarTypes) {
			foreach my $pos ( sort {$a<=>$b} keys %{ $similarTypes{$type} } ) {
				print ERRORTYPE "$type\t$pos\t$similarTypes{$type}{$pos}\n";
			}
		}
		close ERRORTYPE;
		print STDERR "$errorTypeFile\n";
	}

	# print existing tags
	if ($existingTags and $printExistingTags) {
		open EXISTINGTAGS, ">$existingTagsFile" or die "could not write $existingTagsFile: $!";
		print EXISTINGTAGS $existingTags;
		close EXISTINGTAGS;
		print EXISTINGTAGS "$existingTagsFile\n";
	}

	# print recovery metrics
	if ($printRecovery) {
		open RECOVERY, ">$recoveryFile" or die "could not write $recoveryFile: $!";
		print RECOVERY "RECOVERY\tCOUNTS\n";
		foreach my $recoveryTime (sort {$a<=>$b} keys %recoverReadCount) {
			print RECOVERY "$recoveryTime\t$recoverReadCount{$recoveryTime}\n" if $recoverReadCount{$recoveryTime};
		}
		close RECOVERY;
		print STDERR "$recoveryFile\n";
	}

	# print chimeras
	if ($chimeraReads and $printChimeras) {
		open CHIMERAS, "| sort -k" . ($cycles + 2) . "rn >$chimFile" or die "could not write $chimFile: $!";
		foreach my $cpMatch (keys %chimeras) {
			my ($closingPrimerId, $match) = split ";", $cpMatch;
			print CHIMERAS "$match$closingPrimerId\t$chimeras{$cpMatch}\n";
		}
		close CHIMERAS;
		print STDERR "$chimFile\n";
	}

	if ($matchedReads) {

		# print raw tag counts
		if ($printRawCounts) {
			open TAGCOUNTS, "| sort -k3nr >$tagCountsFile" or die "could not open $tagCountsFile: $!";
			foreach my $cycle ( sort {$a<=>$b} keys %matchedTags ) {
				foreach my $tag ( sort keys %{ $matchedTags{$cycle} } ) {
					print TAGCOUNTS "$cycle\t$tag\t$matchedTags{$cycle}{$tag}\t$tags{$cycle}{$tag}\t$prefix\n";
				}
			}
			foreach my $cycle ( sort {$a<=>$b} keys %unmatchedTags ) {
				foreach my $tag ( sort keys %{ $unmatchedTags{$cycle} } ) {
					print TAGCOUNTS "$cycle\t$tag\t$unmatchedTags{$cycle}{$tag}\t\t$prefix\n";
				}
			}
			close TAGCOUNTS;
			print STDERR "$tagCountsFile\n";
		}

		# print particular tag counts
		if ($tags2focus and $detectDegen) {
			(my $tags2focusString = $tags2focus) =~ s/_/\t/g;
			$tags2focusString .= "\t";
			open DIST, "| sort -k3nr >$distFile" or die "could not write $distFile: $!";
			foreach my $closingPrimerId (sort keys %matchedCPreads) {
				my $cpString = "$closingPrimerId;$tags2focusString";
				if ( my $degenSeqs = $cpMatchesDedupStrings{$cpString} ) {
					chop $degenSeqs;
					my %degenSeqsHash = ();
					foreach my $degenSeq (split ";", $degenSeqs) { $degenSeqsHash{$degenSeq}++ }
					foreach my $degenSeq (keys %degenSeqsHash) {
						print DIST "$closingPrimerId\t$degenSeq\t$degenSeqsHash{$degenSeq}\n";
					}
				}
			}
			close DIST;
			print STDERR "$distFile\n";
		}

		# check error rate looking at static sequences
		if ($cleanDegen) {
			my %beginningSeqsCount = ();
			my %beginningSeqsLD = ();
			foreach my $staticSeq (keys %beginningSeqs) {
				my $staticSeqLength = length $staticSeq;
				foreach my $beginningSeq ( keys %{ $beginningSeqs{$staticSeq} } ) {
					my $ld = minSeqLD($staticSeq, $staticSeqLength, $beginningSeq, length $beginningSeq);
					my $ldCount = $beginningSeqs{$staticSeq}{$beginningSeq};
					$beginningSeqsCount{$staticSeq} += $ldCount;
					$beginningSeqsLD{$staticSeq}{$ld} += $ldCount;
				}
				foreach my $ld ( keys %{ $beginningSeqsLD{$staticSeq} } ) {
					my $beginningSeqsBaseError = $beginningSeqsLD{$staticSeq}{$ld} / ( $beginningSeqsCount{$staticSeq} * $staticSeqLength );
					if ( exists $beginningSeqsBaseErrors{$ld} ) {
						if ( $beginningSeqsBaseError > $beginningSeqsBaseErrors{$ld} ) {
							$beginningSeqsBaseErrors{$ld} = $beginningSeqsBaseError;
						}
					} else {
						$beginningSeqsBaseErrors{$ld} = $beginningSeqsBaseError;
					}
				}
			}
		}

		# look for overrepresented tags and clean degen
		foreach my $cpMatch (keys %cpMatchesCounts) {
			my $dedupCount = 1;
			if ($cleanDegen) {
				if ( my $degenSeqs = $cpMatchesDedupStrings{$cpMatch} ) {
					chop $degenSeqs;
					my %degenSeqsHash = ();
					foreach my $degenSeq (split ";", $degenSeqs) { $degenSeqsHash{$degenSeq}++ }
					$dedupCount = keys %degenSeqsHash;
					if ($dedupCount < 10000) {
						my @degenSeqsArray = sort { $degenSeqsHash{$b} <=> $degenSeqsHash{$a} or $a cmp $b } keys %degenSeqsHash;
						my @degenSeqsRev = reverse @degenSeqsArray;
						pop @degenSeqsArray;
						pop @degenSeqsRev;
						my %degenSeqsRemoved = ();
						foreach my $degenSeq1 (@degenSeqsArray) {
							unless ( exists $degenSeqsRemoved{$degenSeq1} ) {
								my $degenSeqLength1 = length $degenSeq1;
								if ($degenSeqLength1 > $maxDegenErrors) {
									foreach my $degenError (1..$maxDegenErrors) {
										if ( exists $beginningSeqsBaseErrors{$degenError} ) {
											my $degenSeqErrorProb = $degenSeqsHash{$degenSeq1} * $degenSeqLength1 * $beginningSeqsBaseErrors{$degenError};
											foreach my $degenSeq2 (@degenSeqsRev) {
												unless ( exists $degenSeqsRemoved{$degenSeq2} ) {
													if ($degenSeqsHash{$degenSeq2} < $degenSeqErrorProb) {
														if (minSeqLD($degenSeq1, $degenSeqLength1, $degenSeq2, length $degenSeq2, $degenError) <= $degenError) {
															$degenSeqsRemoved{$degenSeq2} = 1;
															$dedupCount--;
														}
													} else { last }
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			$cpMatchesDedups{$cpMatch} = $dedupCount;
			my $matchCount = $cpMatchesCounts{$cpMatch};
			my ($closingPrimerId, $match) = split ";", $cpMatch;
			if ($printOver) {
				$countsUniq{$closingPrimerId}++;
				$countsAverage{$closingPrimerId}{raw} += $matchCount;
				$countsAverage{$closingPrimerId}{dedup} += $dedupCount;
			}
			my @matchedTags = split "\t", $match;
			for (my $cycleIndex1 = 0; $cycleIndex1 <= $#matchedTags; $cycleIndex1++) {
				(my $tag1 = $matchedTags[$cycleIndex1]) =~ s/.+-//;
				$tagsFound{$closingPrimerId}{"$cycleIndex1;$tag1"} = 1;
				if ($printOver) {
					$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{raw} += $matchCount;
					$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{dedup} += $dedupCount;
					$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{unique}++;
					for (my $cycleIndex2 = $cycleIndex1 + 1; $cycleIndex2 <= $#matchedTags; $cycleIndex2++) {
						(my $tag2 = $matchedTags[$cycleIndex2]) =~ s/.+-//;
						$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{raw} += $matchCount;
						$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{dedup} += $dedupCount;
						$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{unique}++;
					}
				}
			}
		}
		undef %cpMatchesDedupStrings;
		if ($printOver) {
			my %tempAverage = ();
			my %tempStdDev = ();
			my %tempCount = ();
			foreach my $closingPrimerId (keys %tagsFound) {
				foreach my $type ("raw", "dedup") {
					$countsAverage{$closingPrimerId}{$type} /= $countsUniq{$closingPrimerId};
					$tempAverage{$closingPrimerId}{$type} = $countsAverage{$closingPrimerId}{$type};
				}
			}
			foreach my $cpMatch (keys %cpMatchesCounts) {
				my ($closingPrimerId, $match) = split ";", $cpMatch;
				my $matchCount = $cpMatchesCounts{$cpMatch};
				my $dedupCount = $cpMatchesDedups{$cpMatch};
				foreach my $type ("raw", "dedup") {
					my $tempValue = 0;
					if ($type eq "raw") { $tempValue = $matchCount } else { $tempValue = $dedupCount }
					$tempStdDev{$closingPrimerId}{$type} += ( $tempAverage{$closingPrimerId}{$type} - $tempValue ) ** 2;
					$tempCount{$closingPrimerId}{$type}++;
				}
			}
			foreach my $closingPrimerId (keys %tagsFound) {
				foreach my $type ("raw", "dedup") {
					$countsStdDev{$closingPrimerId}{$type} = sqrt( $tempStdDev{$closingPrimerId}{$type} / $tempCount{$closingPrimerId}{$type} );
				}
			}
			open OVER, ">$overTagsFile" or die "could not open $overTagsFile: $!";
			print OVER "#TYPE\tCOUNT\tSDCOUNT\tCP\tTAGS\n";
			foreach my $type (@tagsOverTypes) {
				foreach my $closingPrimerId (keys %structuresFound) {
					foreach my $structure ("planes", "lines") {
						my $cpCycleTagInfo = $structuresFound{$closingPrimerId}{$structure};
						my @cpCycleTags = keys %$cpCycleTagInfo;
						my $cpCycleTagStdDev = 0;
						my $cpCycleTagCount = 0;
						my $cpCycleTagAverage = 0;
						foreach my $cpCycleTag (@cpCycleTags) {
							$cpCycleTagCount++;
							$cpCycleTagAverage += $cpCycleTagInfo->{$cpCycleTag}{$type};
						}
						$cpCycleTagAverage /= $cpCycleTagCount;
						foreach my $cpCycleTag (@cpCycleTags) {
							$cpCycleTagStdDev += ( $cpCycleTagAverage - $cpCycleTagInfo->{$cpCycleTag}{$type} ) ** 2;
						}
						$cpCycleTagStdDev = sqrt($cpCycleTagStdDev / $cpCycleTagCount);
						my $cpCycleTagThreshold = $cpCycleTagAverage + $cpCycleTagStdDev;
						foreach my $cpCycleTag (@cpCycleTags) {
							my $tempValue = $cpCycleTagInfo->{$cpCycleTag}{$type};
							if ($tempValue > $cpCycleTagThreshold) {
								my ($cycleIndex1, $tag1, $cycleIndex2, $tag2) = split ";", $cpCycleTag;
								my $tag1Code = $tags{ $sortedCycles[$cycleIndex1] }{$tag1};
								my $tagCodeTag = "";
								if ($tag2) {
									my $tag2Code = $tags{ $sortedCycles[$cycleIndex2] }{$tag2};
									$tagCodeTag = "$tag1Code-$tag1;$tag2Code-$tag2";
								} else {
									$tagCodeTag = "$tag1Code-$tag1";
								}
								my $tempCount = 1;
								my $tempStdDev = $cpCycleTagThreshold + $cpCycleTagStdDev;
								while ($tempValue > $tempStdDev) {
									$tempCount++;
									$tempStdDev += $cpCycleTagStdDev;
								}
								$tagsOver{$closingPrimerId}{$type}{$tagCodeTag} = $tempCount;
								print OVER "$type\t$tempValue\t$tempCount\t$closingPrimerId\t$tagCodeTag\n"
							}
						}
					}
				}
			}
			foreach my $cpMatch (keys %cpMatchesCounts) {
				my ($closingPrimerId, $match) = split ";", $cpMatch;
				my $matchCount = $cpMatchesCounts{$cpMatch};
				my $dedupCount = $cpMatchesDedups{$cpMatch};
				my $sdString = "";
				foreach my $type ("raw", "dedup") {
					my $tempValue = 0;
					my $tempCount = 0;
					if ($type eq "raw") {
						$tempValue = $matchCount;
					} else {
						$tempValue = $dedupCount;
					}
					while ($tempValue > $countsAverage{$closingPrimerId}{$type} + ($tempCount + 1) * $countsStdDev{$closingPrimerId}{$type} ) { $tempCount++ }
					$sdString .= "\t$tempCount";
					if ($tempCount) {
						chop $match;
						$match =~ s/\t/;/g;
						print OVER "$type\t$tempValue\t$tempCount\t$closingPrimerId\t$match\n";
					}
				}
				$cpMatchesSD{$cpMatch} = $sdString;
			}
			close OVER;
			if ($cpMatchesCount < $outSortLimit) {
				system "head -n1 $overTagsFile > $overTagsFile.sortTemp";
				system "tail -n +2 $overTagsFile | sort -k1,1 -k3,3nr -k4,4V -k3,3nr >> $overTagsFile.sortTemp";
				rename "$overTagsFile.sortTemp", $overTagsFile;
			}
			print STDERR "$overTagsFile\n";
		}

		# print expected tags
		if ($checkValidTags) {
			if ($printExpected) {
				open EXPECTED, ">$expectedFile" or die "could not open $expectedFile: $!";
			}
			my $tagsFoundCount = 0;
			my $tagsMissingCount = 0;
			for (my $cycleIndex = 0; $cycleIndex <= $#sortedCycles; $cycleIndex++) {
				my $cycle = $sortedCycles[$cycleIndex];
				foreach my $tag ( keys %{ $tags{$cycle} } ) {
					my $tagCode = $tags{$cycle}{$tag};
					print EXPECTED "$tagCode\t$tag" if $printExpected;
					foreach my $closingPrimerId (sort keys %matchedCPreads) {
						my $tagFlag = $validTagCodesCP{$closingPrimerId}{$tagCode} ? 1 : 0;
						if ( exists $tagsFound{$closingPrimerId}{"$cycleIndex;$tag"} ) {
							$tagFlag = "+" . $tagFlag;
							$tagsFoundCount++;
						} else {
							$idTagsMissing{$closingPrimerId}{$cycleIndex}{"$tagCode-$tag"} = 1;
							$tagFlag = "-" . $tagFlag;
							$tagsMissingCount++;
						}
						print EXPECTED "\t$tagFlag" if $printExpected;
					}
					print EXPECTED "\n" if $printExpected;
				}
			}
			if ($tagsMissingCount > $tagsFoundCount) { $filterMatches = 1 }
			if ($printExpected) {
				close EXPECTED;
				open EXPECTED, ">$expectedFile.temp" or die "could not open $expectedFile.temp: $!";
				print EXPECTED "ID\tTAG";
				foreach my $closingPrimerId (sort keys %matchedCPreads) {
					my $closingPrimerHeader = $closingPrimerId ? $closingPrimerId : "EMPTY";
					print EXPECTED "\t$closingPrimerHeader";
				}
				print EXPECTED "\n";
				close EXPECTED;
				if ($cpMatchesCount < $outSortLimit) {
					system "sort -k1n $expectedFile >> $expectedFile.temp";
					rename "$expectedFile.temp", $expectedFile;
				} else {
					system "cat $expectedFile.temp $expectedFile > $expectedFile.temp2";
					rename "$expectedFile.temp2", $expectedFile;
					unlink "$expectedFile.temp";
				}
				print STDERR "$expectedFile\n";
			} else {
				unlink $expectedFile;
			}
		}
	}
	$elapsed = time - $start;
	print STDERR "Elapsed time pre-printing: $elapsed seconds\n";
}

# deal with matches
if ($matchedReads) {
	if ($threaded) {
		# open MATCHES, "| gzip >$txtFileDWallTags" or die "could not create $txtFileDWallTags";
		open MATCHES, ">$txtFileDWallTags" or die "could not write $txtFileDWallTags: $!";
		foreach my $cpMatch (keys %cpMatchesCounts) {
			my ($closingPrimerId, $match) = split ";", $cpMatch;
			my $degenSeqs = "";
			if ( $degenSeqs = $cpMatchesDedupStrings{$cpMatch} ) { chop $degenSeqs } else { $degenSeqs = "" }
			print MATCHES "$match,$closingPrimerId,$cpMatchesCounts{$cpMatch},$cpMatchesStrand{$cpMatch},$degenSeqs,\n";
		}
		close MATCHES;
	} else {
		open MATCHES, ">$txtFileDWallTags" or die "could not write $txtFileDWallTags: $!";
		my $directWrite = 0;
		if ($cpMatchesCount > $outSortLimit) { $directWrite = 1 }
		my ($myheader, $dwheader, $dwfooter) = printHeaderFooter($prefix, $cycles, $directWrite);
		if ($filterMatches) { open MATCHES2, ">$txtFileDWfilter" or die "could not write $txtFileDWfilter: $!" }
		if ($directWrite) {
			print MATCHES $dwheader; print MATCHES $myheader;
			if ($filterMatches) {
				print MATCHES2 $dwheader; print MATCHES2 $myheader;
			}
		}
		my %librarySizeNorm = ();
		foreach my $closingPrimerId (@closingPrimerIds) {
			if ( $librarySizeCP{$closingPrimerId} and $matchedCPreads{$closingPrimerId} ) {
				$librarySizeNorm{$closingPrimerId} = $librarySizeCP{$closingPrimerId} / $matchedCPreads{$closingPrimerId};
			} else {
				$librarySizeNorm{$closingPrimerId} = 0;
			}
		}
		foreach my $cpMatch (keys %cpMatchesCounts) {
			my ($closingPrimerId, $match) = split ";", $cpMatch;
			my $matchCount = $cpMatchesCounts{$cpMatch};
			my $dedupCount = 0;
			unless ( $dedupCount = $cpMatchesDedups{$cpMatch} ) { $dedupCount = 1 }
			$dedupedReads += $dedupCount;
			my $strandCount = sprintf( "%.3f", abs( $cpMatchesStrand{$cpMatch} / $matchCount ) );
			my $cpLibrarySizeNorm = $librarySizeNorm{$closingPrimerId};
			my $expected = 0;
			if ($checkValidTags) {
				$expected = 1;
				while ($match =~ /(\S+)-/g) {
					unless ( exists $validTagCodesCP{$closingPrimerId}{$1} ) {
						$expected = 0;
						last;
					}
				}
			}
			if ($expected == 1) {
				$expectedMatchUniq{$closingPrimerId}++;
				$expectedMatchCounts{$closingPrimerId} += $matchCount;
				$expectedDedupMatchCounts{$closingPrimerId} += $dedupCount;
			} elsif ($expected == 0) {
				next if $invalidTagsOut;
				$nonexpectedMatchUniq{$closingPrimerId}++;
				$nonexpectedMatchCounts{$closingPrimerId} += $matchCount;
				$nonexpectedDedupMatchCounts{$closingPrimerId} += $dedupCount;
			}
			my $overTagString = "";
			if ($printOver) {
				$overTagString = $cpMatchesSD{$cpMatch};
				my @matchedTags = split "\t", $match;
				foreach my $type (@tagsOverTypes) {
					my $overTagValue = 0;
					for (my $cycleIndex1 = 0; $cycleIndex1 < $#matchedTags; $cycleIndex1++) {
						my $matchedTag1 = $matchedTags[$cycleIndex1];
						for (my $cycleIndex2 = $cycleIndex1 + 1; $cycleIndex2 <= $#matchedTags; $cycleIndex2++) {
							my $matchedTag2 = $matchedTags[$cycleIndex2];
							if ( exists $tagsOver{$closingPrimerId}{$type}{"$matchedTag1;$matchedTag2"} ) {
								$overTagValue += 1;
								if ($tagsOver{$closingPrimerId}{$type}{"$matchedTag1;$matchedTag2"} > 1) {
									$overTagValue += 0.1;
								}
							}
						}
					}
					$overTagValue .= ".0" if length $overTagValue == 1;
					$overTagString .= "\t$overTagValue";
				}
				foreach my $type (@tagsOverTypes) {
					my $overTagValue = 0;
					foreach (@matchedTags) {
						if ( exists $tagsOver{$closingPrimerId}{$type}{$_} ) {
							$overTagValue += 1;
							if ($tagsOver{$closingPrimerId}{$type}{$_} > 1) {
								$overTagValue += 0.1;
							}
						}
					}
					$overTagValue .= ".0" if length $overTagValue == 1;
					$overTagString .= "\t$overTagValue";
				}
			}
			my $matchCountNorm = sprintf "%.5f", $matchCount * $cpLibrarySizeNorm;
			my $dedupCountNorm = sprintf "%.5f", $dedupCount * $cpLibrarySizeNorm;
			if ($filterMatches) {
				my $outLine = "$match$closingPrimerId\t$matchCount\t$dedupCount\t$strandCount\t$matchCountNorm\t$dedupCountNorm\t$expected$overTagString\n";
				print MATCHES $outLine; print MATCHES2 $outLine;
			} else {
				print MATCHES "$match$closingPrimerId\t$matchCount\t$dedupCount\t$strandCount\t$matchCountNorm\t$dedupCountNorm\t$expected$overTagString\n";
			}
		}
		if ($filterMatches) {
			if ($directWrite) { print MATCHES2 $dwfooter }
			close MATCHES2;
		}
		foreach my $closingPrimerId (sort keys %idTagsMissing) {
			foreach my $cycleIndex ( sort {$a<=>$b} keys %{ $idTagsMissing{$closingPrimerId} } ) {
				foreach my $tagString ( sort keys %{ $idTagsMissing{$closingPrimerId}{$cycleIndex} } ) {
					my $match = "";
					for (0..$cycles-1) {
						if ($cycleIndex == $_) { $match .= "$tagString\t" } else { $match .= "\t" }
					}
					my $expected = 0;
					if ($tagString =~ /(^\S*?[0-9]+[.-][0-9]+)-[ACGT]+/) {
						$expected = $validTagCodesCP{$closingPrimerId}{$1} ? 1 : 0;
					}
					print MATCHES "$match$closingPrimerId\t0\t0\t0\t0\t0\t$expected$tagsOverTabBlanks\n";
				}
			}
		}
		if ($directWrite) { print MATCHES $dwfooter }
		close MATCHES;
		unless ($directWrite) {
			system "sort -k" . ($cycles + 2) . "rn $txtFileDWallTags > $txtFileDWallTags.sortTemp";
			system "cat $prefix.dwheader $prefix.myheader $txtFileDWallTags.sortTemp $prefix.dwfooter > $txtFileDWallTags";
			unlink "$txtFileDWallTags.sortTemp";
		}
		print STDERR "$txtFileDWallTags\n";
		if ($filterMatches) {
			unless ($directWrite) {
				system "sort -k" . ($cycles + 2) . "rn $txtFileDWfilter > $txtFileDWfilter.sortTemp";
				system "cat $prefix.dwheader $prefix.myheader $txtFileDWfilter.sortTemp $prefix.dwfooter > $txtFileDWfilter";
				unlink "$txtFileDWfilter.sortTemp";
			}
			print STDERR "$txtFileDWfilter\n";
		}
		foreach (".myheader", ".dwheader", ".dwfooter") { unlink "$prefix$_" if -e "$prefix$_" }
	}
}

} # end of the tag files loop

# print stats
if ($totalReads) {
	my $end = time - $start;
	open STATS, ">$logFile" or die "could not write $logFile: $!";
	print STATS "Reads in $prefix\n";
	print STATS "\n";
	print STATS "Total:   $totalReads\n";
	print STATS "Valid:   $validReads\n" if $validReads;
	print STATS "Almost:  $maxLengthSurround\n" if $validReads and $maxLengthSurround;
	print STATS "MaxTagLength: $maxTagLength\n" if $validReads and $maxTagLength and $maxLengthCount != $validReads and not $similarTags and $tagFileCounter == 1;
	print STATS "MaxLengthCount: $maxLengthCount\n" if $validReads and $maxLengthCount and $maxLengthCount != $validReads and not $similarTags and $tagFileCounter == 1;
	print STATS "Opened:  $openedReads\n" if $openedReads;
	print STATS "Forward: $forwardReads\n" if $forwardReads;
	print STATS "Reverse: $reverseReads\n" if $reverseReads;
	print STATS "Similar: $similarReads\n" if $similarReads;
	print STATS "Shorter: $shorterReads\n" if $shorterReads;
	print STATS "Longer:  $longerReads\n" if $longerReads;
	print STATS "Reduced: $reducedReads\n" if $reducedReads;
	print STATS "LowQual: $lowQualReads\n" if $lowQualReads;
	print STATS "Invalid: $invalidReads\n" if $invalidReads;
	print STATS "Recover: $recoverReads\n" if $recoverReads;
	print STATS "\n";
	if ($similarTags) {
		my $errorsDel = $similarCount{"del"} ? $similarCount{"del"} : 0;
		print STATS "Deletions:  $errorsDel\n";
		my $errorsIns = $similarCount{"ins"} ? $similarCount{"ins"} : 0;
		print STATS "Insertions: $errorsIns\n";
		my $errorsVar = $similarCount{"var"} ? $similarCount{"var"} : 0;
		print STATS "1bVariants: $errorsVar\n";
		print STATS "\n";
	}
	print STATS "Chimera: $chimeraReads\n" if $chimeraReads;
	print STATS "Unfound: $unfoundReads\n" if $unfoundReads;
	print STATS "Undedup: $undedupReads\n" if $undedupReads;
	print STATS "Matched: $matchedReads\n";
	print STATS "Deduped: $dedupedReads\n";
	print STATS "Matched recovered: $matchedReadsRecovered\n" if $matchedReadsRecovered;
	print STATS "\n";
	foreach my $closingPrimerId (sort keys %librarySizeCP) {
		my $closingPrimerIdString = "";
		if ($closingPrimerId) {
			$closingPrimerIdString = "$closingPrimerId ";
		} elsif (scalar keys %librarySizeCP > 1) {
			$closingPrimerIdString = "Empty_CP ";
		}
		print STATS $closingPrimerIdString . "Expected uniq: $expectedMatchUniq{$closingPrimerId}\n" if $expectedMatchUniq{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected uniq: $nonexpectedMatchUniq{$closingPrimerId}\n" if $nonexpectedMatchUniq{$closingPrimerId};
		print STATS $closingPrimerIdString . "Expected counts: $expectedMatchCounts{$closingPrimerId}\n" if $expectedMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected counts: $nonexpectedMatchCounts{$closingPrimerId}\n" if $nonexpectedMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Expected dedup counts: $expectedDedupMatchCounts{$closingPrimerId}\n" if $expectedDedupMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected dedup counts: $nonexpectedDedupMatchCounts{$closingPrimerId}\n" if $nonexpectedDedupMatchCounts{$closingPrimerId};
		my $librarySizeString = 0;
		if ( $librarySizeCP{$closingPrimerId} ) {
			$librarySizeString = "$librarySizeCP{$closingPrimerId} (";
			foreach my $cycle (sort {$a<=>$b} keys %{ $libraryTagsCPstats{$closingPrimerId} }) {
				$librarySizeString .= "$libraryTagsCPstats{$closingPrimerId}{$cycle} x ";
			}
			$librarySizeString =~ s/ x $/)/;
		}
		print STATS $closingPrimerIdString . "Library size: $librarySizeString\n";
		print STATS "\n";
	}
	print STATS "Command:\n$commandString\n";
	print STATS "\n";
	print STATS "Running time: $end seconds\n";
	close STATS;
	print STDERR "$logFile\n" unless $threaded;
	print STDERR "Total running time: $end seconds\n" unless $threaded;
}

# end
exit 0;

sub printHeaderFooter {
	my ($prefix, $cycles, $directWrite) = @_;
	my $dwheader = '<datawarrior-fileinfo>
<version="3.1">
</datawarrior-fileinfo>
';
	my $dwfooter = '<datawarrior properties>
<chartType_DEDUP="scatter">
<chartType_RAW="scatter">
<colorColumn_DEDUP="DEDUP">
<colorColumn_RAW="RAW">
<colorListMode_DEDUP="HSBLong">
<colorListMode_RAW="HSBLong">
<color_DEDUP_0="-65536">
<color_DEDUP_1="-16776961">
<color_RAW_0="-65536">
<color_RAW_1="-16776961">
<mainSplitting="0.86587">
<mainView="DEDUP">
<mainViewCount="3">
<mainViewDockInfo0="root">
<mainViewDockInfo1="Table	bottom	0.09">
<mainViewDockInfo2="RAW	right	0.5">
<mainViewName0="Table">
<mainViewName1="RAW">
<mainViewName2="DEDUP">
<mainViewType0="tableView">
<mainViewType1="3Dview">
<mainViewType2="3Dview">
<masterView_DEDUP="RAW">
<rightSplitting="0.68844">
<rotationMatrix_RAW00="-0.87159">
<rotationMatrix_RAW01="-0.16133">
<rotationMatrix_RAW02="0.46298">
<rotationMatrix_RAW10="0.49019">
<rotationMatrix_RAW11="-0.26747">
<rotationMatrix_RAW12="0.82962">
<rotationMatrix_RAW20="-0.010004">
<rotationMatrix_RAW21="0.94997">
<rotationMatrix_RAW22="0.31218">
<markersize_DEDUP="0.6">
<markersize_RAW="0.6">
<sizeColumn_DEDUP="DEDUP">
<sizeColumn_RAW="RAW">
<sizeProportional_DEDUP="true">
<sizeProportional_RAW="true">
<shapeColumn_DEDUP="EXPECTED">
<shapeColumn_RAW="EXPECTED">
<suppressGrid_DEDUP="true">
<suppressGrid_RAW="true">
<fastRendering_DEDUP="true">
<fastRendering_RAW="true">
';
	my $myheader = "";
	my $i = 0;
	for ($i = 0; $i < $cycles; $i++) {
		$myheader .= "TAG" . ($i + 1) . "\t";
		$dwfooter .= "<axisColumn_RAW_$i=\"TAG" . ($i + 1) . "\">\n";
		$dwfooter .= "<filter$i=\"#string#\tTAG" . ($i + 1) . "\">\n";
	}
	$myheader .= "CP\tRAW\tDEDUP\tSTRANDBIAS\tRAW_NORM\tDEDUP_NORM\tEXPECTED$tagsOverHeader\n";
	my $j = $i - 1;
	$j++; $dwfooter .= "<filter$j=\"#string#\tCP\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tRAW\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tDEDUP\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSTRANDBIAS\">\n";
	$j++; $dwfooter .= "<filter$j=\"#category#\tEXPECTED\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSDCOUNT_RAW\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSDCOUNT_DEDUP\">\n";
	foreach my $overString (@overStrings) {
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_RAW_$overString\">\n";
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_DEDUP_$overString\">\n";
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_UNIQUE_$overString\">\n";
	}
	$dwfooter .= "</datawarrior properties>\n";
	unless ($directWrite) {
		open OUT, ">$prefix.myheader" or die;
		print OUT $myheader; close OUT;
		open OUT, ">$prefix.dwheader" or die;
		print OUT $dwheader; close OUT;
		open OUT, ">$prefix.dwfooter" or die;
		print OUT $dwfooter; close OUT;
	}
	return $myheader, $dwheader, $dwfooter;
}

sub split_and_fork {
	my ($fastqFile, $prefix, $threads) = @_;
	my $start = time;
	my $maxThreads = `getconf _NPROCESSORS_ONLN`; # or `nproc --all`, or `grep -c ^processor /proc/cpuinfo`
	if ($threads > $maxThreads) {
		print STDERR "\nThreads limited to $maxThreads\n";
		$threads = $maxThreads;
	}
	if (-f $fastqFile) {
		open FASTQ, "seqtk seq $fastqFile |"
		or open FASTQ, "gunzip -fc $fastqFile |"
		or die "could not open $fastqFile";
	} else {
		die "could not find $fastqFile";
	}
	my $splitFile = "";
	my @splitFiles = ();
	my @filehandles = ();
	my @pids = ();
	foreach my $thread (1..$threads) {
		local *FH;
		$splitFile = "$fastqFile.$thread";
		# open FH, "| gzip >$splitFile" or die "could not create $splitFile";
		open FH, ">$splitFile" or die "could not create $splitFile";
		push @splitFiles, $splitFile;
		push @filehandles, *FH;
	}
	my $reads2split = 0;
	while (my $lines = <FASTQ>) {
		$reads2split++;
		for (1..3) { $lines .= <FASTQ> }
		my $filehandle = $filehandles[$reads2split % $threads];
		print $filehandle $lines;
	}
	close FASTQ;
	foreach my $filehandle (@filehandles) { close $filehandle }
	$elapsed = time - $start;
	print STDERR "Elapsed time post-split: $elapsed seconds\n";
	foreach my $splitFile (@splitFiles) {
		(my $fileCommand = $command) =~ s/$fastqFile/$splitFile/;
		my $childpid = fork() or exec $fileCommand;
		push @pids, $childpid;
	}
	waitpid($_,0) for @pids;
}

# http://www.perlmonks.org/?node=Levenshtein%20distance%3A%20calculating%20similarity%20of%20strings
sub levenshtein {
	my ($s1, $s2) = @_;
	my ($len1, $len2) = (length $s1, length $s2);
	return $len2 if ($len1 == 0);
	return $len1 if ($len2 == 0);
	my %mat;
	for (my $i = 0; $i <= $len1; ++$i) {
		for (my $j = 0; $j <= $len2; ++$j) {
			$mat{$i}{$j} = 0;
			$mat{0}{$j} = $j;
		}
		$mat{$i}{0} = $i;
	}
	my @ar1 = split(//, $s1);
	my @ar2 = split(//, $s2);
	for (my $i = 1; $i <= $len1; ++$i) {
		for (my $j = 1; $j <= $len2; ++$j) {
			my $cost = ($ar1[$i-1] eq $ar2[$j-1]) ? 0 : 1;
			$mat{$i}{$j} = min([
				$mat{$i-1}{$j} + 1,
				$mat{$i}{$j-1} + 1,
				$mat{$i-1}{$j-1} + $cost]);
		}
	}
	return $mat{$len1}{$len2};
}
sub min {
	my @list = @{$_[0]};
	my $min = $list[0];
	foreach my $i (@list) {
		$min = $i if ($i < $min);
	}
	return $min;
}

sub minSeqLD {
	my ($seq1, $len1, $seq2, $len2, $maxIndel) = @_;
	unless (defined $maxIndel) { $maxIndel = $len1 < $len2 ? $len1 : $len2 }
	my $minSeqLD = 0;
	for (my $extraLen = 0; $extraLen <= $maxIndel; $extraLen++) {
		my $ld = 0;
		if ($extraLen) {
			my $lddel = levenshtein($seq1, $seq2 . substr($seq1, $len1 - $extraLen, $extraLen));
			my $ldins = levenshtein($seq1, substr($seq2, 0, $len2 - $extraLen));
			$ld = $lddel < $ldins ? $lddel : $ldins;
			last if $ld > $minSeqLD;
		} else {
			$ld = levenshtein($seq1, $seq2);
		}
		$minSeqLD = $ld;
	}
	return $minSeqLD;
}
