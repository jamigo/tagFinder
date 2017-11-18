#!/usr/bin/perl -w

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
  note1:  define multiple headpieces, overhangs and/or closing primers splitting them by comma (",")
  note2:  valid and invalid tags accepts comma (",") separated lists and regex patterns

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
  -x  number of threads (default: 1) [ NOT FULLY FUNCTIONAL ]
  -h  this (help) message

EOF
exit;
}

$start = time;
$commandOptions .= "$_ " foreach (@ARGV);
if ($commandOptions) {
	$command = "$0 $commandOptions";
	($commandString = $command) =~ s/ -/ \\\n-/g;
}

# get options
use Getopt::Std;
my $opt_string = 'f:t:h:o:p:a:q:lsSiNODIXLEcewv:V:Wd:rRT:x:h';
getopts("$opt_string", \%opt) or usage();
usage() if $opt{h};
$tagFile = $opt{t} ? $opt{t} : "";
$headpieces = $opt{h} ? $opt{h} : "";
$overhangs = $opt{o} ? $opt{o} : "";
$closingprimers = $opt{p} ? $opt{p} : "";
$fastqFile = $opt{f} ? $opt{f} : "";
$mbq = defined $opt{q} ? $opt{q} : 0;
$anchorSize = $opt{a} ? $opt{a} : 7;
$leftAnchored = $opt{l} ? 1 : 0;
$similarTags = $opt{s} ? 1 : 0;
$similarTagsStrict = $opt{S} ? 1 : 0; if ($similarTagsStrict) { $similarTags = 1 }
$reverseCycles = $opt{i} ? 1 : 0;
$detectDegen = $opt{N} ? 0 : 1;
$printOver = $opt{O} ? 0 : 1;
$cleanDegen = $opt{D} ? 0 : 1; $cleanDegen = 0 unless $detectDegen;
$printInvalid = $opt{I} ? 1 : 0;
$printExistingTags = $opt{X} ? 1 : 0;
$printLengths = $opt{L} ? 1 : 0;
$printErrorTypes = $opt{E} ? 1 : 0;
$printChimeras = $opt{c} ? 1 : 0;
$printExpected = $opt{e} ? 1 : 0;
$printRawCounts = $opt{w} ? 1 : 0;
$validTags = $opt{v} ? $opt{v} : "";
$invalidTags = $opt{V} ? $opt{V} : "";
$invalidTagsOut = $opt{W} ? 1 : 0;
$tags2focus = $opt{d} ? $opt{d} : "";
$recoveryMode = $opt{r} ? 1 : 0;
$printRecovery = $opt{R} ? 1 : 0;
$testReads = $opt{T} ? $opt{T} : 0;
$threads = $opt{x} ? $opt{x} : 1;

# print degen regions counts
$printDegen = 0;

# max number of compounds to sort at output
$outSortLimit = 100000;

usage() unless $fastqFile;

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

if ($threads > 1) {
	if ($fastqFile =~ /\.\d+$/) {
		$threaded = 1;
	} else {
		splitJoin($fastqFile, $prefix, $threads);
	}
}

# fastq qualities
if ($mbq) {
	# !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
	$fastqQuals = '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~';
	$lowQualChars = substr($fastqQuals, 0, $mbq + 1);
	$lowQualRegex = qr/[$lowQualChars]/
}

# file names
($configFile = $0) =~ s/([^\/]+)\.pl$/$1.config.ini/;
($tagsPath = $0) =~ s/bins\/[^\/]+$/tags\//;
$txtFileDWallTags = "tags_$prefix.allTags.dwar";
$txtFileDWfilter = "tags_$prefix.filtered.dwar";
$invalidFile = "tags_$prefix.invalid.txt";
$chimFile = "tags_$prefix.chimeras.txt";
$distFile = "tags_$prefix.$tags2focus.txt";
$overTagsFile = "tags_$prefix.over.txt";
$existingTagsFile = "tags_$prefix.existingtags.txt";
$lengthsFile = "tags_$prefix.lengths.txt";
$degenFile = "tags_$prefix.degen.txt";
$errorTypeFile = "tags_$prefix.errors.txt";
$recoveryFile = "tags_$prefix.recovery.txt";
$expectedFile = "tags_$prefix.expected.txt";
$tagCountsFile = "tags_$prefix.tagcounts.txt";
$logFile = "tags_$prefix.log";

$noFwAnchorTag = "";

$tagsOverHeader = "\tSDCOUNT_RAW\tSDCOUNT_DEDUP";
$tagsOverTabBlanks = "";
if ($printOver) {
	@tagsOverTypes = ("raw", "dedup", "unique");
	@overStrings = ("LINES", "PLANES");
	foreach $overString (@overStrings) {
		foreach $tagsOverType (@tagsOverTypes) {
			$tagsOverHeader .= "\tOVER_".uc($tagsOverType)."_$overString";
			$tagsOverTabBlanks .= "\t0";
		}
	}
}

# read config file
if (open CONFIG, $configFile) {
	<CONFIG>;
	while (<CONFIG>) {
		if (/^[^#;]\S/) {
			s/[\r\n]//g;
			($f, $i, $t, $h, $o, $p, $v, $V) = split "\t";
			if ($fastqFile =~ /$f/i or $f =~ /$fastqFile/i) {
				$reverseCycles = $reverseCycles ? $reverseCycles : $i;
				$tagFile = $tagFile ? $tagFile : $t;
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

usage() unless $tagFile;

# look for libraries in tag file name
if ($tagFile =~ /:/) {
	while ($tagFile =~ /:([^:,\s]+)/g) { $libSelect{$1} = 1 }
	$tagFile =~ s/:.+//;
}

# convert dots to pattern
$validTags =~ s/\./\\./g;
$invalidTags =~ s/\./\\./g;
$checkValidTags = $validTags or $invalidTags ? 1 : 0;

# check for everything needed
@headpiecesArray = split ",", $headpieces;
@overhangsArray = split ",", $overhangs;
@closingprimersArray = split ",", $closingprimers;
@validTagsArray = split ",", $validTags;
@invalidTagsArray = split ",", $invalidTags;

# read tag file
if ( ! -f $tagFile and -f "$tagsPath$tagFile") { $tagFile = "$tagsPath$tagFile" }
open TAGS, $tagFile or die "could not open $tagFile: $!";
while (<TAGS>) {
	if (/^#ID\tSEQUENCE\t/) {
		s/[\r\n]//g;
		@cols = split "\t";
		$checkValidTags = 1;
		for ($i = 2; $i <= $#cols; $i++) {
			$headerIndex{$i} = $cols[$i];
		}
	} elsif (/^CPL/) {
		s/[\r\n]//g;
		@cols = split "\t";
		$closingPrimer = "$cols[0]-$cols[1]";
		($closingPrimerId = $closingPrimer) =~ s/N[^-]+$//;
		for ($i = 2; $i <= $#cols; $i++) {
			if ( $cols[$i] and ( not %libSelect or exists $libSelect{ $headerIndex{$i} } ) ) {
				$libraryCPs{$closingPrimer} = 1;
				$libraryCPId{ $headerIndex{$i} }{$closingPrimerId} = 1;
			}
		}
	} elsif (/(^\S*?([0-9]+)[.-][0-9]+)\t([ACGT]+)/) {
		$tagCode = $1;
		$cycle = $2;
		$tag = $3;
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
				@cols = split "\t";
				for ($i = 2; $i <= $#cols; $i++) {
					if ( $cols[$i] and ( not %libSelect or exists $libSelect{ $headerIndex{$i} } ) ) {
						foreach $closingPrimerId ( keys %{ $libraryCPId{ $headerIndex{$i} } } ) {
							$libraryTagsCP{$closingPrimerId}{$cycle}++;
							$validTagCodesCP{$closingPrimerId}{$tagCode} = 1;
						}
					}
				}
			}
			if ($checkValidTags) {
				$validFlag = 0;
				foreach (@validTagsArray) {
					@cols = split ";";
					$validTag = $cols[$#cols];
					if ($tagCode =~ /$validTag/) {
						if ($#cols) {
							for ($i = 0; $i < $#cols; $i++) {
								$closingPrimerId = $cols[$i];
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
					@cols = split ";";
					$invalidTag = $cols[$#cols];
					if ($tagCode =~ /$invalidTag/) {
						if ($#cols) {
							for ($i = 0; $i < $#cols; $i++) {
								$closingPrimerId = $cols[$i];
								$closingPrimerId =~ s/N[^-]+$//;
								if ( exists $validTagCodesCP{$closingPrimerId}{$tagCode} ) {
									$libraryTagsCP{$closingPrimerId}{$cycle}--;
									delete $validTagCodesCP{$closingPrimerId}{$tagCode};
								} else {
									foreach $validCpId (keys %validTagCodesCP) {
										if ($validCpId =~ /$closingPrimerId/) {
											if ( exists $validTagCodesCP{$validCpId}{$tagCode} ) {
												$libraryTagsCP{$validCpId}{$cycle}--;
												delete $validTagCodesCP{$validCpId}{$tagCode};
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

# check closing primers if needed
if (%libraryCPs) {
	foreach $closingPrimer (@closingprimersArray) {
		die "closing primer $closingPrimer not found in $tagFile" if not exists $libraryCPs{$closingPrimer};
	}
	@closingprimersArray = sort keys %libraryCPs;
}

# calculate the expected size for a tag combination
foreach $cycle (sort {$a<=>$b} keys %cycleLengths) {
	$cycles++;
	push @sortedCycles, $cycle;
	$cycleLength = $cycleLengths{$cycle};
	push @sortedCyclesLengths, $cycleLength;
	$tagsLength += $cycleLength;
}
foreach $overhang (@overhangsArray) {
	$overhangLength = length($overhang);
	push @overhangsLengths, $overhangLength;
	$tagsLength += $overhangLength;
}
$tagsLengthTop = $tagsLength + 1;
$tagsLengthLow = $tagsLength - 1;
$anchoredTagsLength = $tagsLength + $anchorSize;
if ($similarTags) { $anchoredTagsLength++ }

# deal with head pieces
unless ( @headpiecesArray ) { die "headpieces needed" }
foreach $headPiece (@headpiecesArray) {
	$primerNumber++;
	$headPieceNumber++;
	if ($anchorSize > length($headPiece)) {
		die "anchor size must be smaller than " . length($headPiece) . ", cause head piece";
	}
	# store minimum primer size
	$minPrimerLength = length($headPiece) unless $minPrimerLength;
	$minPrimerLength = length($headPiece) if length($headPiece) < $minPrimerLength;
	# extract anchors from headPieces
	($headPieceRev = reverse $headPiece) =~ tr/ACGT/TGCA/;
	$fwAnchor5 = substr($headPiece, length($headPiece) - $anchorSize, length($headPiece));
	$rvAnchor3 = substr($headPieceRev, 0, $anchorSize);
	# store anchor patterns
	push @anchor5regexes, qr/$fwAnchor5(.+)/;
	push @anchor3regexes, qr/(.+?)$rvAnchor3/;
	# store static sequences
	push @staticSeqs, $headPiece;
	push @staticSeqLengths, length($headPiece);
	# deal with similar patterns
	if ($similarTags) {
		$fwAnchor5pattern = "";
		for ($i = length($fwAnchor5) - 1; $i >= 0; $i--) {
			$anchorPattern = $fwAnchor5;
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
foreach $closingPrimer (@closingprimersArray) {
	$primerNumber++;
	# extract id
	if ($closingPrimer =~ /^(\S*-)(\S+)/) {
		$closingPrimerId = $1;
		$closingPrimer = $2;
	} else {
		$closingPrimerId = "";
	}
	# check for degenerated bases on closing primer
	if ($closingPrimer =~ /^([ACGT]+)(N+)(.+)/) {
		$anchorN = $1;
		$staticSeq = $3;
		if ($anchorSize > length($anchorN)) {
			die "anchor size must be smaller than " . length($anchorN) . ", cause non degen region";
		}
		$lengthN = length($2);
		$closingPrimerId .= $anchorN;
		$regexNcp{$closingPrimerId} = qr/$anchorN(.{$lengthN})/;
		if ($maxDegenErrors) {
			if ( $lengthN > $maxDegenErrors + 1 ) { $maxDegenErrors = $lengthN - 1 }
		} else {
			$maxDegenErrors = $lengthN - 1;
		}
	} else {
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
	($closingPrimerRev = reverse $closingPrimer) =~ tr/ACGT/TGCA/;
	$fwAnchor3 = substr($closingPrimer, 0, $anchorSize);
	$rvAnchor5 = substr($closingPrimerRev, length($closingPrimer) - $anchorSize, length($closingPrimer));
	# store anchor patterns
	push @anchor3regexes, qr/(.+?)$fwAnchor3/;
	push @anchor5regexes, qr/$rvAnchor5(.+)/;
	# store static sequences
	($staticSeq = reverse $staticSeq) =~ tr/ACGT/TGCA/;
	push @staticSeqs, $staticSeq;
	push @staticSeqLengths, length($staticSeq);
	# deal with similar patterns
	if ($similarTags) {
		$rvAnchor5pattern = "";
		for ($i = length($rvAnchor5) - 1; $i >= 0; $i--) {
			$anchorPattern = $rvAnchor5;
			substr($anchorPattern, $i, 1) = substr($anchorPattern, $i, 1)."?.";
			$rvAnchor5pattern .= $anchorPattern . "|";
		}
		chop $rvAnchor5pattern;
		push @anchor3regexesSimilar, qr/(.{$tagsLengthLow,$tagsLengthTop})$fwAnchor3/;
		push @anchor5regexesSimilar, qr/($rvAnchor5pattern)(.+)/;
	}
}
push @closingPrimerIds, $noFwAnchorTag;
$minSeqLength = $minPrimerLength + $anchoredTagsLength;

# make sure tag info depends on closing primer
foreach $closingPrimerId (@closingPrimerIds) {
	foreach $cycle (keys %libraryTags) { $libraryTagsCP{$closingPrimerId}{$cycle} += $libraryTags{$cycle} }
	foreach $tagCode (keys %validTagCodes) { $validTagCodesCP{$closingPrimerId}{$tagCode} = 1 }
}

# get library size
foreach $closingPrimerId (keys %libraryTagsCP) {
	$librarySizeCP{$closingPrimerId} = 1;
	foreach $cycle ( keys %{ $libraryTagsCP{$closingPrimerId} } ) {
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

# get reads
$matchedReads = 0;
$dedupedReads = 0;
if ($printInvalid) { open INVALID, ">$invalidFile" or die "could not write $invalidFile: $!" }
open FASTQ, "seqtk seq $fastqFile |"
	or open FASTQ, "gunzip -fc $fastqFile |"
	or die "could not open $fastqFile";
while (<FASTQ>) {
	$seq = <FASTQ>;
	<FASTQ>;
	$qual = <FASTQ>;
	$totalReads++;
	$recoverRead = 0;
	if ($testReads) {
		if ($totalReads > $testReads) {
			print STDERR "\n* TEST RUN: $testReads reads analyzed *\n\n";
			$totalReads--;
			last;
		}
	}
RECOVERY:
	if (length($seq) < $minSeqLength) {
		$shorterReads++ unless $recoverRead;
		next;
	} else {
		$forward = "";
		$tagString = "";
		$staticSeq = "";
		$shortRead = 0;
		$openedRead = 0;
		$tagPosMatch = 0;
		$closingPrimerId = "";
		ANCHORLOOP: for ($i = 0; $i < $primerNumber; $i++) {
			# qr/$fwAnchor5(.+)/ qr/$rvAnchor5(.+)/
			if ($seq =~ $anchor5regexes[$i]) {
				$openedRead = 1;
				$anchoredSeq = $1;
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
					for ($j = $jMin; $j < $jMax; $j++) {
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
			SIMILARLOOP: for ($i = 0; $i < $primerNumber; $i++) {
				# qr/($fwAnchor5pattern)(.+)/ qr/($rvAnchor5pattern)(.+)/
				if ($seq =~ $anchor5regexesSimilar[$i]) {
					$openedRead = 1;
					$anchoredSeq = $2;
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
						for ($j = $jMin; $j < $jMax; $j++) {
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
		$longer = 0;
		$shorter = 0;
		$similar = "";
		$searchChimeras = 0;
		@tagStrings = ();
		@tagStringsErrorPos = ();
		$tagLength = length($tagString);
		$tagLengths{$tagLength}++;
		if ($tagLength < $tagsLength) {
			if ($similarTags and $tagLength == $tagsLength - 1) {
				for ($baseno = $tagLength - 1; $baseno >= 0; $baseno--) {
					foreach ("A", "C", "G", "T") {
						$tempTagString = $tagString;
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
				for ($baseno = $tagLength - 1; $baseno >= 0; $baseno--) {
					$tempTagString = $tagString;
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
		$tagStringIndex = 0;
		TAGSTRINGS: foreach $tagString (@tagStrings) {
			$match = "";
			$matched = 0;
			$chimera = 0;
			$tagStart = 0;
			for ($cycleIndex = 0; $cycleIndex < $cycles; $cycleIndex++) {
				$cycle = $sortedCycles[$cycleIndex];
				$cycleLength = $sortedCyclesLengths[$cycleIndex];
				$tag = substr($tagString, $tagStart, $cycleLength);
				if ($overhang = $overhangsArray[$cycleIndex]) {
					$overhangLength = $overhangsLengths[$cycleIndex];
					$postTag = substr($tagString, $tagStart + $cycleLength, $overhangLength);
				} else {
					$overhang = "";
					$overhangLength = 0;
					$postTag = "";
				}
				if ($searchChimeras) {
					if ($tagCode = $tags{$cycle}{$tag}) {
						$tagRepeats = $tagString =~ s/$tag$overhang//g;
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
					if ($tagCode = $tags{$cycle}{$tag}) {
						if ($checkValidTags and $similar) { next TAGSTRINGS unless exists $validTagCodesCP{$closingPrimerId}{$tagCode} }
						$matched++;
						$match .= "$tagCode-$tag\t";
						$tagStart += $cycleLength + $overhangLength;
						$matchedTags{$cycle}{$tag}++;
					} elsif ($similarTags) {
						$similarTagUnmatch = 1;
						unless ($similarTagsStrict and $similar) {
							for ($tagPos = $cycleLength - 1; $tagPos >= 0; $tagPos--) {
								$tagPattern = $tag;
								substr($tagPattern, $tagPos, 1) = "[ACGT]";
								if ($cycleTags{$cycle} =~ /;($tagPattern);/) {
									$tagCode = $tags{$cycle}{$1};
									if ($checkValidTags) { next unless exists $validTagCodesCP{$closingPrimerId}{$tagCode} }
									$errorPos = $tagStart + $tagPos + 1;
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
					foreach $errorTypePos (split ";", $similar) {
						if ($errorTypePos) {
							($errorType, $errorPos) = split ",", $errorTypePos;
							$errorPos = $tagStringsErrorPos[$tagStringIndex] unless $errorPos;
							$similarTypes{$errorType}{$errorPos}++;
							$similarCount{$errorType}++;
						}
					}
				}
				$matchedCPreads{$closingPrimerId}++;
				if ( $regexN = $regexNcp{$closingPrimerId} ) {
					$cpMatch = "$closingPrimerId;$match";
					if ($detectDegen) {
						if ($forward) {
							$seq4degen = substr($seq, $tagPosMatch + $tagsLengthLow);
						} else {
							($seq4degen = reverse substr($seq, 0, $tagPosMatch)) =~ tr/ACGT/TGCA/;
						}
						if ($seq4degen =~ $regexN) {
							$indexMatches{$cpMatch}{dedups}{$1}++;
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
				$indexMatches{$cpMatch}{count}++;
				if ($forward) {
					$indexMatches{$cpMatch}{strand}++;
				} else {
					$indexMatches{$cpMatch}{strand}--;
				}
			} else {
				$unfoundReads++ unless $recoverRead;
			}
		}
	} else {
		$invalidReads++ unless $recoverRead;
		if ($printInvalid) { print INVALID "invalid\t$seq" unless $recoverRead }
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

$elapsed = time - $start;
print STDERR "Elapsed time post-patterns: $elapsed seconds\n";

# check error rate looking at statis sequences
if ($cleanDegen) {
	foreach $staticSeq (keys %beginningSeqs) {
		$staticSeqLength = length $staticSeq;
		foreach $beginningSeq ( keys %{ $beginningSeqs{$staticSeq} } ) {
			$ld = minSeqLD($staticSeq, $staticSeqLength, $beginningSeq, length $beginningSeq);
			$ldCount = $beginningSeqs{$staticSeq}{$beginningSeq};
			$beginningSeqsCount{$staticSeq} += $ldCount;
			$beginningSeqsLD{$staticSeq}{$ld} += $ldCount;
		}
		foreach $ld ( keys %{ $beginningSeqsLD{$staticSeq} } ) {
			$beginningSeqsBaseError = $beginningSeqsLD{$staticSeq}{$ld} / ( $beginningSeqsCount{$staticSeq} * $staticSeqLength );
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

# get headers and footer
printHeaderFooter($prefix, $cycles) unless $threaded;

# deal with tag lengths
$maxLengthCount = 0;
foreach $tagLength (sort {$a<=>$b} keys %tagLengths) {
	$tagLengthCount = $tagLengths{$tagLength};
	$tagLengthString .= "$tagLength\t$tagLengthCount\n";
	if ($tagLengthCount > $maxLengthCount) {
		$maxTagLength = $tagLength;
		$maxLengthCount = $tagLengthCount;
	}
}
$maxLengthSurround = $tagLengths{$maxTagLength - 1} + $tagLengths{$maxTagLength + 1};
if ($printLengths) {
	open LENGTHS, ">$lengthsFile" or die "could not write $lengthsFile: $!";
	print LENGTHS "#LENGTHS\tCOUNTS\n$tagLengthString";
	close LENGTHS;
	print STDERR "$lengthsFile\n";
}

# print degen counts
if ($printDegen) {
	open DEGEN, "| sort -k1,1 -k3nr >$degenFile" or die "could not write $degenFile: $!";
	print DEGEN "#CP\tDEGENSEQ\tCOUNT\n";
	foreach $closingPrimerId (keys %degenCounts) {
		foreach $degenSeq ( keys %{ $degenCounts{$closingPrimerId} } ) {
			print DEGEN "$closingPrimerId\t$degenSeq\t$degenCounts{$closingPrimerId}{$degenSeq}\n";
		}
	}
	close DEGEN;
	print STDERR "$degenFile\n";
}

# print error types
if ($printErrorTypes) {
	open ERRORTYPE, ">$errorTypeFile" or die "could not write $errorTypeFile: $!";
	foreach $type (sort keys %similarTypes) {
		foreach $pos ( sort {$a<=>$b} keys %{ $similarTypes{$type} } ) {
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
	foreach $recoveryTime (sort {$a<=>$b} keys %recoverReadCount) {
		print RECOVERY "$recoveryTime\t$recoverReadCount{$recoveryTime}\n" if $recoverReadCount{$recoveryTime};
	}
	close RECOVERY;
	print STDERR "$recoveryFile\n";
}

# print chimeras
if ($chimeraReads and $printChimeras) {
	open CHIMERAS, "| sort -k" . ($cycles + 2) . "rn >$chimFile" or die "could not write $chimFile: $!";
	foreach $cpMatch (keys %chimeras) {
		($closingPrimerID, $match) = split ";", $cpMatch;
		print CHIMERAS "$match$closingPrimerID\t$chimeras{$cpMatch}\n";
	}
	close CHIMERAS;
	print STDERR "$chimFile\n";
}

# print raw tag counts
if ($matchedReads and $printRawCounts) {
	open TAGCOUNTS, "| sort -k3nr >$tagCountsFile" or die "could not open $tagCountsFile: $!";
	foreach $cycle ( sort {$a<=>$b} keys %matchedTags ) {
		foreach $tag ( sort keys %{ $matchedTags{$cycle} } ) {
			print TAGCOUNTS "$cycle\t$tag\t$matchedTags{$cycle}{$tag}\t$tags{$cycle}{$tag}\t$prefix\n";
		}
	}
	foreach $cycle ( sort {$a<=>$b} keys %unmatchedTags ) {
		foreach $tag ( sort keys %{ $unmatchedTags{$cycle} } ) {
			print TAGCOUNTS "$cycle\t$tag\t$unmatchedTags{$cycle}{$tag}\t\t$prefix\n";
		}
	}
	close TAGCOUNTS;
	print STDERR "$tagCountsFile\n";
}

# print particular tag counts
if ($matchedReads and $tags2focus and $detectDegen) {
	($tags2focusString = $tags2focus) =~ s/_/\t/g;
	$tags2focusString .= "\t";
	open DIST, "| sort -k3nr >$distFile" or die "could not write $distFile: $!";
	foreach $closingPrimerId (sort keys %matchedCPreads) {
		$cpString = "$closingPrimerId;$tags2focusString";
		if ( exists $indexMatches{$cpString} ) {
			if ( exists $indexMatches{$cpString}{dedups} ) {
				foreach $degenSeq ( keys %{ $indexMatches{$cpString}{dedups} } ) {
					print DIST "$closingPrimerId\t$degenSeq\t$indexMatches{$cpString}{dedups}{$degenSeq}\n";
				}
			}
		}
	}
	close DIST;
	print STDERR "$distFile\n";
}

# get the number of compounds
$cpMatchesCount = keys %indexMatches;

# look for overrepresented tags
if ($matchedReads) {
	foreach $cpMatch (keys %indexMatches) {
		$dedupCount = scalar keys %{ $indexMatches{$cpMatch}{dedups} };
		$dedupCount = 1 unless $dedupCount;
		if ($cleanDegen and $dedupCount > 1 and $dedupCount < 10000) {
			@degenSeqs = sort { $indexMatches{$cpMatch}{dedups}{$b} <=> $indexMatches{$cpMatch}{dedups}{$a} or $a cmp $b } keys %{ $indexMatches{$cpMatch}{dedups} };
			@degenSeqsRev = reverse @degenSeqs;
			pop @degenSeqs;
			pop @degenSeqsRev;
			%degenSeqsRemoved = ();
			foreach $degenSeq1 (@degenSeqs) {
				unless ( exists $degenSeqsRemoved{$degenSeq1} ) {
					$degenSeqLength1 = length $degenSeq1;
					if ($degenSeqLength1 > $maxDegenErrors) {
						for $degenError (1..$maxDegenErrors) {
							if ( exists $beginningSeqsBaseErrors{$degenError} ) {
								$degenSeqErrorProb = $indexMatches{$cpMatch}{dedups}{$degenSeq1} * $degenSeqLength1 * $beginningSeqsBaseErrors{$degenError};
								foreach $degenSeq2 (@degenSeqsRev) {
									unless ( exists $degenSeqsRemoved{$degenSeq2} ) {
										if ($indexMatches{$cpMatch}{dedups}{$degenSeq2} < $degenSeqErrorProb) {
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
		delete $indexMatches{$cpMatch}{dedups} unless $threaded;
		$indexMatches{$cpMatch}{dedup} = $dedupCount;
		$matchCount = $indexMatches{$cpMatch}{count};
		($closingPrimerId, $match) = split ";", $cpMatch;
		if ($printOver) {
			$countsUniq{$closingPrimerId}++;
			$countsAverage{$closingPrimerId}{raw} += $matchCount;
			$countsAverage{$closingPrimerId}{dedup} += $dedupCount;
		}
		@matchedTags = split "\t", $match;
		for ($cycleIndex1 = 0; $cycleIndex1 <= $#matchedTags; $cycleIndex1++) {
			($tag1 = $matchedTags[$cycleIndex1]) =~ s/.+-//;
			$tagsFound{$closingPrimerId}{"$cycleIndex1;$tag1"} = 1;
			if ($printOver) {
				$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{raw} += $matchCount;
				$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{dedup} += $dedupCount;
				$structuresFound{$closingPrimerId}{planes}{"$cycleIndex1;$tag1"}{unique}++;
				for ($cycleIndex2 = $cycleIndex1 + 1; $cycleIndex2 <= $#matchedTags; $cycleIndex2++) {
					($tag2 = $matchedTags[$cycleIndex2]) =~ s/.+-//;
					$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{raw} += $matchCount;
					$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{dedup} += $dedupCount;
					$structuresFound{$closingPrimerId}{lines}{"$cycleIndex1;$tag1;$cycleIndex2;$tag2"}{unique}++;
				}
			}
		}
	}
	if ($printOver) {
		foreach $closingPrimerId (keys %tagsFound) {
			foreach $type ("raw", "dedup") {
				$countsAverage{$closingPrimerId}{$type} /= $countsUniq{$closingPrimerId};
				$tempAverage{$closingPrimerId}{$type} = $countsAverage{$closingPrimerId}{$type};
			}
		}
		foreach $cpMatch (keys %indexMatches) {
			($closingPrimerId, $match) = split ";", $cpMatch;
			$matchCount = $indexMatches{$cpMatch}{count};
			$dedupCount = $indexMatches{$cpMatch}{dedup};
			foreach $type ("raw", "dedup") {
				if ($type eq "raw") { $tempValue = $matchCount } else { $tempValue = $dedupCount }
				$tempStdDev{$closingPrimerId}{$type} += ( $tempAverage{$closingPrimerId}{$type} - $tempValue ) ** 2;
				$tempCount{$closingPrimerId}{$type}++;
			}
		}
		foreach $closingPrimerId (keys %tagsFound) {
			foreach $type ("raw", "dedup") {
				$countsStdDev{$closingPrimerId}{$type} = sqrt( $tempStdDev{$closingPrimerId}{$type} / $tempCount{$closingPrimerId}{$type} );
			}
		}
		open OVER, ">$overTagsFile" or die "could not open $overTagsFile: $!";
		print OVER "#TYPE\tCOUNT\tSDCOUNT\tCP\tTAGS\n";
		foreach $type (@tagsOverTypes) {
			foreach $closingPrimerId (keys %structuresFound) {
				foreach $structure ("planes", "lines") {
					$cpCycleTagInfo = $structuresFound{$closingPrimerId}{$structure};
					@cpCycleTags = keys %$cpCycleTagInfo;
					$cpCycleTagStdDev = 0;
					$cpCycleTagCount = 0;
					$cpCycleTagAverage = 0;
					foreach $cpCycleTag (@cpCycleTags) {
						$cpCycleTagCount++;
						$cpCycleTagAverage += $cpCycleTagInfo->{$cpCycleTag}{$type};
					}
					$cpCycleTagAverage /= $cpCycleTagCount;
					foreach $cpCycleTag (@cpCycleTags) {
						$cpCycleTagStdDev += ( $cpCycleTagAverage - $cpCycleTagInfo->{$cpCycleTag}{$type} ) ** 2;
					}
					$cpCycleTagStdDev = sqrt($cpCycleTagStdDev / $cpCycleTagCount);
					$cpCycleTagThreshold = $cpCycleTagAverage + $cpCycleTagStdDev;
					foreach $cpCycleTag (@cpCycleTags) {
						$tempValue = $cpCycleTagInfo->{$cpCycleTag}{$type};
						if ($tempValue > $cpCycleTagThreshold) {
							($cycleIndex1, $tag1, $cycleIndex2, $tag2) = split ";", $cpCycleTag;
							$tag1Code = $tags{ $sortedCycles[$cycleIndex1] }{$tag1};
							if ($tag2) {
								$tag2Code = $tags{ $sortedCycles[$cycleIndex2] }{$tag2};
								$tagCodeTag = "$tag1Code-$tag1;$tag2Code-$tag2";
							} else {
								$tagCodeTag = "$tag1Code-$tag1";
							}
							$tempCount = 1;
							$tempStdDev = $cpCycleTagThreshold + $cpCycleTagStdDev;
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
		foreach $cpMatch (keys %indexMatches) {
			($closingPrimerId, $match) = split ";", $cpMatch;
			$matchCount = $indexMatches{$cpMatch}{count};
			$dedupCount = $indexMatches{$cpMatch}{dedup};
			$sdString = "";
			foreach $type ("raw", "dedup") {
				if ($type eq "raw") {
					$tempValue = $matchCount;
				} else {
					$tempValue = $dedupCount;
				}
				$tempCount = 0;
				while ($tempValue > $countsAverage{$closingPrimerId}{$type} + ($tempCount + 1) * $countsStdDev{$closingPrimerId}{$type} ) { $tempCount++ }
				$sdString .= "\t$tempCount";
				if ($tempCount) {
					chop $match;
					$match =~ s/\t/;/g;
					print OVER "$type\t$tempValue\t$tempCount\t$closingPrimerId\t$match\n";
				}
			}
			$indexMatches{$cpMatch}{sd} = $sdString;
		}
		close OVER;
		if ($cpMatchesCount < $outSortLimit) {
			system "sort -k1,1 -k3,3nr -k4,4V -k3,3nr $overTagsFile > $overTagsFile.sortTemp";
			rename "$overTagsFile.sortTemp", $overTagsFile;
		}
		print STDERR "$overTagsFile\n";
	}
}

# print expected tags
if ($matchedReads and $checkValidTags) {
	if ($printExpected) {
		open EXPECTED, ">$expectedFile" or die "could not open $expectedFile: $!";
	}
	for ($cycleIndex = 0; $cycleIndex <= $#sortedCycles; $cycleIndex++) {
		$cycle = $sortedCycles[$cycleIndex];
		foreach $tag ( keys %{ $tags{$cycle} } ) {
			$tagCode = $tags{$cycle}{$tag};
			print EXPECTED "$tagCode\t$tag" if $printExpected;
			foreach $closingPrimerId (sort keys %matchedCPreads) {
				$tagFlag = $validTagCodesCP{$closingPrimerId}{$tagCode} ? 1 : 0;
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
		open EXPECTED, ">$expectedFile.sortTemp" or die "could not open $expectedFile.sortTemp: $!";
		print EXPECTED "ID\tTAG";
		foreach $closingPrimerId (sort keys %matchedCPreads) {
			$closingPrimerHeader = $closingPrimerId ? $closingPrimerId : "EMPTY";
			print EXPECTED "\t$closingPrimerHeader";
		}
		print EXPECTED "\n";
		close EXPECTED;
		if ($cpMatchesCount < $outSortLimit) {
			system "sort -k1n $expectedFile >> $expectedFile.sortTemp";
		} else {
			writeFiles($expectedFile, "$expectedFile.sortTemp", 1);
		}
		rename "$expectedFile.sortTemp", $expectedFile;
		print STDERR "$expectedFile\n";
	} else {
		unlink $expectedFile;
	}
}

# deal with matches
if ($matchedReads) {
	open MATCHES, ">$txtFileDWallTags" or die "could not write $txtFileDWallTags: $!";
	foreach $closingPrimerId (@closingPrimerIds) {
		if ( $librarySizeCP{$closingPrimerId} and $matchedCPreads{$closingPrimerId} ) {
			$librarySizeNorm{$closingPrimerId} = $librarySizeCP{$closingPrimerId} / $matchedCPreads{$closingPrimerId};
		} else {
			$librarySizeNorm{$closingPrimerId} = 0;
		}
	}
	foreach $cpMatch (keys %indexMatches) {
		($closingPrimerId, $match) = split ";", $cpMatch;
		$matchCount = $indexMatches{$cpMatch}{count};
		$dedupCount = $indexMatches{$cpMatch}{dedup};
		$strandCount = sprintf( "%.3f", abs( $indexMatches{$cpMatch}{strand} / $matchCount ) );
		$dedupedReads += $dedupCount;
		$cpLibrarySizeNorm = $librarySizeNorm{$closingPrimerId};
		if ($checkValidTags) {
			$expected = 1;
			while ($match =~ /(\S+)-/g) {
				unless ( exists $validTagCodesCP{$closingPrimerId}{$1} ) {
					$expected = 0;
					last;
				}
			}
		} else {
			$expected = 0;
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
		if ($threaded) {
			$degenSeq = join ",", keys %{ $indexMatches{$cpMatch}{dedups} };
			print MATCHES "$match$closingPrimerId\t$matchCount\t$expected\t$degenSeq\n";
		} else {
			if ($printOver) {
				$overTagString = $indexMatches{$cpMatch}{sd};
				@matchedTags = split "\t", $match;
				foreach $type (@tagsOverTypes) {
					$overTagValue = 0;
					for ($cycleIndex1 = 0; $cycleIndex1 < $#matchedTags; $cycleIndex1++) {
						$matchedTag1 = $matchedTags[$cycleIndex1];
						for ($cycleIndex2 = $cycleIndex1 + 1; $cycleIndex2 <= $#matchedTags; $cycleIndex2++) {
							$matchedTag2 = $matchedTags[$cycleIndex2];
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
				foreach $type (@tagsOverTypes) {
					$overTagValue = 0;
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
			} else {
				$overTagString = "";
			}
			$matchCountNorm = sprintf "%.5f", $matchCount * $cpLibrarySizeNorm;
			$dedupCountNorm = sprintf "%.5f", $dedupCount * $cpLibrarySizeNorm;
			print MATCHES "$match$closingPrimerId\t$matchCount\t$dedupCount\t$strandCount\t$matchCountNorm\t$dedupCountNorm\t$expected$overTagString\n";
		}
	}
	close MATCHES;
	if ( ($cpMatchesCount < $outSortLimit) and not $threaded) {
		system "sort -k" . ($cycles + 2) . "rn $txtFileDWallTags > $txtFileDWallTags.sortTemp";
	} else {
		rename $txtFileDWallTags, "$txtFileDWallTags.sortTemp";
	}
	if ($filterMatches) {
		system "cp $txtFileDWallTags.sortTemp $txtFileDWfilter.sortTemp";
	}
	open MATCHES, ">>$txtFileDWallTags.sortTemp" or die "could not write $txtFileDWallTags.sortTemp: $!";
	foreach $closingPrimerId (sort keys %idTagsMissing) {
		foreach $cycleIndex ( sort {$a<=>$b} keys %{ $idTagsMissing{$closingPrimerId} } ) {
			foreach $tagCode ( sort keys %{ $idTagsMissing{$closingPrimerId}{$cycleIndex} } ) {
				$match = "";
				for (0..$cycles-1) {
					if ($cycleIndex == $_) { $match .= "$tagCode\t" } else { $match .= "\t" }
				}
				print MATCHES "$match$closingPrimerId\t0\t0\t0\t0\t0\t$idTagsMissing{$closingPrimerId}{$cycleIndex}{$tagCode}\t0\t0$tagsOverTabBlanks\n";
			}
		}
	}
	close MATCHES;
	if ($threaded) {
		rename "$txtFileDWallTags.sortTemp", $txtFileDWallTags;
		rename "$txtFileDWfilter.sortTemp", $txtFileDWfilter;
	} else {
		writeFiles("$prefix.dwheader", $txtFileDWallTags, 0);
		writeFiles("$prefix.header", $txtFileDWallTags, 1);
		writeFiles("$txtFileDWallTags.sortTemp", $txtFileDWallTags, 1);
		writeFiles("$prefix.dwfooter", $txtFileDWallTags, 1);
		unlink "$txtFileDWallTags.sortTemp";
		print STDERR "$txtFileDWallTags\n";
		if ($filterMatches) {
			writeFiles("$prefix.dwheader", $txtFileDWfilter, 0);
			writeFiles("$prefix.header", $txtFileDWfilter, 1);
			writeFiles("$txtFileDWfilter.sortTemp", $txtFileDWfilter, 1);
			writeFiles("$prefix.dwfooter", $txtFileDWfilter, 1);
			unlink "$txtFileDWfilter.sortTemp";
			print STDERR "$txtFileDWfilter\n";
		}
	}
}

# print stats
if ($totalReads) {
	$end = time - $start;
	open STATS, ">$logFile" or die "could not write $logFile: $!";
	print STATS "Reads in $prefix\n";
	print STATS "\n";
	print STATS "Total:   $totalReads\n" if $totalReads;
	print STATS "Valid:   $validReads\n" if $validReads;
	print STATS "Almost:  $maxLengthSurround\n" if $validReads and $maxLengthSurround;
	print STATS "MaxTagLength: $maxTagLength\n" if $validReads and $maxTagLength and $maxLengthCount != $validReads and not $similarTags;
	print STATS "MaxLengthCount: $maxLengthCount\n" if $validReads and $maxLengthCount and $maxLengthCount != $validReads and not $similarTags;
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
		$errorsDel = $similarCount{"del"} ? $similarCount{"del"} : 0;
		print STATS "Deletions:  $errorsDel\n";
		$errorsIns = $similarCount{"ins"} ? $similarCount{"ins"} : 0;
		print STATS "Insertions: $errorsIns\n";
		$errorsVar = $similarCount{"var"} ? $similarCount{"var"} : 0;
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
	foreach $closingPrimerId (sort keys %matchedCPreads) {
		if ($closingPrimerId) {
			$closingPrimerIdString = "$closingPrimerId ";
		} elsif (scalar keys %matchedCPreads > 1) {
			$closingPrimerIdString = "Empty_CP ";
		} else {
			$closingPrimerIdString = "";
		}
		print STATS $closingPrimerIdString . "Expected uniq: $expectedMatchUniq{$closingPrimerId}\n" if $expectedMatchUniq{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected uniq: $nonexpectedMatchUniq{$closingPrimerId}\n" if $nonexpectedMatchUniq{$closingPrimerId};
		print STATS $closingPrimerIdString . "Expected counts: $expectedMatchCounts{$closingPrimerId}\n" if $expectedMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected counts: $nonexpectedMatchCounts{$closingPrimerId}\n" if $nonexpectedMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Expected dedup counts: $expectedDedupMatchCounts{$closingPrimerId}\n" if $expectedDedupMatchCounts{$closingPrimerId};
		print STATS $closingPrimerIdString . "Non expected dedup counts: $nonexpectedDedupMatchCounts{$closingPrimerId}\n" if $nonexpectedDedupMatchCounts{$closingPrimerId};
		if ( $librarySizeCP{$closingPrimerId} ) {
			$librarySizeString = "$librarySizeCP{$closingPrimerId} (";
			foreach $cycle (sort {$a<=>$b} keys %{ $libraryTagsCP{$closingPrimerId} }) {
				$librarySizeString .= "$libraryTagsCP{$closingPrimerId}{$cycle} x ";
			}
			$librarySizeString =~ s/ x $/)/;
		} else {
			$librarySizeString = 0;
		}
		print STATS $closingPrimerIdString . "Library size: $librarySizeString\n";
		print STATS "\n";
	}
	print STATS "Command:\n$commandString\n";
	print STATS "\n";
	print STATS "Running time: $end seconds\n";
	close STATS;
	print STDERR "$logFile\n";
	print STDERR "Total running time: $end seconds\n" unless $threaded;
}

# end
unlink "$prefix.header";
unlink "$prefix.dwheader";
unlink "$prefix.dwfooter";
exit;

sub writeFiles {
	my ($inFile, $outFile, $append) = @_;
	if ($append) { $appendString = ">>" } else { $appendString = ">" }
	open IN, $inFile or die "could not read $inFile: $!";
	open OUT, $appendString, $outFile or die "could not write $outFile: $!";
	while (<IN>) { print OUT } close OUT; close IN;
}

sub printHeaderFooter {
	my ($prefix, $cycles) = @_;
	$dwheader = '<datawarrior-fileinfo>
<version="3.1">
</datawarrior-fileinfo>
';
	$dwfooter = '<datawarrior properties>
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
	$header = "";
	for ($i = 0; $i < $cycles; $i++) {
		$header .= "TAG" . ($i + 1) . "\t";
		$dwfooter .= "<axisColumn_RAW_$i=\"TAG" . ($i + 1) . "\">\n";
		$dwfooter .= "<filter$i=\"#string#\tTAG" . ($i + 1) . "\">\n";
	}
	$header .= "CP\tRAW\tDEDUP\tSTRANDBIAS\tRAW_NORM\tDEDUP_NORM\tEXPECTED$tagsOverHeader\n";
	$j = $i - 1;
	$j++; $dwfooter .= "<filter$j=\"#string#\tCP\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tRAW\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tDEDUP\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSTRANDBIAS\">\n";
	$j++; $dwfooter .= "<filter$j=\"#category#\tEXPECTED\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSDCOUNT_RAW\">\n";
	$j++; $dwfooter .= "<filter$j=\"#double#\tSDCOUNT_DEDUP\">\n";
	foreach $overString (@overStrings) {
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_RAW_$overString\">\n";
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_DEDUP_$overString\">\n";
		$j++; $dwfooter .= "<filter$j=\"#category#\tOVER_UNIQUE_$overString\">\n";
	}
	$dwfooter .= "</datawarrior properties>\n";
	open OUT, ">$prefix.header" or die;
	print OUT $header; close OUT;
	open OUT, ">$prefix.dwheader" or die;
	print OUT $dwheader; close OUT;
	open OUT, ">$prefix.dwfooter" or die;
	print OUT $dwfooter; close OUT;
}

sub splitJoin {
	my ($fastqFile, $prefix, $threads) = @_;
	$start = time;
	$maxThreads = `getconf _NPROCESSORS_ONLN`; # or `nproc --all`, or `grep -c ^processor /proc/cpuinfo`
	if ($threads > $maxThreads) {
		print STDERR "\nThreads limited to $maxThreads\n";
		$threads = $maxThreads;
	}
	open FASTQ, "seqtk seq $fastqFile |"
		or open FASTQ, "gunzip -fc $fastqFile |"
		or die "could not open $fastqFile";
	foreach $thread (1..$threads) {
		local *FH;
		$splitFile = "$fastqFile.$thread";
		# open(FH, "| gzip >$splitFile") || die "could not create $splitFile";
		open(FH, ">$splitFile") || die "could not create $splitFile";
		push @splitFiles, $splitFile;
		push @filehandles, *FH;
	}
	while ($lines = <FASTQ>) {
		$totalReads++;
		for (1..3) { $lines .= <FASTQ> }
		$filehandle = $filehandles[$totalReads % $threads];
		print $filehandle $lines;
	}
	close FASTQ;
	foreach $filehandle (@filehandles) { close $filehandle }
	$elapsed = time - $start;
	print STDERR "Elapsed time post-split: $elapsed seconds\n";
	# use Parallel::ForkManager;
	# $pm = new Parallel::ForkManager($threads);
	foreach $splitFile (@splitFiles) {
		($fileCommand = $command) =~ s/$fastqFile/$splitFile/;
		my $childpid = fork() or exec $fileCommand;
		push @pids, $childpid;
		# my $pid = $pm->start and next;
		# $forkResults{$splitFile} = `$fileCommand`;
		# $pm->finish;
	}
	# $pm->wait_all_children;
	# foreach (@splitFiles) { print $forkResults{$_} }
	waitpid($_,0) for @pids;
	foreach (@splitFiles) { unlink }
	foreach (< tags_$prefix.*.log >) { unlink }
	foreach (< tags_$prefix.*.dwar >) { unlink }
	$end = time - $start;
	print STDERR "Total running time: $end seconds\n";
	exit;
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
	for ($extraLen = 0; $extraLen <= $maxIndel; $extraLen++) {
		if ($extraLen) {
			$lddel = levenshtein($seq1, $seq2 . substr($seq1, $len1 - $extraLen, $extraLen));
			$ldins = levenshtein($seq1, substr($seq2, 0, $len2 - $extraLen));
			$ld = $lddel < $ldins ? $lddel : $ldins;
			last if $ld > $minSeqLD;
		} else {
			$ld = levenshtein($seq1, $seq2);
		}
		$minSeqLD = $ld;
	}
	return $minSeqLD;
}
