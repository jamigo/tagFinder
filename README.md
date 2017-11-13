tagFinder - efficient DEL results analysis
==========================================

``tagFinder`` is a method for fast tag detection and thorough characterization of DNA encoded libraries sequencing results which requires minimal hardware resources, scales linearly, and does not introduce any analytical error. It can even deal with sequencing errors and PCR duplicates on single or double stranded DNA, enhancing the analytical detection and quantification of molecules and the informativeness of the entire process.

Citation
--------
Please cite the following article if you use ``tagFinder`` in your research:
  * Amigo et al., 2017. tagFinder: A Novel Tag Analysis Methodology that Enables Detection of Molecules from DNA-Encoded Chemical Libraries.

Repository files
----------------
  * DEL-A.txt contains tags described by Luca Mannocci in 2009, for example purposes
  * DEL-B.txt contains tags described by Clark et al. in 2009, for example purposes
  * README.md contains this informative content being displayed
  * tagFinder.config.ini contains the schematic format of the data the script requires
  * tagFinder.pl is the main script, to be run on any unix or linux machine

Output: Log file
----------------
Not everything is detecting valid tags from valid reads. Instead, ``tagFinder`` is capable of disecting the input so that not only valid results are perfectly described, but also invalid results are throroughly inspected. It provides a much better understanding of the experiment, which is of great help specially when developing a new library or testing a new technology.
1. General descriptive stats
* **Total**: total number of reads in the experiment
* **Valid**: reads with minimum size and content requirements to be processed
* **Almost**: reads with tag region 1bp more or less than expected (these can be recovered with ``-S`` option)
* **MaxTagLength**: most frequent tag region length (if different from expected)
* **MaxLengthCount**: reads with most frequent tag region length (if different from expected)
* **Opened**: reads with left anchor and no right anchor (with headpiece and no closing primer, if forward)
* **Forward**: reads found in forward direction
* **Reverse**: reads found in reverse direction
* **Similar**: reads containing 1 error (single base change, insertion or deletion)
* **Shorter**: reads without minimum length (tag region + minimum anchoring)
* **Longer**:  reads with bigger tag region than expected (could be chimeras)
* **Reduced**: reads with smaller tag region than expected
* **LowQual**: reads filtered out by base quality
* **Invalid**: reads where tag region could not be found
* **Recover**: reads containing multiple sections of valid data
2. General matching stats
* **Chimera**: reads with unexpectedly repeated tags inside tag region
* **Unfound**: valid reads inspected with no tag match
* **Undedup**: valid reads where degenerated region could not be found (no duplicates check performed)
* **Matched**: valid reads inspected with tag match
* **Deduped**: valid reads inspected with tag match, excluding duplicates
3. Per library stats
* **Expected uniq**: expected unique matches
* **Non expected uniq**: not expected unique matches
* **Expected counts**: total expected matches
* **Non expected counts**: total not expected matches
* **Expected dedup counts**: expected matches, excluding duplicates
* **Non expected dedup counts**: not expected matches, excluding duplicates
* **Library size**: maximum estimation from the number of tags expected per cycle

Output: Over file
-----------------
When it comes to analyze DEL results, finding overrepresented patterns becomes handy. Big beautiful dots standing alone in the resulting plot highlighting single high-affinity compounds would be desireable, but DEL libraries are usually designed to detect sets of compounds based on common characteristics. These can be seen as patterns (lines, planes) in the resulting plot, but they can also be presented as text listing the different tag combinations per library with the number of standard deviations from the average count of all same type of combinations.

Output: Results file
--------------------
Final results are provided in tabular text file that can directly load into [DataWarrior](http://www.openmolecules.org/datawarrior/), a free cheminformatics program for data visualization and analysis. Results include found and not found tags, raw and deduped counts if available, strand bias check if input is double stranded, normalized values by library size, expected and unexpected information, plus pattern discovery through deviations from average.
