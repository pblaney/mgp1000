#!/usr/bin/env perl

# This script was developed by Igor Dolgalev https://github.com/igordot/sns/blob/main/gather-fastqs
# It has been adapted here with modifications.

use strict;
use warnings;

main();

sub main {

	# Set user-defined input directory as directory to search
	my $input_dir = $ARGV[0];
	chomp($input_dir);

	# Find all the FASTQs in the input directory
	my $fastq_name_globs = "-name '*_R1_0*.fastq.gz' -or -name '*_R1.fastq.gz'";
	my $fastq_find_command = "find -L $input_dir -maxdepth 8 -type f $fastq_name_globs | LC_ALL=C sort";
	my @input_fastqs = `$fastq_find_command`;

	# Output sample file with paired FASTQs
	my $output_filename = $ARGV[1];
	open(my $output_file, ">>", $output_filename);

	# Add header to output file
	print $output_file "sample_id\tread_1\tread_2\n";

	# Using all found R1 FASTQs in the input directory, process each pair
	while (my $read1_fastq = shift(@input_fastqs)) {
		chomp($read1_fastq);

		# Generate the R2 FASTQ filename based on R1
		my $read2_fastq = $read1_fastq;
		$read2_fastq =~ s/(.*)_R1_0([0-9]+.fastq.gz)/${1}_R2_0${2}/;
		$read2_fastq =~ s/(.*)_R1.fastq.gz/${1}_R2.fastq.gz/;

		# Remove directory structure from reads
		$read1_fastq =~ s/.*\///;
		$read2_fastq =~ s/.*\///;

		# Extract sample ID to be used as identifier for downstream processes
		my $sample_id = $read1_fastq;
		# First remove the trailing .fastq.gz
		$sample_id =~ s/.fastq.gz//;

		# bcl2fastq2 format (with S sample number)
		$sample_id =~ s/_S[0-9]{1,4}_L00[0-9]_R1.*//;
		# bcl2fastq format with 2 barcodes
		$sample_id =~ s/_[ACTG]{6,}-[ACTG]{6,}_L00[0-9]_R1.*//;
		# bcl2fastq format with 1 barcode
		$sample_id =~ s/_[ACTG]{4,}_L00[0-9]_R1.*//;
		# no barcodes
		$sample_id =~ s/_L00[0-9]_R[12].*//;
		# no barcodes or lane
		$sample_id =~ s/_R[12]//;

		my $output_line = "$sample_id\t$read1_fastq\t$read2_fastq\n";
		print $output_file "$output_line";
	}

	close($output_file);
}

