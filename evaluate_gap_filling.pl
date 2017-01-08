use strict;
	
	my $scaffold_file = shift;
	my $scaffold_filling_file = shift;
	my $reference_file = shift;
	my $output_directory = shift;
	
	if(-e $output_directory){
		unlink glob "$output_directory/*";
	}else{
		mkdir($output_directory);
	}
	
	my $contig_set_file = "$output_directory/contig_set.fa";
	my @temp;
	@temp = ("./getContigInScaffold $scaffold_file $contig_set_file");
	system(@temp) == 0 or die "When running getContigInScaffold, there is an error occuring. Please check the file $scaffold_file.";
	
	my $nucmer_out_prefix = "$output_directory/referece.nucmer";
	my $nucmer_out_delta = $nucmer_out_prefix.".delta";
	my $nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
	my $reference_nucmer_out_coords = $nucmer_out_filter.".coords";
	@temp = ("nucmer -c 50 -p $nucmer_out_prefix $reference_file $contig_set_file 2>&-");
	system(@temp) == 0 or die "Please check whether MUMer has been installed, and check the file $reference_file.";
	@temp = ("delta-filter -i 97 -q $nucmer_out_delta > $nucmer_out_filter");
	`@temp`;
	@temp = ("show-coords -dTlro $nucmer_out_filter > $reference_nucmer_out_coords");
	`@temp`;
	#unlink glob $nucmer_out_delta;
	#unlink glob $nucmer_out_filter;
	
	$nucmer_out_prefix = "$output_directory/scaffold_filling.nucmer";
	$nucmer_out_delta = $nucmer_out_prefix.".delta";
	$nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
	my $scaffold_filling_nucmer_out_coords = $nucmer_out_filter.".coords";
	@temp = ("nucmer -f -c 50 -p $nucmer_out_prefix $scaffold_filling_file $contig_set_file 2>&-");
	system(@temp) == 0 or die "Please check whether MUMer has been installed, and check the file $scaffold_filling_file.";
	@temp = ("delta-filter -i 97 $nucmer_out_delta > $nucmer_out_filter");
	`@temp`;
	@temp = ("show-coords -dTlro $nucmer_out_filter > $scaffold_filling_nucmer_out_coords");
	`@temp`;
	#unlink glob $nucmer_out_delta;
	#unlink glob $nucmer_out_filter;
	
	my $gap_in_reference_file = "$output_directory/gap_in_reference.fa";
	my $gap_in_filling_file = "$output_directory/gap_in_filling.fa";
	my $result_file = "$output_directory/result.fa";
	
	my $gap_infor = "$output_directory/gap_infor.txt";
	@temp = ("./getGapRegion $scaffold_file $reference_nucmer_out_coords $reference_file $gap_in_reference_file $scaffold_filling_nucmer_out_coords $scaffold_filling_file $gap_in_filling_file $gap_infor");
	system(@temp) == 0 or die print OUT "When running getGapRegion, there is an error occuring.\n";

	#unlink glob $reference_nucmer_out_coords;
	#unlink glob $scaffold_filling_nucmer_out_coords;
	
	
	@temp = ("./needleman_wunsch $gap_in_reference_file $gap_in_filling_file $result_file");
	system(@temp) == 0 or die "When running needleman_wunsch, there is an error occuring.";
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
