use strict;
use Getopt::Long;
	
	my $scaffold_file = "";
	my $scaffold_filling_file = "";
	my $reference_file = "";
	my $output_directory = "";
	my $min_contig_length = 100;
	my $min_contig_end_distance = 200;
	my $max_distance_normal_gap = 3000;
	my $times_gap_distance = 3;
	
	
	GetOptions('s=s'=>\$scaffold_file,'f=s'=>\$scaffold_filling_file,'r=s'=>\$reference_file,'c=i'=>\$min_contig_length,'m=i'=>\$min_contig_end_distance,'d=i'=>\$max_distance_normal_gap,'t=i'=>\$times_gap_distance,'o=s'=>\$output_directory)
	or die("Error in command line arguments\n");
	
	if(-e $output_directory){
		unlink glob "$output_directory/*";
	}else{
		mkdir($output_directory);
	}
	
	my $contig_set_file = "$output_directory/contig_set.fa";
	my @temp;
	@temp = ("./getContigInScaffold $scaffold_file $contig_set_file $min_contig_length");
	system(@temp) == 0 or die "When running getContigInScaffold, there is an error occuring. Please check the file $scaffold_file.";
	
	my $nucmer_out_prefix = "$output_directory/referece.nucmer";
	my $nucmer_out_delta = $nucmer_out_prefix.".delta";
	my $nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
	my $reference_nucmer_out_coords = $nucmer_out_filter.".coords";
	@temp = ("nucmer --maxmatch -l 90 -p $nucmer_out_prefix $reference_file $contig_set_file 2>&-");
	system(@temp) == 0 or die "Please check whether MUMer has been installed, and check the file $reference_file.";
	@temp = ("delta-filter -i 97 -q $nucmer_out_delta > $nucmer_out_filter");
	`@temp`;
	@temp = ("show-coords -dTlro $nucmer_out_filter > $reference_nucmer_out_coords");
	`@temp`;
	unlink glob $nucmer_out_delta;
	unlink glob $nucmer_out_filter;
	
	$nucmer_out_prefix = "$output_directory/scaffold_filling.nucmer";
	$nucmer_out_delta = $nucmer_out_prefix.".delta";
	$nucmer_out_filter = $nucmer_out_prefix.".delta-filter";
	my $scaffold_filling_nucmer_out_coords = $nucmer_out_filter.".coords";
	@temp = ("nucmer -f --maxmatch -l 90 -p $nucmer_out_prefix $scaffold_filling_file $contig_set_file 2>&-");
	system(@temp) == 0 or die "Please check whether MUMer has been installed, and check the file $scaffold_filling_file.";
	@temp = ("delta-filter -i 97 $nucmer_out_delta > $nucmer_out_filter");
	`@temp`;
	@temp = ("show-coords -dTlro $nucmer_out_filter > $scaffold_filling_nucmer_out_coords");
	`@temp`;
	unlink glob $nucmer_out_delta;
	unlink glob $nucmer_out_filter;
	
	my $gap_in_reference_file = "$output_directory/gap_in_reference.txt";
	my $gap_in_filling_file = "$output_directory/gap_in_filling.txt";
	my $result_file = "$output_directory/result.txt";
	
	my $gap_infor = "$output_directory/gap_infor.txt";
	@temp = ("./getGapRegion $scaffold_file $reference_nucmer_out_coords $reference_file $gap_in_reference_file $scaffold_filling_nucmer_out_coords $scaffold_filling_file $gap_in_filling_file $gap_infor $min_contig_length $min_contig_end_distance $max_distance_normal_gap $times_gap_distance");
	system(@temp) == 0 or die print OUT "When running getGapRegion, there is an error occuring.\n";

	unlink glob $reference_nucmer_out_coords;
	unlink glob $scaffold_filling_nucmer_out_coords;
	
	@temp = ("./needleman_wunsch $gap_in_reference_file $gap_in_filling_file $result_file");
	system(@temp) == 0 or die "When running needleman_wunsch, there is an error occuring.";
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
