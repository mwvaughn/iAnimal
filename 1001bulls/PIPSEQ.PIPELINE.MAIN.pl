#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $readType = ""; #variable to capture read type input
my $input1 = ""; #variable to capture first read file
my $input2 = ""; #variable to capture second read file, if there is one
my $reference = ""; #variable to capture the reference to use for alignment
my $platform = ""; #variable to capture platform information for picard sofware
my $animalNumber = ""; #variable to capture animal id information for picard software
my $coverage = 0; #initialization of the coverage variable for capturing the calculated coverage
my $bwa = ""; #pathway to bwa software
my $samtools = ""; #pathway to samtools software
my $aorrg = ""; #pathway to AddOrReplaceReadGroups.jar file of the picard suite
my $gatk = ""; #pathway to GenomeAnalysisTK.jar file of the GenomeAnalysisTK suite
my @prestack; #stack for pre-coverage file read in system calls
my @poststack; #stack for post-coverage file read in system calls
my $time = time; #setting of unique portion of temporary file names

#creation of variables with temporary and output file names
my $temp1sai = $time.".temp1.sai"; 
my $temp2sai = $time.".temp2.sai";
my $outsam = "bwa_output.".$time.".sam";
my $outbam = "bwa_output.".$time.".bam";
my $outsorted = "bwa_output.".$time.".sorted.bam";
my $picardbam = "picard_output.".$time.".bam";
my $precoverage = $time.".coverage";
my $coveragefile = $precoverage.".sample_summary";
my $gatkSNP = "GATK.SNPS.".$time.".vcf";
my $gatkINDEL = "GATK.INDELS.".$time.".vcf";

#import the user imput and assign to variables
GetOptions('readtype|t=s' => \$readType, 'reference|r=s' => \$reference, 'input1|1=s' => \$input1, 'input2|2=s' => \$input2, 'platform|p=s' => \$platform, 'animalnumber|n=s' => \$animalNumber);
	
push(@prestack, "$bwa aln $reference $input1 > $temp1sai"); #add the first bwa command to the stack

if($readType eq "SE")
{
	push(@prestack, "$bwa samse $reference $temp1sai $input1 > $outsam"); #add the final bwa command to the stack if single-end reads used
}
elsif($readType eq "PE")
{
	push(@prestack, "$bwa aln $reference $input2 > $temp2sai"); #add the next bwa command to the stack if paired-end reads used
	push(@prestack, "$bwa sampe $reference $temp1sai $temp2sai $input1 $input2 > $outsam"); #add the final bwa command to the stack if paired-end reads used
}
		
push(@prestack, "$samtools view -bT $reference $outsam > $outbam"); #add the samtools conversion from sam to bam to the stack
push(@prestack, "$samtools sort $outbam bwa_output.sorted"); #add the samtools sort of bam to the stack
push(@prestack, "$samtools index $outsorted"); #add the samtools index of bam to the stack
push(@prestack, "java -jar $aorrg INPUT=$outsorted OUTPUT=$picardbam RGLB=1 RGPL=$platform RGPU=allruns RGSM= VALIDATION_STRINGENCY=SILENT"); #add the picard AddOrReplaceReadGroups.jar command to the stack
push(@prestack, "$samtools index $picardbam"); #add the samtools index of the corrected bam to the stack
push(@prestack, "java -jar $gatk -R $reference -T DepthOfCoverage -o $precoverage -I $picardbam --omitDepthOutputAtEachBase --omitIntervalStatistics --omitLocusTable"); #add the GenomeAnalysisTK.jar coverage calculation to the stack

#loop through and run the pre-coverage commands
foreach my $cmd (@prestack)
{
	report($cmd);
	system($cmd);
}
				
open (GENE, "<$coveragefile"); #open and loop through the coverage file
while(my $finput = <GENE>)
{
  	chomp $finput;
  	@inputArray = split(/\n/, $finput);
  	foreach $inputLine(@inputArray)
  	{
    	@lineSplit = split(/\t/, $inputLine);
    	if($lineSplit[0] eq "Total")
    		{$coverage = int($lineSplit[2]);} #find and assign the coverage calculation to the coverage variable
  	}
}
			
if($coverage >= 10) #if the coverage calculated is higher than 10 use the following commands
{
	#adds the SNP and the INDEL commands to the stack
	push(@poststack, "java -jar $gatk -R $reference -T UnifiedGenotyper -I $picardbam -o $gatkSNP -stand_call_conf 30.0 -stand_emit_conf 0.0 -dcov $coverage");
	push(@poststack, "java -jar $gatk -R $reference -T UnifiedGenotyper -I $picardbam -o $gatkINDEL -glm INDEL -stand_call_conf 30.0 -stand_emit_conf 0.0 -dcov $coverage");
}
else
{
	#adds the SNP and the INDEL commands to the stack
	push(@poststack, "java -jar $gatk -R $reference -T UnifiedGenotyper -I $picardbam -o $gatkSNP -stand_call_conf 4.0 -stand_emit_conf 0.0 -dcov $coverage");
	push(@poststack, "java -jar $gatk -R $reference -T UnifiedGenotyper -I $picardbam -o $gatkINDEL -glm INDEL -stand_call_conf 4.0 -stand_emit_conf 0.0 -dcov $coverage");
}

#loop through and run the post-coverage commands
foreach my $cmd2 (@poststack)
{
	report($cmd2);
	system($cmd2);
}