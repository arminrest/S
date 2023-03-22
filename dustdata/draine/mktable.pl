#!/usr/bin/perl

# This script creates the 2-dimensional tables <dusttype>.<parameter>.2dtab for
# <dusttype> in ['Gra_81','PAHion_30','PAHneu_30','suvSil_81'] from Draine
# http://www.astro.princeton.edu/~draine/dust/dust.diel.html
# <parameter> in [' Q_ext','Q_abs','Q_sca','g']
# x-dimension: dust grain radius a in microns
# y-dimension: wavelength in microns
# values: <parameter>(x,y)

push(@INC, $ENV{"PIPE_PERL"});
require "genericprocs.pl";

my %statustrings = (-1 => '?', 0 => 'BAD', 1 => '#', 2 => '-');

%colnames = (1=>'Q_ext',
	     2=>'Q_abs',
	     3=>'Q_sca',
	     4=>'g');

#&mktable("Gra_81");
#exit(0);
 
&mktable("PAHion_30");
&mktable("PAHneu_30");
%colnames = (1=>'Q_abs',
	     2=>'Q_sca',
	     3=>'g',
	     4=>'delme');
&mktable("suvSil_81");
&mktable("Gra_81");

sub mktable{
    my ($filename) =@_;
    print STDERR "$filename\n";
    my @lines = &LoadFromFile($filename);
    my $istart = -1;
    for ($i=0;$i<@lines;$i++){
	if ($lines[$i] =~ /radius\(micron\)/){
	    $istart = $i;
	    last;
	}
    }
    if ($istart<0){
	print STDERR "ERROR: Could not determine istart!\n";
	exit(0);
    }
    print STDERR "istart: $istart\n";
    
    undef my %output;
    for (my $j=1;$j<=4;$j++){
	$output[$j]=[];
    }
    my @lambda = [];

    my $a_micron = 0.0;
    my $a_counter = 0;
    my $lambdacounter = 0;
    my $Nlambda = -1;
    for (my $i=$istart;$i<@lines;$i++){
	my $line = $lines[$i];
	next if ($line eq "");
	next if ($line =~ /^w/);
	$line =~ s/^\s+//;
	if ($line =~ /radius\(micron\)/){
	    ($a_micron) = $line =~ /(\S+)\s+\=\s+radius\(micron\)/;
	    if ($a_micron eq ""){
		print STDERR "ERROR: Could not determine a_micron ($line)!\n";
		exit(0);		
	    }
	    $a_counter ++;
	    # just to be sure...
	    if ($a_counter>2){
		if ($Nlambda != $lambdacounter){
		    print STDERR "Bug! inconsistent number of lambdas ($Nlambda!=$lambdacounter), a_counter=$a_counter\n";
		    exit(0);		
		}
	    }
	    $Nlambda = $lambdacounter;
	    $lambdacounter = 0;

	    if ($a_counter == 1){
		for (my $j=1;$j<=4;$j++){
		    $output[$j][$lambdacounter] = sprintf("#%11s %10.3f","lambda/loga",log10($a_micron));
		    $lambda[$lambdacounter] = 0.0;
		}	    
	    } else {
		for (my $j=1;$j<=4;$j++){
		    $output[$j][$lambdacounter] .= sprintf(" %10.3f",log10($a_micron));
		}	    	    
	    }
	    print STDERR "$a_counter: radius(microns) = $a_micron\n";

	    $lambdacounter++;
	    next;
	}	
	if ($a_counter<1){
	    print STDERR "ERROR: The first line does not contain radius(micron)!\n";
	    exit(0);		
	}
	my @data = split(/\s+/,$line);
	#print STDERR "TEST!!!!: <$line> $data[4]\n";
	if ($a_counter == 1){
	    for (my $j=1;$j<=4;$j++){
		$output[$j][$lambdacounter] = sprintf("%12.3e %10.3e",$data[0],$data[$j]);
		$lambda[$lambdacounter] = $data[0];
	    }	    
	} else {
	    for (my $j=1;$j<=4;$j++){
		$output[$j][$lambdacounter] .= sprintf(" %10.3e",$data[$j]);
		#print STDERR "FFF: $a_counter $lambdacounter $data[0] $output[$j][$lambdacounter]\n";
		if ($lambda[$lambdacounter] != $data[0]){
		    print STDERR "BUG: $lambda[$lambdacounter] != $data[0], $lambdacounter, <@data>!\n";
		    exit(0);				    
		}
	    }	    	    
	}
	$lambdacounter++;
    }
    for (my $j=1;$j<=4;$j++){
	open (DUMMYFILE, ">$filename.$colnames{$j}.2dtab");
	print DUMMYFILE "# x-dimension: dust grain radius 'a' in microns\n";
	print DUMMYFILE "# y-dimension: wavelength 'lambda' in microns\n";
	print DUMMYFILE "# values: $colnames{$j}(x,y)\n";
	my $dloga = &calclogspacing($output[$j][0]);
	print DUMMYFILE "# dloga: $dloga\n";	
	for (my $l=0;$l<$Nlambda;$l++){	   
	    print DUMMYFILE "$output[$j][$l]\n";
	}
	close (DUMMYFILE);
    }
}

sub calclogspacing{
    my ($logastring) =@_;    
    $logastring =~ s/^\#\s*lambda\/loga\s+//;
    #print STDERR "VVVV:$logastring\n";
    my @loga = split(/\s+/,$logastring);
    undef my @dloga;
    for(my $i=0;$i<@loga;$i++){
	if ($i>0){
	    $dloga[$i] = $loga[$i] - $loga[$i-1];
	} else {
	    $dloga[$i] = 0.0;
	}
	$dloga[$i] = int($dloga[$i]*100.0+0.5)/100.0;
	#print "$loga[$i] $dloga[$i]\n";
    }
    $dloga[0] = $dloga[1];
    @tmp = &CleanMultipleEntries(@dloga);
    if (scalar(@tmp)!=1){
	print STDERR "ERROR: More than one dloga! @tmp\n";
	exit(0);			
    }
    return($tmp[0]);
}

