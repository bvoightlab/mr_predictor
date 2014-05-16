#!/usr/bin/perl -w
#mr_predictor.pl
#
#a program to simulate multilocus genotypes of a intermediate phenotype (related to a terminal phenotype), to construct effect those loci have
# on the destination phenotype, given a known relationship between intermediate and destination phenotypes 
#
#BF Voight
#created: 01.09.09

#INFO ON PROGRAM 
sub VERSION() { "v0.027" } ; #Current version of program
sub MODIFIED() { "05.16.14" } ; #Date of current version

#outline for script

#1. Read in a 'locus' file which contains: 
#(1) SNP 
#(2) risk and non-risk allele
#(3) Additive Variance term (lnOR Effect sizes, or beta coefficients for additive term)
#(4) SE on Va (-9 if you assume a model without Va (Va = 0 always)  
#(5) Dominance Variance term (lnOR effect size, or beta coefficients for domdev term)
#(6) SE on Vd (-9 if you assume a model without Vd (Vd = 0 always)
#(7) frequency of risk allele

#2. Read in phenotype info file which contains:
#(1) phenotype label
#(2) (I)ntermediate or (D)estination phenotype flag
#if (I), then
#(3a) effect size of intermediate to destination (units of beta or ln(OR))
#(4a) SE on effect size
#if (D), then
#(3b) (L)iability or (Q)uantitative model
#(4b) overall prevalence of destination phenotype (liability model) (L) or -9

#3. Modelling: intermediate phenotype
#(a) calculate sigma_r, the variance contributed by the loci that you are simulating.
#sigma_r for a single locus should be:
#
# sig_r = 2*p*q*(a + d(q-p))^2 + (2pqd)^2
#
#(b) calculate sigma_t, which is sig_assumed - sigma_r. in most cases, I'm assuming a model of mu_assumed=0, sig_assumed=1, so this is just = 1-sigma_r
#--make sure you check that sigma_t is not negative (i.e., sigma_r > sigma_assumed).

#4. Modelling: destination phenotype
#(a) figure out if the destination is a quantiative or boolean trait.
#(b1_1) if boolean, calculate the prevalence of the given intermediate trait distribution assuming the relationship of intermediate to destination
#will have to use numerical integration (approximation hack, say).
#if Kp > Kassume, then assert Kp rather than Kass and warn, extra_k = 0;
#if Kp > 1 or Kaassume > 1, die
#if if Kp < Kassume, add the difference in terms to probability of extra_k = Kassume-Kp to liability model. 
#(b1_2) model prob(aff) ~ intermediate + other predictors + extra_k. 
#(b2_1) if quantiative trait, you can draw model specifically:
#
# Dpheno = N(0,sigma_r2) + sqrt(r2)*Ipheno
#
#where sigma_r2 = 1-r^2, and r^2 is the correlation between Dpheno and Ipheno [drawn from the distribution of correlation]

#5. Simulations

#to do list
# add check for exclude-mean-adj to make sure traits are present (and give warning if not)
# note that some characters in phenotype names when shifted from perl to R change -- e.g. "-" character changed to "." e.g. HDL-G to HDL.G
# could make machinery faster by not iterating over all SNP counts when tabulating ipheno or dpheno distortions from genetics
# bug handling when "-" characters in the scoresheet are long dash characters (excel template issue)
# check to make sure check all SNPs listed in scoresheet are present in the infofile (otherwise die).


#Things to look at/Options to add:
#2. continuous traits that are not normally distributed: e.g. "age"
#3. implement a "basic" hazard - time to event model
#5. multiple ascertainment conditions for multiple endpoints (e.g. condition on collecting T2D _AND_ MI, say)
#6. ascertainment condition on intermediate traits (e.g. patients with LDL > some value)
#
#7. family-based designs
#8. fix correlation/covariance calculation issue.
#10. add more verbose output to scorefile input (e.g. how many SNPs are read in per trait, how many SNPs overlap, intersect with endpoint(s), etc.)
#11. cleanup of some revised old function (get_ipheno_chunk is obsolete, some bits in mk_varcovfile)
#12. need more robust testing of the simulated intermediates given genetic simulations (with mean adjustments). Does this all work precisely right?

#revision history

#v0.026 
#05.16.14: 
# commented out reference to .forR for var/covar matrix generation (file name convention still maintained)
# added a static (but variable) default for the version and date of the program             

#v0.024, v0.025: no major changes; just setup of git repository

#v0.023:
#11.19.13: dropped the use of an R calling function and revert explicitly to Math::Random functions, now a required installed module
#### download  
#          Math::Random 0.71 http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.71.tar.gz
#
#### installation (without root access)
#          perl Makefile.PL INSTALL_BASE=~/bin/Perllibs
#          make
#          make test
#          make install
#
#### setup for .bash_profile           
#          set PERL5LIB=$PERL5LIB:$HOME/bin/Perllibs/lib/perl5/x86_64-linux-thread-multi
#          export PERL5LIB
#
#11.20.13: moved evaluations of phenofile type "d/i/c" to case insensitive
#12.01.13: added covariates with effects to prediction model (e.g. SMOKING, T2D, ) see sim_covars
#12.02.13: added an option to --cc_nsamp. If ncases -1, then simulate a static number of samples (regardless of number of cases or controls obtained). e.g. pseudo odds of '10-year hazard' informally.
#12.03.13: added an option to skip association testing in plink (--skip-plink)

#v0.022
#
#10.31.13: added a 'phenotype-verbose' argument: outputs each intermediate phenotype perturbation to file 
#10.31.13: modified R function to report 6 significant digits (previously was 5)
#          also modified calc_mean_adj to report 6 sig digits (previously was no limit)
#          also modified the sig digits from the step generating the overall phenotypic distortion
#10.31.13: moved score calculation to "print scores" condition (speedup)
#10.31.13: added the expected contribution due to genetics (EXPvG) "boost" to the score stat file
#10.31.13: Modified do_pheno_distort to save ipheno contribution that stem from SNPs which modulate an ipheno and an endpoint
#          modified sim_dpheno to address when genetic variant is specified for both an intermediate and an outcome (uses endpoint for disease by default).
#
#10.30.13: added an explicit call to srand() to ensure results are reproducible
#10.30.13: also added explicited seed calling to R; updated mr_predictor_sampler.R to mr_predictor_sampler_v2.R
#10.30.13: fixed a bug in get_iiphenofile where default variances were not being initialized.
#10.30.13: renovated output data files to explicitly to function calls. 
#10.30.13: added argument to explicitly calculate (and print) score files; default is not to print scores
#10.30.13: fixed a couple of bugs argument parsing (copy errors for hazard function)
#10.30.13: improved phenotypic reporting, and implemented check for SNP overlap and 'dump' of phenotypic (genetic/otherwise) perturbation
#10.30.13: moved definition of "ipheno_chunk" (number of random multivar norm samples) drawn when R is called; added parameter to 
#          govern how many samples are drawn at once (potential speedup, depending).

#v0.021
#10.03.13: added a --set-ipheno-var flag + error handling in iiphenofile to allow for specification of arbitrary variance in intermediate traits
#default behavior is DEFAULT_VAR constant fucntion, which is assumed to be 1
#
#this means that the ii-rel file is now a variance/correlation matrix (which is a variance/covariance matrix under default conditions). 
#I still need to fix the calculation when variance != 1 to convert correlations to covariances
#
#10.09.13: fixed a bug in the iiphenofile function when there is only one trait file and no entries in the iifile (default variance not specified)

#10.08.13: converted sim_dpheno function to a formal logit model for simulations (see Wu et al, AJHG 2011 89(1):82-93)

#v0.02
#
#tested variance/covariance calculations for additive models (assuming no dominance). Need to check for models which include a domdev term (and se)
#added score calculation and reporting framework
#06.27.09: fixed a bug with how I was advancing random draws of intermediate -> trait effects
#01.24.12: fixed a bug when outputting the var/cov matrix R file, if phenotypes are named similarly (e.g. fHDL and HDL).
#02.18.13: allow for ACGT instead of 1234 as code for alleles
#02.19.13: fixed a bug which led to an undef for the mean_pheno_adj when there are no snps underlying a trait
#02.19.13: added a new argument --mean-adj trait1,trait2,trait gives a lists of traits which a mean adjustment will be applied.

#v0.01
#
#initial code and framework (see above).

###for outputting purposes
$| = 1;

########
#MODULES INCLUDED

use Math::Random;

# @mean = (0,0);
# push @vcovm, [ ( 1, 0.8 ) ];
# push @vcovm, [ ( 0.8, 1 ) ];
#
# push @ref_arr, \@{$vcovm[0]};
# push @ref_arr, \@{$vcovm[1]};
#
# @results = Math::Random::random_multivariate_normal(1000,@mean, @vcovm);
#
# foreach $res (@results) {
#     print "$res->[0] $res->[1]\n";
# } 
#
# exit();

sub string_numerically { $a cmp $b; }

sub print_header {
    my ($filehandle) = @_;
    print $filehandle "\n";
    print $filehandle "#----------------------------------#\n";
    print $filehandle "# mr_predictor # " . VERSION() . " # " . MODIFIED() . " #\n";
    print $filehandle "#----------------------------------#\n";
    print $filehandle "#      (c) Benjamin F. Voight      #\n";
    print $filehandle "#----------------------------------#\n";
    print $filehandle "\n";
}

sub print_string {
    my ($string, $out_fh, $log_fh) = @_;
    print $out_fh "$string";
    print $log_fh "$string";
}

sub print_covars {
    my ($ref_sim_covars, $ref_nlabel, $out_fh) = @_;

    if ($$ref_nlabel == 1) { #print header
	print $out_fh "FID IID";
	foreach my $covar (sort string_numerically keys %$ref_sim_covars) {
	    print $out_fh " $covar";
	}
	print $out_fh "\n";
    }

    #print labels
    print $out_fh "$$ref_nlabel $$ref_nlabel";

    #print covars
    foreach my $covar (sort string_numerically keys %$ref_sim_covars) {
	print $out_fh " $$ref_sim_covars{$covar}";
    }
    print $out_fh "\n";
}

sub print_simpheno {
    my ($ref_ipheno_order, $ref_dpheno_order, $ref_dis_ipheno, $ref_sim_dpheno_data, $ref_nlabel, $out_fh) = @_;
    
    if ($$ref_nlabel == 1) { #print header
	print $out_fh "FID IID";
	for (my $p=0; $p<scalar(@{$ref_ipheno_order}); $p++) {
	    print $out_fh " $$ref_ipheno_order[$p]";
	}
	for (my $p=0; $p<scalar(@{$ref_dpheno_order}); $p++) {
	    print $out_fh " $$ref_dpheno_order[$p]";
	}
	print $out_fh "\n";
    }

    #print phenotypes
    print $out_fh "$$ref_nlabel $$ref_nlabel";
    
    #print intermediate traits first
    for (my $p=0; $p<scalar(@{$ref_dis_ipheno}); $p++) {
	print $out_fh " $$ref_dis_ipheno[$p]";
    }

    #do destination next
    for (my $p=0; $p<scalar(@{$ref_sim_dpheno_data}); $p++) {
	print $out_fh " $$ref_sim_dpheno_data[$p]";
    }
    print $out_fh "\n";
}

sub print_simverbpheno {
    my ($ref_ipheno_order, $ref_dpheno_order, $ref_dis_ipheno, $ref_sim_dpheno_data, $ref_verbpheno, $ref_nlabel, $out_fh) = @_;

    if ($$ref_nlabel == 1) { #print header
	print $out_fh "FID IID";
	for (my $p=0; $p<scalar(@{$ref_ipheno_order}); $p++) {
	    print $out_fh " $$ref_ipheno_order[$p]";
	    
	    #print out the env and genetic contributions separately (order is env, genetic)
	    print $out_fh " vE_@{[$$ref_ipheno_order[$p]]} vG_@{[$$ref_ipheno_order[$p]]}";
	}
	
	for (my $p=0; $p<scalar(@{$ref_dpheno_order}); $p++) {
	    print $out_fh " $$ref_dpheno_order[$p]";
	}
	print $out_fh "\n";
    }

    #print phenotypes
    print $out_fh "$$ref_nlabel $$ref_nlabel";

    #print intermediate traits first
    for (my $p=0; $p<scalar(@{$ref_dis_ipheno}); $p++) {
	print $out_fh " $$ref_dis_ipheno[$p]";
	print $out_fh " $$ref_verbpheno{$$ref_ipheno_order[$p]}[0] $$ref_verbpheno{$$ref_ipheno_order[$p]}[1]";
    }

    #do destination next
    for (my $p=0; $p<scalar(@{$ref_sim_dpheno_data}); $p++) {
	print $out_fh " $$ref_sim_dpheno_data[$p]";
    }
    print $out_fh "\n";  
}

sub print_simscore {
    my ($ref_ipheno_order, $ref_dpheno_order, $ref_scorestats, $asc_pheno_lab, $asc_pheno_val, $ref_nlabel, $out_fh) = @_;
    
    #note that scores are based on the score data file provided (i.e., expected value), not what is drawn for a given simulation.

    if ($$ref_nlabel == 1) { #print header
	print $out_fh "FID IID $asc_pheno_lab";
	for (my $p=0; $p<scalar(@{$ref_ipheno_order}); $p++) {
	    print $out_fh " SCORE_@{[$$ref_ipheno_order[$p]]}";
	    print $out_fh " EXPvG_@{[$$ref_ipheno_order[$p]]}"; #this is the expected addition to the trait due to genetics
	}
	for (my $p=0; $p<scalar(@{$ref_dpheno_order}); $p++) {
	    print $out_fh " SCORE_@{[$$ref_dpheno_order[$p]]}";
	}
	print $out_fh "\n";
    }

    #print scores
    print $out_fh "$$ref_nlabel $$ref_nlabel $asc_pheno_val"; #print the ascertainment phenotype status for ease

    #do intermediate traits first
    for (my $p=0; $p<scalar(@{$ref_ipheno_order}); $p++) {
	if ($$ref_scorestats{$$ref_ipheno_order[$p]}[1] == 0) {
	    print $out_fh " -9 -9";
	} else {
	    #this is the score normalized for nloci;
	    printf $out_fh " %1.5f", $$ref_scorestats{$$ref_ipheno_order[$p]}[0]/$$ref_scorestats{$$ref_ipheno_order[$p]}[1];
	    printf $out_fh " %1.5f", $$ref_scorestats{$$ref_ipheno_order[$p]}[0];
	}
    }

    #do destination next
    for (my $p=0; $p<scalar(@{$ref_dpheno_order}); $p++) {
	if ($$ref_scorestats{$$ref_dpheno_order[$p]}[1] == 0) {
	    print $out_fh " -9";
	} else {
	    printf $out_fh " %1.5f", $$ref_scorestats{$$ref_dpheno_order[$p]}[0]/$$ref_scorestats{$$ref_dpheno_order[$p]}[1];
	}
    }

    print $out_fh "\n";
}

sub print_simped {
    my ($ref_genodata, $ref_sim_covars, $ref_this_pheno, $ref_nlabel, $out_fh) = @_;
    my $snp;
   
    print $out_fh "$$ref_nlabel $$ref_nlabel 0 0 $$ref_sim_covars{'SEX'} $$ref_this_pheno ";
    foreach $snp (sort string_numerically keys %$ref_genodata) {
	print $out_fh " $$ref_genodata{$snp}[0] $$ref_genodata{$snp}[1]";
    }

    print $out_fh "\n";
}

sub check_file_exists {
    my ($file, $out_fh, $log_fh) = @_;
    my $myprint;

    stat($file);
    if ( !(-e _) ) {
        my $myprint = "ERROR: Can't locate " . $file . ".\n";
        print_string($myprint, $out_fh, $log_fh);
        exit();
    }
}

sub gen_Rseed {
    return(int(rand(1000000000))); #random seed feed to R are 1 to a billion
}

sub normal { #This distribution checks out
    my ($mu, $sigma) = @_;
    if ($sigma < 0) {
        print "Improper Variance.\n";
        exit();
    }
    my ($p1, $p2, $p);
    do {
        $p1 = -1 + 2*rand();
        $p2 = -1 + 2*rand();
        $p = $p1 * $p1 + $p2 * $p2;
    } while ( $p >= 1. );
    return $mu + $sigma * $p1 * sqrt(-2 * log ($p) / $p );
}

sub get_phenofile {
    my ($file, $ref_phenodata, $ref_covars_exist, $out_fh, $log_fh) = @_;
    my ($readline, $myprint);
    my @entry;
    my $ntot_i = 0;
    my $ntot_d = 0;
    my $ntot_c = 0;

    open PHENO, "<$file" or die "Can't open $file!\n";
    $myprint = "Reading in phenotype info from [ " . $file . " ]\n";
    print_string($myprint, $out_fh, $log_fh);
    while ($readline = <PHENO>) {
        @entry = split '\s+', $readline;
        if (scalar(@entry) != 4) {
            $myprint = "Warning: Different number of entries found in [ " . $file . " ]. Expected 4, but found " . scalar(@entry) . "\n";
            print_string($myprint, $out_fh, $log_fh);
            $myprint = "Last line read: " . $readline . "\n";
            print_string($myprint, $out_fh, $log_fh);
            exit();
        }

	if (!defined($$ref_phenodata{$entry[0]}[0])) {
	    @{$$ref_phenodata{$entry[0]}} = ($entry[1], $entry[2], $entry[3]); #key based on name; store type, trait model and Kp
	    if ($entry[1] =~ m/i/i) {
		$ntot_i++;
	    } elsif ($entry[1] =~ m/d/i) {
		$ntot_d++;

		if ($entry[2] !~ m/b/) {
		    $myprint = "Warning: improper entry found for type of phenotype for $entry[0] [ $entry[2] ]. Requires 'b'.\n";
		    print_string($myprint, $out_fh, $log_fh);
		    exit();
		}

	    } elsif ($entry[1] =~ m/c/i) {
		$ntot_c++;
		$$ref_covars_exist = 1; #you found covariates, make sure these data are reported.

		if ($entry[2] !~ m/b/) {
		    $myprint = "Warning: improper entry found for type of phenotype for $entry[0] [ $entry[2] ]. Requires 'b' ['q' not implemented].\n";
		    print_string($myprint, $out_fh, $log_fh);
		    exit();
		}

	    } else {
		$myprint = "Warning: improper entry found for type of phenotype for $entry[0] [ $entry[1] ]. Requires 'i', 'c', or 'd'.\n";
		print_string($myprint, $out_fh, $log_fh);
		exit();
	    }

	} else {
	    $myprint = "Warning: $entry[0] found twice in phenotype file! Exiting.\n";
            print_string($myprint, $out_fh, $log_fh);
            exit();
	}
    }
    close(PHENO);

    $myprint = $ntot_i . " intermediate phenotypes read from [ " . $file . " ]\n";
    print_string($myprint, $out_fh, $log_fh);
    $myprint = $ntot_c . " covariates read from [ " . $file . " ]\n";
    print_string($myprint, $out_fh, $log_fh);
    $myprint = $ntot_d . " destination phenotypes read from [ " . $file . " ]\n";
    print_string($myprint, $out_fh, $log_fh);

    #check if prev of sex is defined. If not, set to default
    if (!defined($$ref_phenodata{'SEX'}[0])) {
	@{$$ref_phenodata{'SEX'}} = ('c', 'b', SEX_PREV()); 
	$myprint = "Prevalence of SEX not defined. Using default ratio of men/women of [ " . SEX_PREV() . " ]\n";
	print_string($myprint, $out_fh, $log_fh);
    }

}

sub get_iiphenofile {
    my ($file, $ref_iiphenodata, $ref_phenolist, $outfix, $flag_specvar, $out_fh, $log_fh) = @_;
    my ($readline, $myprint, $n, $exp_tally, $unspecfile);
    my @entry;
    my $ntot = 0;

    #count the number of intermediate phenotypes.
    foreach my $pheno (keys %$ref_phenolist) {
	if ($$ref_phenolist{$pheno}[0] =~ m/i/i) { #this is an intermediate phenotype
	    $n++;
	}
    }
    $exp_tally = ($n * ($n-1))/2; 
    
    #determine if variances are specified in the file or not
    if ($flag_specvar == 0) {
	$myprint = "Assuming intermediate traits all have Variance = @{[ DEFAULT_VAR() ]}\n";
	print_string($myprint, $out_fh, $log_fh);
	foreach my $i (keys %$ref_phenolist) {
	    $$ref_iiphenodata{$i}{$i} = DEFAULT_VAR();
	}
    } else {
	$myprint = "Obtaining variances for intermediate traits from intermediate trait file.\n";
	print_string($myprint, $out_fh, $log_fh);
	$exp_tally += $n;
    }    

    if ($exp_tally == 0) { #only one trait specified
	$myprint = "One intermediate trait specified (default variance assumed), no phenotype relationships expected/read in from [ " . $file . " ]\n";
	print_string($myprint, $out_fh, $log_fh);
	foreach my $i (keys %$ref_phenolist) {
	    $$ref_iiphenodata{$i}{$i} = DEFAULT_VAR();
	}
    } else { #more than one trait specified; expecting covariance terms (and possibly variance terms)
	open PHENO, "<$file" or die "Can't open $file!\n";
	#$myprint = "Reading in intermediate->intermediate phenotype relationships from [ " . $file . " ]\n";
	#print_string($myprint, $out_fh, $log_fh);
	while ($readline = <PHENO>) {
	    @entry = split '\s+', $readline;
	    if (!defined($$ref_phenolist{$entry[0]}[0])) { 
		$myprint = "ERROR: Found unknown phenotype [ $entry[0] ] in relationships file!\n";
		print_string($myprint, $out_fh, $log_fh);
		$myprint = "Last line read: " . $readline . "\n";
		print_string($myprint, $out_fh, $log_fh);
		exit();
	    } elsif (!defined($$ref_phenolist{$entry[1]}[0])) {
		$myprint = "ERROR: Found unknown phenotype [ $entry[1] ] in relationships file!\n";
		print_string($myprint, $out_fh, $log_fh);
		$myprint = "Last line read: " . $readline . "\n";
		print_string($myprint, $out_fh, $log_fh);
		exit();
	    }
	    $$ref_iiphenodata{$entry[0]}{$entry[1]} = $entry[2]; 
	    $$ref_iiphenodata{$entry[1]}{$entry[0]} = $entry[2]; #symmetric
	    $ntot++;
	}
	
	$myprint = "Expecting $exp_tally intermediate phenotype relationship(s) from [ " . $file . " ]\n";    
	print_string($myprint, $out_fh, $log_fh);        
	$myprint = $ntot . " intermediate->intermediate phenotype relationships read from [ " . $file . " ]\n";    
	print_string($myprint, $out_fh, $log_fh);
	
	if ($exp_tally > $ntot) {
	    $unspecfile = $outfix . "_unspec.err";
	    $myprint = "Logging all unspecified relationships that were set to 0 (correlation) or @{[ DEFAULT_VAR() ]} (variance) to [ $unspecfile ]\n";
	    print_string($myprint, $out_fh, $log_fh);
	    open UNSPEC, ">$unspecfile" or die "Can't open $unspecfile!\n";
	    foreach my $i (keys %$ref_phenolist) {
		foreach my $j (keys %$ref_phenolist) {
		    if ($$ref_phenolist{$i}[0] =~ m/i/i && $$ref_phenolist{$j}[0] =~ m/i/i) { #insist on these all being intermediates
			if (!defined($$ref_iiphenodata{$i}{$j}) && !defined($$ref_iiphenodata{$j}{$i})) {
			    print UNSPEC "$i $j\n";
			    if ($i !~ m/$j/) {
				$$ref_iiphenodata{$i}{$j} = 0;
				$$ref_iiphenodata{$j}{$i} = 0;
			    } else {
				$$ref_iiphenodata{$i}{$j} = DEFAULT_VAR(); 
			    }
			} elsif (!defined($$ref_iiphenodata{$i}{$j}) && defined($$ref_iiphenodata{$j}{$i})) {
			    if ($i !~ m/$j/) {
				$$ref_iiphenodata{$i}{$j} = 0;
			    } else {
				$$ref_iiphenodata{$i}{$j} = DEFAULT_VAR(); 
			    }
			} elsif (defined($$ref_iiphenodata{$i}{$j}) && !defined($$ref_iiphenodata{$j}{$i})) {
			    if ($i !~ m/$j/) {
				$$ref_iiphenodata{$j}{$i} = 0;
			    } else {
				$$ref_iiphenodata{$i}{$j} = DEFAULT_VAR(); 
			    }
			}
		    }
		}
	    }
	    close(UNSPEC);
	} elsif ($exp_tally < $ntot) {
	    $myprint = "ERROR: More intermediate relationships than expected.\n";
	    print_string($myprint, $out_fh, $log_fh);
	    exit();
	}

    } #END number of results in ii file > 1

    return($ntot);
}

sub get_idphenofile {
    my ($file, $ref_idphenodata, $ref_phenolist, $ref_asc_pheno, $out_fh, $log_fh) = @_;
    my ($readline, $myprint, $gotit);
    my @entry;
    my $ntot = 0;
    my %lookup_check;
    
    open PHENO, "<$file" or die "Can't open $file!\n";
    #$myprint = "Reading in intermediate->destination phenotype relationships from [ " . $file . " ]\n";
    #print_string($myprint, $out_fh, $log_fh);
    while ($readline = <PHENO>) {
        @entry = split '\s+', $readline;
        if (!defined($$ref_phenolist{$entry[0]}[0])) {
            $myprint = "ERROR: Found unknown phenotype [ $entry[0] ] in relationships file!\n";
            print_string($myprint, $out_fh, $log_fh);
            $myprint = "Last line read: " . $readline . "\n";
            print_string($myprint, $out_fh, $log_fh);
            exit();
        } elsif (!defined($$ref_phenolist{$entry[1]}[0])) {
            $myprint = "ERROR: Found unknown phenotype [ $entry[1] ] in relationships file!\n";
            print_string($myprint, $out_fh, $log_fh);
            $myprint = "Last line read: " . $readline . "\n";
            print_string($myprint, $out_fh, $log_fh);
            exit();
        }
        @{$$ref_idphenodata{$entry[0]}{$entry[1]}} = ($entry[2], $entry[3]); #keyed on D, I; store beta and SE for relationship.
	$lookup_check{$entry[1]} = 1; #TRUE if intermediate hooks up with a destination.
        $ntot++;
    }

    $myprint = $ntot . " intermediate->destination phenotype relationships read from [ " . $file . " ]\n";
    print_string($myprint, $out_fh, $log_fh);

    #Check for specified effect of SEX. if does not exist, define as no effect
    if (!defined($$ref_idphenodata{$$ref_asc_pheno}{'SEX'}[0])) {
	@{$$ref_idphenodata{$entry[0]}{'SEX'}} = (SEX_EFFECT(), 0);
	$myprint = "No effect of sex specified. Assuming Beta = [ " . SEX_EFFECT() . " ]\n";
	print_string($myprint, $out_fh, $log_fh);
    }
    
    #check if all intermediate phenotypes relate to a destination phenotype. Warn if failure (but don't kill).
    for my $i (keys %$ref_phenolist) {
	if ($$ref_phenolist{$i}[0] =~ m/i/i) {
	    if (!defined($lookup_check{$i})) {
		$gotit = 0;
	    } else {
		$gotit = 1;
	    }
	} else { #otherwise you are a destination trait, and thus fine.
	    $gotit = 1; 
	}

	if ($gotit == 0) {
	    $myprint = "Warning: intermediate phenotype $i does not relate to a destination phenotype!\n";
	    print_string($myprint, $out_fh, $log_fh);
	}
    }
}

sub get_infofile {
    my ($infofile, $ref_infodata, $out_fh, $log_fh) = @_;
    my ($readline, $myprint);
    my @entry;
    my $ntot = 0;
    
    open INFO, "<$infofile" or die "Can't open $infofile!\n";
    #$myprint = "Reading in SNP base information from [ " . $infofile . " ]\n";
    #print_string($myprint, $out_fh, $log_fh);
    while ($readline = <INFO>) {
        @entry = split '\s+', $readline;
        if (scalar(@entry) != 4) {
            $myprint = "Warning: Different number of entries found in [ " . $infofile . " ]. Expected 4, but found " . scalar(@entry) . "\n";
            print_string($myprint, $out_fh, $log_fh);
            $myprint = "Last line read: " . $readline . "\n";
            print_string($myprint, $out_fh, $log_fh);
            exit();
        }
        #keyed by SNP label, then risk, nrisk, freq
	if (!defined($$ref_infodata{$entry[0]}[0])) {
	    @{$$ref_infodata{$entry[0]}} = ($entry[1], $entry[2], $entry[3]);
	    $ntot++;
	} else {
	    $myprint = "Warning: $entry[0] found twice in info file! Exiting.\n";
	    print_string($myprint, $out_fh, $log_fh);
	    exit();
	}
    }
    close(INFO);

    $myprint = $ntot . " SNPs read from [ " . $infofile . " ]\n";
    print_string($myprint, $out_fh, $log_fh);
    return($ntot);
}

sub get_scorefile {
    my ($scorefile, $ref_scoredata, $out_fh, $log_fh) = @_;
    my ($readline, $myprint);
    my @entry;
    my $ntot = 0;

    open SCORE, "<$scorefile" or die "Can't open $scorefile!\n";
    #$myprint = "Reading in SNP score info from [ " . $scorefile . " ]\n";
    #print_string($myprint, $out_fh, $log_fh);
    while ($readline = <SCORE>) {
        @entry = split '\s+', $readline;
	if ($entry[0] =~ m/\#/) {
	    #skip comments.
	} else {
	    if (scalar(@entry) != 6) {
		$myprint = "Warning: Different number of entries found in [ " . $scorefile . " ]. Expected 6, but found " . scalar(@entry) . "\n";
		print_string($myprint, $out_fh, $log_fh);
		$myprint = "Last line read: " . $readline . "\n";
		print_string($myprint, $out_fh, $log_fh);
		exit();
	    }
	    
	    #keyed by SNP rs number and phenotype; then add_fx, add_fx_se, dom_fx, dom_fx_se
	    if (!defined($$ref_scoredata{$entry[0]}{$entry[1]})) {
		@{$$ref_scoredata{$entry[0]}{$entry[1]}} = ($entry[2], $entry[3], $entry[4], $entry[5]);
		$ntot++;
	    } else {
		$myprint = "Warning: $entry[0] found twice for the same phenotype in score file! Exiting.\n";
		print_string($myprint, $out_fh, $log_fh);
		exit();
	    }
	}
    }
    close(SCORE);

    $myprint = $ntot . " SNP->Phenotype relationships read from [ " . $scorefile . " ]\n";
    print_string($myprint, $out_fh, $log_fh);
    return($ntot);    
}

sub get_ipheno_chunk{ 
    my ($ipheno_chunk, $ref_ipheno_order, $ref_sim_ipheno_data, $out_fh, $log_fh) = @_;
    my ($readline, $myprint);
    my @entry;

    open CHUNK, "<$ipheno_chunk" or die "Can't open $ipheno_chunk!\n";
    #$myprint = "Reading in simulated intermediate trait values from [ " . $ipheno_chunk . " ]\n";
    #print_string($myprint, $out_fh, $log_fh);

    #get the header, reserve the phenotype order
    $readline = <CHUNK>; #header;
    chomp($readline);
    @entry = split '\s+', $readline;
    @{$ref_ipheno_order} = @entry;

    while ($readline = <CHUNK>) {
	chomp($readline);
	@entry = split '\s+', $readline;
	push @{$ref_sim_ipheno_data}, [ @entry ];
    }
    close(CHUNK);    
}

sub mk_mapfile {
    my ($ref_scoredata, $outfix, $out_fh, $log_fh) = @_;
    my ($snp, $a, $b);
    my $pos = 1000;

    open MAP, ">@{[$outfix]}.map" or die "Can't open @{[$outfix]}.map!\n";
    $myprint = "Outputting generic map file to [ " . $outfix . ".map ]\n";
    print_string($myprint, $out_fh, $log_fh);
    
    foreach $snp (sort string_numerically keys %$ref_scoredata) {
	print MAP "1 $snp 0 $pos\n";
	$pos += 1000;
    } 
    close(MAP);
}

sub mk_basedata {
    my ($ref_scoredata, $ref_basedata, $ref_infodata, $ref_phenolist, $rand_flag, $out_fh, $log_fh) = @_;
    my $a_fx = 0;
    my $d_fx = 0;
    my %phenolist;
    my %snplist;
    
    foreach my $snp (keys %$ref_infodata) {
	foreach my $pheno (keys %{$$ref_scoredata{$snp}}) {
	    if (!defined($$ref_scoredata{$snp}{$pheno}[0])) { #no score value for this SNP for this phenotype. set to no effect.
		$a_fx = 0;
		$d_fx = 0;
	    } else {
	    	#draw an additive term if applicable.
		if ($$ref_scoredata{$snp}{$pheno}[1] != -9) { 
		    if ($rand_flag == 1) { #random draw from fx distribution 
			$a_fx = normal($$ref_scoredata{$snp}{$pheno}[0], $$ref_scoredata{$snp}{$pheno}[1]);
		    } else { #do not randomize
			$a_fx = $$ref_scoredata{$snp}{$pheno}[0];
		    }
		} else {
		    $a_fx = 0;
		}
	    
		if ($$ref_scoredata{$snp}{$pheno}[3] != -9) {
		    if ($rand_flag == 1) {
			$d_fx = normal($$ref_scoredata{$snp}{$pheno}[2], $$ref_scoredata{$snp}{$pheno}[3]);
		    } else {
			$d_fx = $$ref_scoredata{$snp}{$pheno}[2];
		    }
		} else {
		    $d_fx = 0;
		}
	    }
	    @{$$ref_basedata{$snp}{$pheno}} = ($a_fx, $d_fx, $$ref_infodata{$snp}[2]);
	}
	#record the SNP.
	$snplist{$snp} = 1;   
    }

    #go through the list and zero out all SNPs and pheno   
    foreach $snp (keys %snplist) {
	foreach $pheno (keys %$ref_phenolist) {
	    if (!defined($$ref_basedata{$snp}{$pheno}[0])) {
		@{$$ref_basedata{$snp}{$pheno}} = (0, 0, $$ref_infodata{$snp}[2]);
	    }
	}
    }

}

sub mk_varcovfile {
    my ($forRfile, $ref_vcovm, $ref_iiphenodata, $ref_phenolist, $ref_basedata, $ref_cov_adj_matrix, $out_fh, $log_fh) = @_;    
    my $print_endline;
    my @row;
    
    #open RVARCOV, ">$forRfile" or die "Can't open $forRfile!\n";
    #$myprint = "Outputting variance/covariance matrix to [ " . $forRfile . " ]\n";
    #print_string($myprint, $out_fh, $log_fh);
    
    #print RVARCOV "X";
    #foreach my $i (sort string_numerically keys %$ref_phenolist) {
    #	if ($$ref_phenolist{$i}[0] =~ m/i/i) {
    #	    print RVARCOV " $i";
    #	}
    #}
    #print RVARCOV "\n";

    foreach my $i (sort string_numerically keys %$ref_phenolist) {
	@row = ();
	if ($$ref_phenolist{$i}[0] =~ m/i/i) {
	    #print RVARCOV "$i";
	    foreach my $j (sort string_numerically keys %$ref_phenolist) {
		if ($$ref_phenolist{$j}[0] =~ m/i/i) {
		    if ($i =~ m/^$j$/) { #get the adjusted variance for this phenotype
			$val = sprintf("%1.4f", $$ref_iiphenodata{$i}{$j}-get_sigma_r($ref_basedata, $i, $out_fh, $log_fh));
			#print RVARCOV " $val";
			push @row, $val;
		    } else {
			#just a test
			#printf RVARCOV " %1.4f", ($$ref_iiphenodata{$i}{$j});
			$val = sprintf("%1.4f", $$ref_iiphenodata{$i}{$j} - $$ref_cov_adj_matrix{$i}{$j});
			#print RVARCOV " $val";
			push @row, $val;
		    }
		}
	    }
	    #print RVARCOV "\n";
	    push @$ref_vcovm, [ @row ] ;
	}
       
    }
    
    #close(RVARCOV);
}

sub get_sigma_r {
    my ($ref_basedata, $this_pheno, $out_fh, $log_fh) = @_;
    my ($a, $d, $p, $q) = 0;
    my $sig_r = 0;

    #3. Modelling: intermediate phenotype
    #(a) calculate sigma_r, the variance contributed by the loci that you are simulating.
    #sigma_r for a single locus should be:
    #
    # sig_r = 2*p*q*(a + d(q-p))^2 + (2pqd)^2

    #assume base data has the structure
    #add_fx, dom_fx, freq

    foreach my $snp (keys %$ref_basedata) {	
	$a = $$ref_basedata{$snp}{$this_pheno}[0];
	$d = $$ref_basedata{$snp}{$this_pheno}[1];
	$p = $$ref_basedata{$snp}{$this_pheno}[2];
	$q = 1-$p;
	$sig_r += 2*$p*$q*($a + $d*($q-$p))**2 + (2*$p*$q*$d)**2;
    }
    
    return($sig_r);
}

sub get_cov_adj_matrix {
    my ($ref_basedata, $ref_cov_adj_matrix, $ref_phenodata, $out_fh, $log_fh) = @_;
    my ($a1, $d1, $a2, $d2, $p, $q);
    
    #Modelling: intermediate phenotype covariance 
    #discussion with Ben Neale brings us to
    #
    #Cov(X1,X2) = (2*p*q)*(a1 + d1(q-p))*(a2 + d2(q-p)) + (2*p*q)^2*d1 + (2*p*q)^2*d2
    #
    #@{$$ref_basedata{$snp}{$pheno}} = ($a_fx, $d_fx, $$ref_infodata{$snp}[2]);
    #
    #this works assuming that the total variance = 1, and correlation == covariance. 

    #initialize the matrix.
    foreach my $pheno_1 (keys %$ref_phenodata) {
	if ($$ref_phenodata{$pheno_1}[0] =~ m/i/i) {
	    foreach my $pheno_2 (keys %$ref_phenodata) {
		if ($$ref_phenodata{$pheno_2}[0] =~ m/i/i) {
		    if ($pheno_1 !~ m/$pheno_2/) {
			$$ref_cov_adj_matrix{$pheno_1}{$pheno_2} = 0;
		    } 
		}
	    }
	}
    }

    foreach my $snp (keys %$ref_basedata) {
	foreach my $pheno_1 (keys %$ref_cov_adj_matrix) {
	    $a1 = $$ref_basedata{$snp}{$pheno_1}[0];
	    $d1 = $$ref_basedata{$snp}{$pheno_1}[1];
	    $p = $$ref_basedata{$snp}{$pheno_1}[2];
	    $q = 1-$p;	    
	    foreach my $pheno_2 (keys %{$$ref_cov_adj_matrix{$pheno_1}}) {		
		$a2 = $$ref_basedata{$snp}{$pheno_2}[0];
		$d2 = $$ref_basedata{$snp}{$pheno_2}[1];
		if (!defined($$ref_cov_adj_matrix{$pheno_1}{$pheno_2})) {
		    $myprint = "ERROR: {$pheno_1}{$pheno_2} combination is not defined. Exiting!\n";
		    print_string($myprint, $out_fh, $log_fh);
		    exit();
		} else {
		    $$ref_cov_adj_matrix{$pheno_1}{$pheno_2} += (2*$p*$q) * ($a1 + $d1*($q-$p)) * ($a2 + $d2*($q-$p)) + $d1*(2*$p*$q)**2 + $d2*(2*$p*$q)**2;
		}
	    }
	}
    }
    
}

sub calc_mean_adj {
    my ($ref_basedata, $ref_mean_pheno_adj, $ref_phenolist, $ref_excl_adj_traitlist, $ref_sim_meanvec, $ref_sim_ipheno_order) = @_;
    my ($hz_inc, $het, $hz_dec);

    #02.19.13 - if no SNPs underlying the trait, set to zero 
    #check to make sure each pheno gets an adjustment if necessary
    #initialize all i phenotypes to zero:     
    #key based on name; store type, trait model and Kp
    foreach my $pheno (keys %$ref_phenolist) {
	if ($$ref_phenolist{$pheno}[0] =~ m/i/i) { #intermediate trait 
	    $$ref_mean_pheno_adj{$pheno} = 0;
	}
    }
    
    #@{$$ref_basedata{$snp}{$pheno}} = (0, 0, $$ref_infodata{$snp}[2]);
    foreach my $snp (keys %$ref_basedata) {
	foreach my $pheno (keys %$ref_mean_pheno_adj) {
	    $hz_inc = $$ref_basedata{$snp}{$pheno}[2]**2 * $$ref_basedata{$snp}{$pheno}[0];
	    $het = 2 * $$ref_basedata{$snp}{$pheno}[2] * (1-$$ref_basedata{$snp}{$pheno}[2]) * $$ref_basedata{$snp}{$pheno}[1];
	    $hz_dec = -1 * (1-$$ref_basedata{$snp}{$pheno}[2])**2 * $$ref_basedata{$snp}{$pheno}[0];
	    if (!defined ($$ref_mean_pheno_adj{$pheno})) {
		$$ref_mean_pheno_adj{$pheno} = sprintf("%1.6f", $hz_inc + $het + $hz_dec);
	    } else {
		$$ref_mean_pheno_adj{$pheno} += sprintf("%1.6f", $hz_inc + $het + $hz_dec);
	    }
	}
    }

    #02.19.13: exclude some traits from phenotypic adjustment 
    foreach my $pheno (keys %$ref_phenolist) {
	#print "B-- $pheno -- $$ref_phenolist{$pheno}[0] -- $$ref_mean_pheno_adj{$pheno}\n";
	for (my $i=0; $i<scalar(@{$ref_excl_adj_traitlist}); $i++) {
	    if ($pheno =~ $$ref_excl_adj_traitlist[$i]) {
		$$ref_mean_pheno_adj{$pheno} = 0;
	    }
	}
	#print "A-- $pheno -- $$ref_phenolist{$pheno}[0] -- $$ref_mean_pheno_adj{$pheno}\n";
    }

    #11.19.13: initialize mean adj matrix
    foreach my $pheno (sort string_numerically keys %$ref_mean_pheno_adj) { 
        #it is -1 * this because the simulated data is added to the adjustment, thus bringing the total to a mu=0 trait 
	push @$ref_sim_meanvec, -1 * $$ref_mean_pheno_adj{$pheno};
	push @$ref_sim_ipheno_order, $pheno;
    }

}

sub calc_scorestat {
    my ($ref_simgeno, $ref_scorestat, $ref_infodata, $ref_basedata, $ref_scoredata) = @_;
    my $a1_cnt;
    
    foreach my $snp (keys %$ref_infodata) {
	foreach my $pheno (keys %{$$ref_basedata{$snp}}) {	    
	    $a1_cnt = 0;
	    
	    #get the count of A1 alleles
	    if ($$ref_simgeno{$snp}[0] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    }
	    if ($$ref_simgeno{$snp}[1] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    }

	    if (!defined($$ref_scorestat{$pheno})) {
		@{$$ref_scorestat{$pheno}} = (0, 0); #record score, and number of loci.
	    }

	    #print "$snp $pheno $$ref_scoredata{$snp}{$pheno}[0]\n";
	    #sleep 1;

	    #calculate a score based on the score data file provided (expected value), not what is drawn.
	    ##keyed by SNP rs number and phenotype; then add_fx, add_fx_se, dom_fx, dom_fx_se
	    if (!defined($$ref_scoredata{$snp}{$pheno}[0])) {
		#skip if the SNPxPHENO relationship does not exist in the score data file
	    } else {
		if ($a1_cnt == 0) { #2 non A1 alleles, decrease phenotype value by a
		    $$ref_scorestat{$pheno}[0] += -1*$$ref_scoredata{$snp}{$pheno}[0];
		} elsif ($a1_cnt == 1) { #add dom dev for heterozygote
		    $$ref_scorestat{$pheno}[0] += $$ref_scoredata{$snp}{$pheno}[2]; 
		} elsif ($a1_cnt == 2) { 
		    $$ref_scorestat{$pheno}[0] += $$ref_scoredata{$snp}{$pheno}[0];
		}

		$$ref_scorestat{$pheno}[1]++; #add one to the count.
	    }
	}
    }

}

sub sim_genotypes {
    my ($ref_infodata, $ref_simgeno) = @_;
    my ($allele_a, $allele_b);

    ##keyed by SNP label, then risk, nrisk, freq

    foreach my $snp (sort string_numerically keys %$ref_infodata) {
	if (rand() < $$ref_infodata{$snp}[2]) {
	    $allele_a = $$ref_infodata{$snp}[0];
	} else {
	    $allele_a = $$ref_infodata{$snp}[1];
	}

	if (rand() < $$ref_infodata{$snp}[2]) {
	    $allele_b = $$ref_infodata{$snp}[0];
	} else {
	    $allele_b = $$ref_infodata{$snp}[1];
	}
	@{$$ref_simgeno{$snp}} = ($allele_a, $allele_b);
    }
}

sub sim_covars {
    my ($ref_pheno, $ref_sim_covars) = @_;

    foreach my $covar (sort string_numerically keys %$ref_pheno) { #%pheno{'label'} = (trait class, trait type, prevalence) 
	#print "$covar $$ref_pheno{$covar}[0] $$ref_pheno{$covar}[1] $$ref_pheno{$covar}[2]\n"; 

	if ($$ref_pheno{$covar}[0] =~ m/c/i) {
	    if ($$ref_pheno{$covar}[1] =~ m/b/i) {
		if (rand() < $$ref_pheno{$covar}[2]) {
		    $$ref_sim_covars{$covar} = 1;
		} else {
		    $$ref_sim_covars{$covar} = 0;
		}
	    } elsif ($$ref_pheno{$covar}[1] =~ m/q/i) { 
                #not implemented
		$$ref_sim_covars{$covar} = -9; 
	    }
	}
    }

    #adjust the covar for SEX to a 1/2 (male/female) than 0/1
    $$ref_sim_covars{'SEX'}++;

}

sub do_dpheno_fx {
    my ($ref_simgeno, $ref_basedata, $ref_infodata, $phenolabel, $out_fh, $log_fh) = @_;
    my $fx_dpheno = 0;
    my $a1_cnt;

    #@{$$ref_basedata{$snp}{$pheno}} = ($a_fx, $d_fx, $$ref_infodata{$snp}[2]);
    #@{$$ref_simgeno{$snp}} = ($allele_a, $allele_b);   

    foreach my $snp (sort string_numerically keys %$ref_infodata) {
	$a1_cnt = 0;
	if (!defined($$ref_basedata{$snp}{$phenolabel}[0])) {
	    #for this phenotype, the SNP has no modelled effect
	    #do nothing
	} else {
	    #get the count of A1 alleles
	    if ($$ref_simgeno{$snp}[0] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    } 
	    if ($$ref_simgeno{$snp}[1] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    }

	    #print "$snp $a1_cnt $phenolabel\n";

	    if ($a1_cnt == 0) { #2 non A1 alleles, decrease phenotype value by a
		$fx_dpheno += -1*$$ref_basedata{$snp}{$phenolabel}[0]; 
	    } elsif ($a1_cnt == 1) { #heterozygote, add the domdev term
		$fx_dpheno += $$ref_basedata{$snp}{$phenolabel}[1];
	    } elsif ($a1_cnt == 2) { #2 A1 alleles, increase phenotype value by a
		$fx_dpheno += $$ref_basedata{$snp}{$phenolabel}[0];
	    }
	}
    }

    #print "Genetic log-odds augmented to $phenolabel : $fx_dpheno\n";
    #sleep 1;
    
    return(sprintf("%1.5f", $fx_dpheno));
}

sub do_pheno_distort {
    my ($ref_simgeno, $ref_basedata, $ref_infodata, $phenolabel, $ref_shared_ipheno, $asc_pheno, $out_fh, $log_fh) = @_;
    my $adj_pheno = 0;   
    my $a1_cnt;
    my $this_adj = 0;
    
    $$ref_shared_ipheno = 0;

    #@{$$ref_basedata{$snp}{$pheno}} = ($a_fx, $d_fx, $$ref_infodata{$snp}[2]);
    #@{$$ref_simgeno{$snp}} = ($allele_a, $allele_b);   
    
    #print "here\n";

    foreach my $snp (sort string_numerically keys %$ref_infodata) {
	$a1_cnt = 0;
	if (!defined($$ref_basedata{$snp}{$phenolabel}[0])) {
	    #for this phenotype, the SNP has no modelled effect.
	    #do nothing	
	} else {
	    #get the count of A1 alleles
	    if ($$ref_simgeno{$snp}[0] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    } 
	    if ($$ref_simgeno{$snp}[1] =~ m/^$$ref_infodata{$snp}[0]$/) {
		$a1_cnt++;
	    }
	    
	    #print "$snp $a1_cnt $phenolabel -- $$ref_basedata{$snp}{$phenolabel}[0]\n";
	    #sleep 1;

	    #running tally of total adjustment
	    if ($a1_cnt == 0) { #2 non A1 alleles, decrease phenotype value by a
		$this_adj = -1*$$ref_basedata{$snp}{$phenolabel}[0]; 
	    } elsif ($a1_cnt == 1) { #heterozygote, add the domdev term
		$this_adj = $$ref_basedata{$snp}{$phenolabel}[1];
	    } elsif ($a1_cnt == 2) { #2 A1 alleles, increase phenotype value by a
		$this_adj = $$ref_basedata{$snp}{$phenolabel}[0]; 
	    }
	    $adj_pheno += $this_adj;

	    #for this SNP, determine if also has a non-zero entry for the ascertained phenotype
	    #this is a slight hack (value could be 0 by user spec) - need to fix
	    #if SNP has a non-zero value AND if endpoint has non-zero value
	    if ($$ref_basedata{$snp}{$phenolabel}[0] != 0 ) { #this SNP has an impact on phenolabel
		if ( $$ref_basedata{$snp}{$asc_pheno}[0] == 0 ) { #this SNP has no impact on endpoint
		    #no ascertainment pheno for this trait. do nothing.
		} else {
		    #print "$snp $phenolabel $$ref_basedata{$snp}{$asc_pheno}[0]\n";
		    $$ref_shared_ipheno += $this_adj;
		}
	    }

	}
    }
    
    #print "Adjustment $phenolabel : $adj_pheno\n";
    #print "Shared Adjustment $phenolabel : $$ref_shared_ipheno\n";
    #$adj_pheno += $phenoval;
    #print "$adj_pheno\n";
    #sleep 1;
    
    $$ref_shared_ipheno = sprintf("%1.6f", $$ref_shared_ipheno);
    return(sprintf("%1.6f", $adj_pheno));
}

sub draw_idrel {
    my ($ref_idpheno, $ref_sim_idpheno, $out_fh, $log_fh) = @_;
    foreach my $dpheno (keys %$ref_idpheno) {
	foreach my $ipheno (keys %{$$ref_idpheno{$dpheno}}) {
	    #draw a random ipheno risk from the distribution of risks (with error).
	    $my_ipheno_fx = sprintf("%1.5f",normal($$ref_idpheno{$dpheno}{$ipheno}[0], $$ref_idpheno{$dpheno}{$ipheno}[1]));
	    @{$$ref_sim_idpheno{$dpheno}{$ipheno}} = ($my_ipheno_fx, -9); #keyed on D, I; store beta, -9 for error [simulation instance]
	}
    }
}



sub sim_dpheno {
    my ($ref_phenodata, $ref_sim_idpheno, $ref_ipheno_vals, $ref_ipheno_order, $ref_dpheno_vals, $ref_dpheno_order, $ref_gt_fx_dpheno, $ref_shared_ipheno, $ref_sim_covars, $out_fh, $log_fh) = @_;
    my ($basefactor, $my_ipheno_fx, $fold_adj, $ipheno, $dpheno, $gt_fx, $disprob, $cov_val);

    #for now, assuming type is 'b' for 0/1 dichotomous trait.
    for (my $i=0; $i<scalar(@{$ref_dpheno_order}); $i++) {
	$dpheno = $$ref_dpheno_order[$i];
	$gt_fx = $$ref_gt_fx_dpheno[$i];

	if ($$ref_phenodata{$dpheno}[1] =~ m/b/i) { #binary/dichotomous trait
	    my $Kp = $$ref_phenodata{$dpheno}[2]; #this is the background prevalence.
	    $basefactor = log($Kp) - log(1-$Kp); #logit, call this a0 (see Wu et al); this determined from the background prevalence.

	    #get the fold change to background given intermediate traits which alter the destination trait
	    for (my $i=0; $i<scalar(@{$ref_ipheno_order}); $i++) {
		$ipheno = $$ref_ipheno_order[$i];

		if (!defined($$ref_sim_idpheno{$dpheno}{$ipheno})) {
		    #this ipheno has no influnence on this dpheno
		} else { 		    
		    #draw a ipheno risk from the distribution of risks (with error).
		    $my_ipheno_fx = $$ref_sim_idpheno{$dpheno}{$ipheno}[0];

		    #subtract out SNP contribution which has a specifically enumerated endpoint allotment
		    #print "$ipheno $dpheno $my_ipheno_fx $$ref_ipheno_vals[$i] $$ref_shared_ipheno[$i]\n";
		    #sleep 1;

		    $basefactor += $my_ipheno_fx * ($$ref_ipheno_vals[$i] - $$ref_shared_ipheno[$i]);

		}
	    }	    

	    #add the effect of SNPs with direct impact on the destination trait
	    $basefactor += $gt_fx;

	    #add the effect of SEX and other covariates to the prob of disease
	    foreach my $covar (keys %$ref_sim_covars) {
		$cov_val = $$ref_sim_covars{$covar};
		$my_ipheno_fx = $$ref_sim_idpheno{$dpheno}{$covar}[0];
		if ($covar =~ m/^SEX$/) {
		    $cov_val--;
		} 
		
		$basefactor += $my_ipheno_fx * $cov_val;
		
		#print "$covar $cov_val $my_ipheno_fx\n";
		#sleep 1;
	    }

	    #implement logit function
	    $disprob = 1 / (1 + exp(-1*$basefactor));

	    #print "$Kp $basefactor $gt_fx $disprob\n";
	    #sleep 1;

	    if ($disprob > 1.0) { #should never get here
		$myprint = "WARNING: prob of disease for this simulated individual is greater than 1!\n";
		print_string($myprint, $out_fh, $log_fh);
	    }

	    if ( rand() < $disprob ) {
		push @{$ref_dpheno_vals}, 2; #got a case
	    } else { 
		push @{$ref_dpheno_vals}, 1; #is a control
	    }
	} elsif ($$ref_phenodata{$dpheno}[2] =~ m/q/i) { #outcome trait is quantitative
	    #not implemented.
	} 
    }

}

#Defaults
$specify_seed = 0; #default initialize of srand(), user can set to specific seed on command line.
$covars_exist = 0; #default assumes that no covariates have been included in the phenotype listing
$got_outfix = 0; #assume default output file name
$do_cc_sim = 1; #assume case/control destination variable.
$do_hazard_sim = 0; #default is case/control ascertainment (not prospective cohort)
$do_excl_meanadj = 0; #assume all traits apply a mean adjustment
$do_specify_var = 0; #assume all intermediate traits have DEFAULT_VAR variance associated with them
$do_printscores = 0; #default is not to print score-based statistics
$do_phenoverbose = 0; #default is a simple output of phenotypes
$do_plink_testing = 1; #default is to perform association testing in plink
@excl_meanadj_phenolist = (); #default initialization

$seed = "";

$phenofile = "";
$iirelfile = "";
$idrelfile = "";
$infofile = "";
$scorefile = "";
$asc_pheno = "";

$pheno_sim_warning = 0;

$sim_n_fold_pheno = 3; #fold factor to generate more phenotypes ($sim_n_fold_pheno x ncases) to cover rejection sampling when selecting cases 
$ncases = 1000; #default case sample size
$ncntls = 1000; #default controls sample size
$nsims = 10; #number of simulations to perform

#genetic parameters
sub DEFAULT_VAR() { 1 } ; #default variance of 1 for intermediate traits -- note that if this is modified, the covariance calculations will break down
sub SEX_EFFECT() { 0 } ; #default effect of SEX set to OR=1, ln(OR)=0 
sub SEX_PREV() { 0.5 } ; #default prob that individual is female (female = 1, male = 0; in plink coding, female = 2, male = 1);
sub NUM_INTG_BIN() { 0.01 } ; #default window to perform numerical integration for kp estimation.


#set up the log file
#make sure you have an out file specified. if tried (and failed), die.
for($i=0; $i<scalar(@ARGV); $i++) {
    $arg = $ARGV[$i];
    if ($arg =~ m/--out$/) { #outfile prefix
	$got_outfix = 1;
	$outfix = $ARGV[$i+1];
	if (!defined($outfix)) {
	    $mylogerr = "ERROR: Missing argument for " . $arg. "\n";
	    exit();
	} else {
	    if ($outfix =~ m/--/) {
		$mylogerr = "ERROR: Invalid argument for " . $arg . " [ " . $outfix . " ] (Could an argument be missing?)\n";
		exit();
	    }
	}
    }
}

if ($got_outfix == 0) {
    $outfix = "mr_pred";
}

#open the log file and set output handles.
$mylog = LOG;
$myout = STDOUT;
open $mylog, ">@{[$outfix]}.log" or die "Error: Can't write logfile to @{[$outfix]}.log. Exiting!\n";

#print the header splash.
print_header($myout);
print_header($mylog);

$myprint = "Writing this text to log file [ " . $outfix . ".log ]\n";
print_string($myprint, $myout, $mylog);

#print start time.
$myprint = "Analysis started:\t" . scalar(localtime) . "\n\n";
print_string($myprint, $myout, $mylog);

if (scalar(@ARGV) == 0) {
    $myprint = "No arguments Specified. Exiting!\n";
    print_string($myprint, $myout, $mylog);
    exit();
} else {

    $cmd_string = "Options in effect:\n\n";
    while (scalar(@ARGV) > 0) {
        $arg = shift(@ARGV);
        $cmd_string = $cmd_string . "\t" . $arg . " ";
	if ($arg =~ m/--score$/) { #score sheet
	    $scorefile = shift(@ARGV);
            if (!defined($scorefile)) {
                $myprint = "ERROR: Missing argument for $scorefile\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($scorefile =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $scorefile . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $scorefile . "\n";
                }
            }
	} elsif ($arg =~ m/--info$/) { #info sheet
	    $infofile = shift(@ARGV);
            if (!defined($infofile)) {
                $myprint = "ERROR: Missing argument for $infofile\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($infofile =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $infofile . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $infofile . "\n";
                }
            }
	} elsif ($arg =~ m/--pheno$/) { #pheno sheet
            $phenofile = shift(@ARGV);
            if (!defined($phenofile)) {
                $myprint = "ERROR: Missing argument for $phenofile\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($phenofile =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $phenofile . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $phenofile . "\n";
                }
            }
        } elsif ($arg =~ m/--id-rel$/) { #intermediate->destination phenotype relationship sheet
            $idrelfile = shift(@ARGV);
            if (!defined($idrelfile)) {
                $myprint = "ERROR: Missing argument for $idrelfile\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($idrelfile =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $idrelfile . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $idrelfile . "\n";
                }
            }
        } elsif ($arg =~ m/--ii-rel$/) { #intermediate->intermiediate phenotype relationship sheet
            $iirelfile = shift(@ARGV);
            if (!defined($iirelfile)) {
                $myprint = "ERROR: Missing argument for $iirelfile\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($iirelfile =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $iirelfile . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $iirelfile . "\n";
                }
            }
        } elsif ($arg =~ m/--out$/) {
	    $outfix = shift(@ARGV);
            if (!defined($outfix)) {
                $myprint = "ERROR: Missing argument for " . $arg . "\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($outfix =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $outfix . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } else {
                    $cmd_string = $cmd_string . $outfix . "\n";
                }
            }
	} elsif ($arg =~ m/--exclude-mean-adj$/) { #list of traits which NOT to apply mean adjustment (default is to apply mean adjustment to all)
	    $traits = shift(@ARGV);
	    if (!defined($traits)) {
		$myprint = "ERROR: Missing argument for " . $arg . "\n";
		print_string($myprint, $myout, $mylog);
                exit();
	    } else {
                if ($traits =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $traits . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
		} else {
		    $do_excl_meanadj = 1;
                    $cmd_string = $cmd_string . $traits . "\n";
                }
	    }

	} elsif ($arg =~ m/--nsims$/) {
	    $nsims = shift(@ARGV);
            if (!defined($nsims)) {
                $myprint = "ERROR: Missing argument for " . $arg . "\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($nsims =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $nsims . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } elsif ($nsims !~ m/\d+/ ) {
		    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $nsims . " ] (Argument requires an integer.)\n";
		    print_string($myprint, $myout, $mylog);
		    exit();
		} else {
                    $cmd_string = $cmd_string . $nsims . "\n";
                }
            }
	} elsif ($arg =~ m/--cc_nsamp$/) { 
	    $asc_pheno = shift(@ARGV);
	    $ncases = shift(@ARGV);
	    $ncntls = shift(@ARGV);
	    if (!defined($asc_pheno) || !defined($ncases) || !defined($ncntls)) {
                $myprint = "ERROR: Missing argument for " . $arg . "\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($asc_pheno =~ m/--/ || $ncases =~ m/--/ || $ncntls =~ m/--/) {
                    $myprint = "ERROR: Invalid arguments for " . $arg . " [ " . $asc_pheno . " " . $ncases . " " . $ncntls . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } elsif ($ncases !~ m/\d+/ || $ncntls !~ m/\d+/) {
		    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $asc_pheno . " " . $ncases . " " . $ncntls . " ] (Argument requires an integer for cases/control numbers.)\n";
		    print_string($myprint, $myout, $mylog);
                    exit();
		} elsif ($ncntls < 0 || $ncases < -1) {
		    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $asc_pheno . " " . $ncases . " " . $ncntls . " ] (Argument requires a non-negative integer for controls and cases; -1 tolerated for cases for ascertainment-free sampling.)\n";
		    print_string($myprint, $myout, $mylog);
                    exit();
		} elsif ($do_hazard_sim == 1) {
		    $myprint = "ERROR: Case/Control ascertainment not compatible with prospective cohort (--hazard) simulation.\n";
		    print_string($myprint, $myout, $mylog);
                    exit();		    
		} else {
		    $do_cc_sim = 1;
                    $cmd_string = $cmd_string . $asc_pheno . " " . $ncases . " " . $ncntls . "\n";		    
                }
            }
	} elsif ($arg =~ m/--hazard$/) {
	    $nsamp = shift(@ARGV);
	    if (!defined($nsamp)) {
                $myprint = "ERROR: Missing argument for " . $arg . "\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($nsamp =~ m/--/) {
                    $myprint = "ERROR: Invalid arguments for " . $arg . " [ " . $nsamp . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } elsif ($nsamp !~ m/\d+/) {
		    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $nsamp . " ] (Argument requires an integer for cases/control numbers.)\n";
		    print_string($myprint, $myout, $mylog);
                    exit();
		} elsif ($do_cc_sim == 1) {
		    $myprint = "ERROR: Prospective cohort ascertainment not compatible with case/control ascertainment option (--cc_samp) simulation.\n";
		    print_string($myprint, $myout, $mylog);
                    exit();		    
		} else {
		    $do_hazard_sim = 1;		    
                    $cmd_string = $cmd_string . $nsamp . "\n";
                }
            }
	} elsif ($arg =~ m/--set-ipheno-var$/) {
	    $do_specify_var = 1;
	    $cmd_string = $cmd_string . "\n";
	} elsif ($arg =~ m/--scores$/) {
	    $do_printscores = 1;
	    $cmd_string = $cmd_string . "\n";
	} elsif ($arg =~ m/--pheno-verbose$/) {
	    $do_phenoverbose = 1;
	    $cmd_string = $cmd_string . "\n";
	} elsif ($arg =~ m/--seed$/) {
	    $specify_seed = 1;
	    $seed = shift(@ARGV);
            if (!defined($seed)) {
                $myprint = "ERROR: Missing argument for " . $arg . "\n";
                print_string($myprint, $myout, $mylog);
                exit();
            } else {
                if ($seed =~ m/--/) {
                    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $seed . " ] (Could an argument be missing?)\n";
                    print_string($myprint, $myout, $mylog);
                    exit();
                } elsif ($seed !~ m/\d+/ ) {
		    $myprint = "ERROR: Invalid argument for " . $arg . " [ " . $seed . " ] (Argument requires an integer.)\n";
		    print_string($myprint, $myout, $mylog);
		    exit();
		} else {
                    $cmd_string = $cmd_string . $seed . "\n";
                }
            }
	} elsif ($arg =~ m/--skip-plink$/) {
	    $do_plink_testing = 0;
	    $cmd_string = $cmd_string . "\n";
	} elsif ($arg =~ m/--qtl_nsamp$/) {
	    #not defined yet.
	} else {
	    $myprint = "Improper Argument Specified: " . $arg. "\n";
            print_string($myprint, $myout, $mylog);
            exit();
	}
    }
}


$cmd_string = $cmd_string . "\n";
print_string($cmd_string, $myout, $mylog);

#set seed if applicable
if ($specify_seed == 1) {
    $myprint = "Seeding random number generator with specified seed [ $seed ]\n";
    print_string($myprint, $myout, $mylog);
    srand($seed); #explicitly select the seed.    
} else {
    srand(); #default call;
} 

#initialize the amount of phenotype data to simulate at a time
if ($ncases != -1) {
    $ipheno_chunksize = ($ncases*$sim_n_fold_pheno)+$ncntls;
} else {
    $ipheno_chunksize = $ncntls; 
}

#basic functionality always requires a phenotype file and scoresheet file.
#test if file exists -- die if failure.
check_file_exists($infofile, $myout, $mylog);
check_file_exists($scorefile, $myout, $mylog);	    
check_file_exists($phenofile, $myout, $mylog);	    
check_file_exists($iirelfile, $myout, $mylog);
check_file_exists($idrelfile, $myout, $mylog);

get_infofile($infofile, \%info, $myout, $mylog);
get_phenofile($phenofile, \%pheno, \$covars_exist, $myout, $mylog);
get_iiphenofile($iirelfile, \%iipheno, \%pheno, $outfix, $do_specify_var, $myout, $mylog);
get_idphenofile($idrelfile, \%idpheno, \%pheno, \$asc_pheno, $myout, $mylog);
get_scorefile($scorefile, \%snps, $myout, $mylog);

if ($do_cc_sim == 1) { #perform case/control destination simulation.
    #check that you have specified an the ascertainment phenotype (and that it is of type boolean).
    if (!defined($pheno{$asc_pheno}[0])) {
	$myprint = "ERROR: Ascertainment phenotype $asc_pheno not specified in phenotype file. Exiting!\n";
	print_string($myprint, $myout, $mylog);
	exit();
    } else {
	if ($pheno{$asc_pheno}[0] !~ m/d/i) {
	    $myprint = "ERROR: Ascertainment phenotype $asc_pheno must be a destination phenotype in the phenotype file. Exiting!\n";
	    print_string($myprint, $myout, $mylog);
	    exit();
	} elsif ($pheno{$asc_pheno}[1] !~ m/b/) {
	    $myprint = "ERROR: Ascertainment phenotype $asc_pheno must be set to a dichotomous trait in the phenotype file. Exiting!\n";
	    print_string($myprint, $myout, $mylog);
            exit();
	}
    }

    #log parameters to file.
    if ($ncases != -1) {
	$myprint = "Assuming sample size of $ncases cases and $ncntls controls.\n";
	print_string($myprint, $myout, $mylog);
    } else {
	$myprint = "Assuming sample size of $ncntls.\n";
	print_string($myprint, $myout, $mylog);
    }

    $myprint = "Performing $nsims simulations.\n";
    print_string($myprint, $myout, $mylog);

    #calculate sigma_r baseline.
    mk_basedata(\%snps, \%base, \%info, \%pheno, 0, $myout, $mylog);

    foreach $ipheno (sort string_numerically keys %pheno) {
	if ($pheno{$ipheno}[0] =~ m/i/i) {
	    $sig_r = sprintf("%1.4f", get_sigma_r(\%base, $ipheno));
	    $myprint = "Baseline variance for known loci on intermediate trait $ipheno [baseline]: $sig_r\n";
	    print_string($myprint, $myout, $mylog);
	}
    }

    foreach $dpheno (sort string_numerically keys %pheno) {
	if ($pheno{$dpheno}[0] =~ m/d/i) {
	    $sig_r = sprintf("%1.4f", get_sigma_r(\%base, $dpheno));
	    $myprint = "Baseline variance for known loci on destination trait $dpheno [baseline]: $sig_r\n";
	    print_string($myprint, $myout, $mylog);
	}
    }

    #create a generic map file for use with plink
    mk_mapfile(\%snps, $outfix, $myout, $mylog);

    for ($N=1; $N<=$nsims; $N++) {
	undef(%pheno_mean_adj);	
	@sim_vcovm = ();
	@sim_meanvec = ();
	@sim_ipheno_order = ();

	$myprint = "#### Performing Simulation $N\n";
	print_string($myprint, $myout, $mylog);

	$mysimped = SIMPED;
	$mysimpheno = SIMPHENO;
	open $mysimped, ">@{[$outfix]}_@{[$N]}.ped" or die "Can't open @{[$outfix]}_@{[$N]}.ped!\n";
	open $mysimpheno, ">@{[$outfix]}_@{[$N]}.pheno" or die "Can't open @{[$outfix]}_@{[$N]}.pheno!\n";

	if ($do_printscores == 1) {
	    $mysimscore = SIMSCORE;
	    open $mysimscore, ">@{[$outfix]}_@{[$N]}.score" or die "Can't open @{[$outfix]}_@{[$N]}.score!\n";
	}
	
	if ($covars_exist == 1) {
	    $mysimcovars = SIMCOVARS;
	    open $mysimcovars, ">@{[$outfix]}_@{[$N]}.cov" or die " Can't open @{[$outfix]}_@{[$N]}.cov!\n";
	}

	if ($N == 1) {
	    #Reading in phenotype info from [ " . $phenofile . " ]\n";
	    $myprint = "Outputting simulated pedigree data to [ " . $outfix . "_" . $N . ".ped ]\n";
	    print_string($myprint, $myout, $mylog);

	    $myprint = "Outputting simulated phenotype data to [ " . $outfix . "_" . $N . ".pheno ]\n";
            print_string($myprint, $myout, $mylog);

	    if ($do_printscores == 1) {
		$myprint = "Outputting simulated phenotype score data to [ " . $outfix . "_" . $N . ".score ]\n";
		print_string($myprint, $myout, $mylog);		
	    } 

	    if ($covars_exist == 1) {
		$myprint = "Outputting simulated covariate data to [ " . $outfix . "_" . $N . ".cov ]\n";
		print_string($myprint, $myout, $mylog);
	    }

	    $myprint = "Supressing further output file information...\n";
	    print_string($myprint, $myout, $mylog);	    	    
	} 

	#calculate random sigma_r from set
	mk_basedata(\%snps, \%base, \%info, \%pheno, 1, $myout, $mylog);
	foreach $ipheno (sort string_numerically keys %pheno) {
	    if ($pheno{$ipheno}[0] =~ m/i/i) {
		$sig_r = sprintf("%1.4f", get_sigma_r(\%base, $ipheno));
		$myprint = "Baseline variance for known loci on intermediate trait for $ipheno Sim $N: $sig_r\n";
		print_string($myprint, $myout, $mylog);
	    }
	}

	foreach $dpheno (sort string_numerically keys %pheno) {
	    if ($pheno{$dpheno}[0] =~ m/d/i) {
		$sig_r = sprintf("%1.4f", get_sigma_r(\%base, $dpheno));
		$myprint = "Baseline variance for known loci on destination trait $dpheno [baseline]: $sig_r\n";
		print_string($myprint, $myout, $mylog);
	    }
	}
	
	#calculate the mean adjustment effect
        #list traits excluded from mean adjument
	if ($do_excl_meanadj == 1) { 
	    @excl_meanadj_phenolist = split ",", $traits; 
	    $size = scalar(@excl_meanadj_phenolist);

	    #To do: check if exclude trait lists matches given phenotypes

	    $myprint = "Excluding $size intermediate traits from mean adjustment.\n";
	    print_string($myprint, $myout, $mylog);
	} 
	calc_mean_adj(\%base, \%pheno_mean_adj, \%pheno, \@excl_meanadj_phenolist, \@sim_meanvec, \@sim_ipheno_order);
	#print "Adjustment for LDL: $pheno_mean_adj{LDL}\n";

	#foreach $pheno (sort string_numerically keys %pheno_mean_adj) {
	#    print "$pheno: $pheno_mean_adj{$pheno}\n";
	#}
        #
	#foreach $val (@sim_meanvec) {
	#    print "$val ";
	#}
	#print "\n";
	#exit();
	
	#foreach $val (@sim_ipheno_order) {
	#    print "$val\n";
	#}
	#exit();

	#calculate the covar adjustment matrix for SNPs which contribute to multiple traits.
	$myprint = "Calculating covariance adjustments to baseline for Sim $N...\n";
	print_string($myprint, $myout, $mylog);
	get_cov_adj_matrix(\%base, \%cov_adj, \%pheno, $myout, $mylog);

        #Construct a var/covar file for the I-I, and the reduced variances for each I 	
	$rvarcovfile = $outfix . "_varcov_info_" . $N . ".forR";
	mk_varcovfile($rvarcovfile, \@sim_vcovm, \%iipheno, \%pheno, \%base, \%cov_adj, $myout, $mylog);
	
	#for ($i=0; $i<scalar(@sim_vcovm); $i++) {
	#    for ($j=0; $j<scalar(@{$sim_vcovm[$i]}); $j++) {
	#	print "$sim_vcovm[$i][$j] ";
	#    }
	#    print "\n";
	#}
	
	#construct a random I-D draw from this simulation
	draw_idrel(\%idpheno, \%sim_idpheno, $myout, $mylog);
	if ($N == 1) {
	    #Reading in phenotype info from [ " . $phenofile . " ]\n";
	    foreach $dpheno (keys %idpheno) {
		foreach $ipheno (keys %{$idpheno{$dpheno}}) {
		    $myprint = "Baseline effect of $ipheno -> $dpheno: $idpheno{$dpheno}{$ipheno}[0] $idpheno{$dpheno}{$ipheno}[1]\n";
		    print_string($myprint, $myout, $mylog);
		}
	    }
	    foreach $dpheno (keys %sim_idpheno) {
		foreach $ipheno (keys %{$sim_idpheno{$dpheno}}) {
		    $myprint = "Simulated effect for $ipheno -> $dpheno for Sim $N: $sim_idpheno{$dpheno}{$ipheno}[0]\n";
		    print_string($myprint, $myout, $mylog);
		}
	    }
	    $myprint = "Supressing further simulation data output information...\n";
	    print_string($myprint, $myout, $mylog);
	}

	##################
	#BEGIN SIMULATION#
	##################

	$continue = 1;
	$n_sim_cases = 0;
	$n_sim_controls = 0;
	$nlabel = 0;
	@sim_ipheno_data = ();
	@dis_ipheno = ();
	@shared_ipheno = ();

	while ($continue == 1) {
	    $printme = 0;
	    @sim_dpheno_data = ();
	    @sim_dpheno_order = ();
	    @gt_fx_dpheno = ();
	    undef(%scorestats);
	    undef(%sim_verbosepheno);
	    undef(%sim_covars);

	    #simulate a genotype vector from the given SNP map for the given individual 
	    sim_genotypes(\%info, \%genodata);

	    #simulate covariate vector for the given individual
	    sim_covars(\%pheno, \%sim_covars);

	    #simulate intermediate trait vector for the given individual
	    #simulate a chunk of data to start with [minimize subroutine calls]
	    if (scalar(@sim_ipheno_data) == 0) {
		if ($pheno_sim_warning == 0) {
		    #$myprint = "Simulating intermediate trait data as specified in the Variance/Covariance Matrix file [ $rvarcovfile ]\n";
		    #print_string($myprint, $myout, $mylog);
		    $myprint = "Supressing further intermediate phenotype simulation output...\n";
		    print_string($myprint, $myout, $mylog);
		    $pheno_sim_warning = 1;
		}
		#generate the chunk of data, minimize functional calls.
		@sim_ipheno_data = Math::Random::random_multivariate_normal($ipheno_chunksize, @sim_meanvec, @sim_vcovm);
	    }
	    
	    #take a sample off the simulated phenotype list.
	    $this_data = shift(@sim_ipheno_data);
	    
	    #print "$sim_ipheno_data[0][0] $sim_ipheno_data[0][1]\n";
	    #print "$sim_ipheno_data[10][0] $sim_ipheno_data[10][1]\n";
	    #exit();

	    #output pheno (testing).
	    #print "$this_data->[0] $this_data->[1]\n";
	    #exit();

	    #Modified the phenotype based on genetic data. 
	    for ($p=0; $p<scalar(@sim_ipheno_order); $p++) {
		#get the amount of ipheno contribution where the SNP is also tagged with the outcome
		$gen_distort = do_pheno_distort(\%genodata, \%base, \%info, $sim_ipheno_order[$p], \$shared_ipheno[$p], $asc_pheno, $myout, $mylog); 
		#print "$gen_distort # $this_data->[$p]\n";

		$dis_ipheno[$p] = sprintf("%1.6f", $gen_distort + $this_data->[$p]); # - $pheno_mean_adj{$sim_ipheno_order[$p]});

		if ($do_phenoverbose == 1) {
		    @{$sim_verbosepheno{$sim_ipheno_order[$p]}} = (sprintf("%1.6f",$this_data->[$p]), $gen_distort);
		}
		
		#print "$sim_ipheno_order[$p] # vE $this_data->[$p] vG $gen_distort mu_adj $pheno_mean_adj{$sim_ipheno_order[$p]} shared $shared_ipheno[$p]\n";
		#sleep 1;		
	    }

	    #make a sim_dpheno_order list
	    $iter = 0;
	    foreach $i (keys %pheno) {
		if ($pheno{$i}[0] =~ m/d/i) {
		    $sim_dpheno_order[$iter] = $i;
		    $iter++;
		}
	    }

	    #generate the effect of SNPs which direct impact outcome based on genetic data
	    for ($p=0; $p<scalar(@sim_dpheno_order); $p++) {
		$gt_fx_dpheno[$p] = do_dpheno_fx(\%genodata, \%base, \%info, $sim_dpheno_order[$p], $myout, $mylog);
	    }

	    #now do destination phenotypes based on iphenos, etc.
	    sim_dpheno(\%pheno, \%sim_idpheno, \@dis_ipheno, \@sim_ipheno_order, \@sim_dpheno_data, \@sim_dpheno_order, \@gt_fx_dpheno, \@shared_ipheno, \%sim_covars, $myout, $mylog);

	    #$size = scalar(@sim_dpheno_order);
	    #print "-- $size\n";
	    #sleep 1;

	    #asc_pheno is the name (e.g. MI, CAD, etc.) taken from the command line
	    #this_pheno is the actual value of the endpoint
	    for ($i=0; $i<scalar(@sim_dpheno_order); $i++) {
		if ($sim_dpheno_order[$i] =~ /$asc_pheno/) {
		    $this_pheno = $sim_dpheno_data[$i];
		}
	    }

	    ########why doesn't this just output nsample individuals regardless of trait?

	    #output informtion
	    if ($ncases != -1) {
		if ($n_sim_cases < $ncases && $this_pheno == 2) {
		    #print to .ped file the case info
		    #print to .pheno file the phenotype for this sample
		    $printme = 1;
		    $n_sim_cases++;
		    $nlabel++;		
		} elsif ($n_sim_controls < $ncntls && $this_pheno == 1) {
		    #print to .ped file the control data
		    #print to .pheno file the phenotype for this sample
		    $printme = 1;
		    $n_sim_controls++;
		    $nlabel++;
		}
	    } else {
		$printme = 1;
		$n_sim_controls++;
		$nlabel++;
	    }
	    
	    if ($printme == 1) { #print outputs for this individual
		
		#print pedigree data
		print_simped(\%genodata, \%sim_covars, \$this_pheno, \$nlabel, $mysimped);
		
		#print covariate information if covariates beyond SEX are present
		if (scalar(keys %sim_covars) > 1) {
		    print_covars(\%sim_covars, \$nlabel, $mysimcovars);
		}

		#print phenotype data
		if ($do_phenoverbose == 1) {
		    print_simverbpheno(\@sim_ipheno_order, \@sim_dpheno_order, \@dis_ipheno, \@sim_dpheno_data, \%sim_verbosepheno, \$nlabel, $mysimpheno);
		} else {
		    print_simpheno(\@sim_ipheno_order, \@sim_dpheno_order, \@dis_ipheno, \@sim_dpheno_data, \$nlabel, $mysimpheno);
		}
				
		#print score info in requested
		if ($do_printscores == 1) {
		    calc_scorestat(\%genodata, \%scorestats, \%info, \%base, \%snps);	   
		    print_simscore(\@sim_ipheno_order, \@sim_dpheno_order, \%scorestats, $asc_pheno, $this_pheno, \$nlabel, $mysimscore);
		}

	    } #END printme

            #continue simulating if you haven't met the ascertainment criteria
	    if ($ncases != -1) {
		if ($n_sim_cases == $ncases && $n_sim_controls == $ncntls) { 
		    $continue = 0;
		}
	    } else {
		if ($n_sim_controls == $ncntls) {
		    $continue = 0;
		}
	    }

	} #END get ncases and ncontrols for this SIM

	#close up Files from this simulation.
	close($mysimped);
	close($mysimpheno);
	
	if ($do_printscores == 1) {
	    close($mysimscore);	
	}

	if ($covars_exist == 1) {
	    close($mysimcovars);
	}
	
	#perform testing.
	if ($do_plink_testing == 1) {
	    print "performing testing in plink\n";
	    if ($covars_exist == 1) {
		system "plink --noweb --silent --ped @{[$outfix]}_@{[$N]}.ped --map @{[$outfix]}.map --pheno @{[$outfix]}_@{[$N]}.pheno --all-pheno --covar @{[$outfix]}_@{[$N]}.cov --logistic --ci 0.95 --out @{[$outfix]}_@{[$N]}";
	    } else {
		system "plink --noweb --silent --ped @{[$outfix]}_@{[$N]}.ped --map @{[$outfix]}.map --pheno @{[$outfix]}_@{[$N]}.pheno --all-pheno --assoc --ci 0.95 --out @{[$outfix]}_@{[$N]}";
	    }
	}
       
    } #END each sim

} #END do a Case/control simulation.
elsif ($do_hazard_sim == 1) {
    #not implemented
} 


#End!
$myprint = "\nAnalysis finished:\t" . scalar(localtime). "\n";
print_string($myprint, $myout, $mylog);

#./mr_predictor.pl --score scoresheet_gi_mr.txt --info infosheet_gi_mr.txt --pheno phenofile_gi_mr.txt --ii-rel iifile_mr.txt --id-rel idfile_mr.txt --nsims 2 --cc_nsamp MI 1 4999 --out test/tester

#controls 3055
#cases 2949
