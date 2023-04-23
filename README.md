# README.md

# MR_predictor: A simulation engine for Mendelian Randomization Studies

#########
#### INSTALLATION

I've created a public github repository from which the engine can be obtained.

Assuming one has git installed, on the command line one could obtain a clone via

%> git clone http://github.com/bvoight/mr_predictor

##########
#### REQUIREMENTS

MR_predictor requires:

- PERL (tested with v5.10.1)
- The Math::Random package (tested with v0.71)

To install, one could:

A. Install via CPAN (assuming this is active)

> cpan

> make install GROMMEL/Math-Random-0.71.tar.gz

One may need root access to install. Alternatively, specify a local directory for the installation.

B. MANUALLY INSTALL

1. download http://search.cpan.org/CPAN/authors/id/G/GR/GROMMEL/Math-Random-0.71.tar.gz

2. installation (without root access)

> perl Makefile.PL INSTALL_BASE=~/bin/Perllibs

> make

> make test

> make install

3. Then setup your .bash_profile to look in the correct directory, e.g.

> set PERL5LIB=$PERL5LIB:$HOME/bin/Perllibs/lib/perl5/x86_64-linux-thread-multi

> export PERL5LIB

##########
#### USAGE

The package comes with a file (cmdline) which gives the example of usage.

If everything is working properly, after downloading the git repository and installing required packages, you should be able to invoke:

%> ./cmdline

which will generate a simulated data set based on the included sample files.

#########
#### COMMENTS, FEATURE REQUESTS, AND BUG REPORTING

I will endeavor to keep a running list of changes that come in, but may not immediate get to fixing them instantly!

