The tools currently require the ''AI::Genetic'' Perl module, it can easily
installed like so:

    perl -MCPAN -e 'install AI::Genetic'

The perl script depends on an executable ''compare'', which depends on
the the ''netpbm'' and ''openjpg'' libraries; if you're using OS X,
they can be installed like so:

    $ brew install netpbm openjpeg

The ''compare'' binary can be built with the following command line:

    $ make

To run the search for a good quantization matrix, run the following:

    $ ./quanty 30 lenaofficial.pgm

Where 30 it the JPEG quality, and the image is an image to find a good
quantization matrix for.

That's it.  One day this might actually be useful.

