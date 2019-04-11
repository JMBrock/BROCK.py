def ProgressBar(i,j,percent):
    import sys
    sys.stdout.write('\r')
    # sys.stdout.write("%i of %i[%-50s] %d%%" % (i,j, '='*int(percent/2.0), percent ) )
    sys.stdout.write("[%-50s] %d%%" % ('='*int(percent/2.0), percent ) )
    sys.stdout.flush()