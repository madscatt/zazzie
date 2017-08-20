import math
def get_ranges(ncpus, nfiles):

    if nfiles < ncpus:
        nbatches = 0
        nfiles_per_batch = nfiles
    else:
        nbatches = int(math.floor(float(nfiles) / float(ncpus)))
        nfiles_per_batch = ncpus

    #nfiles_per_batch = nfiles / nbatches

    if nbatches != 1 :
        remainder = nfiles % ncpus
    else:
        remainder = 0

    remainder = nfiles % ncpus

    print 'ncpus = ', ncpus
    print 'nfiles = ', nfiles
    print 'nbatches = ', nbatches
    print 'nfiles_per_batch = ', nfiles_per_batch
    print 'remainder = ', remainder

    first = [] ; last = [] ; number = []

    for i in xrange(nbatches):
        this_first = i * nfiles_per_batch
        first.append(this_first)
        this_last = this_first + nfiles_per_batch - 1
        last.append(this_last)
        number.append(this_last - this_first + 1)
        print 'i = ', i, ' : first = ', this_first, ' : last = ', this_last, ' : number = ', number[-1]

    if nbatches > 1: 
        final_first = this_last + 1
        final_last = final_first + remainder - 1
        final_number = final_last - final_first + 1
        print 'final : first = ', final_first, ' : last = ', final_last, ' : number = ', final_number

    else:
        final_first = 0
        final_last = remainder - 1
        final_number = remainder
        print 'final : first = ', final_first, ' : last = ', final_last, ' : number = ', final_number

    if final_number > 0:
        first.append(final_first)
        last.append(final_last)
        number.append(final_number)

    return first, last, number

if __name__ == "__main__":
    ncpus = 64
    nfiles = 1*64 + 2
    #nfiles = 2

    first, last, number = get_ranges(ncpus, nfiles)

    print '-' * 50
    print '-' * 50

    print 'first = ', first
    print 'last = ', last
    print 'number = ', number ; print 'sum number = ', sum(number)

    print '-' * 50
    print '-' * 50
