import os
import sys
import string
import locale


def write_weight_files(python_basis, x2, rg, weight_file_names, runst):

    # The filter module catches the errors encountered below.
    # The program shouldn't get to this point but error checks are kept just in case.
    # If these errors should ever be encountered, the program would continue
    # without writing the weight file for the errant python_basis

    error = []
    k = 0
    for basis in python_basis:

        basis = basis.replace("'", "")
        print 'basis in write_weight_files: ', basis
        x2_type = False
        rg_type = False

        if "x2[i]" in basis:
            x2_type = True
        elif "rg[i]" in basis:
            rg_type = True
        else:
            error.append('error in variable types: no rg or x2 in basis; weight file ' +
                         weight_file_names[k] + ' not written')
            k += 1
            continue

        mask_array = []
        try:
            for i in xrange(len(x2)):
                if(eval(basis)):
                    mask_array.append(1)
                else:
                    mask_array.append(0)
        except:
            error.append('ERROR: failed to parse basis; weight file ' +
                         weight_file_names[k] + ' not written')
            k += 1
            continue
#            sys.exit()

        if x2_type:
            print 'x2 = ', x2
        elif rg_type:
            print 'rg = ', rg
        print 'mask_array = ', mask_array

        outfile = open(os.path.join(runst, weight_file_names[k]), 'w')
        for i in xrange(0, len(x2)):
            if x2_type:
                outfile.write('%i\t%f\t%f\n' % (i + 1, x2[i], mask_array[i]))
            elif rg_type:
                outfile.write('%i\t%f\t%f\n' % (i + 1, rg[i], mask_array[i]))

        outfile.close()
        k += 1

    if len(error) > 0:
        print error

    return

if __name__ == "__main__":

    python_basis = ["rg[i]  <  '10' ", "x2[i]  >  '1.0' "]

    x2 = [0.1 * float(x) for x in xrange(20)]
    rg = [float(x) for x in xrange(20)]

    weight_file_names = ['w1.txt', 'w2.txt']

    write_weight_files(python_basis, x2, rg, weight_file_names)
