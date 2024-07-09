'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import sys
import string
import locale
import sassie.interface.input_filter as input_filter


def data_check(expdata, dq, maxpoints, datatype):
    #    returns ev=1 &  value=1 if everything is okay
    #    returns ev=0 if data file does not exist
    #    returns ev=2 if data file is a directory
    #    returns value=2 if npoints can not be read from gui
    #    returns value=3 if interpolate range exceeds the range of data
    #    returns value=4 if dq is smaller than the first point

    print('checking data')
    ev = 0
    value = 0
    good = 0
    thisq = 0
    squak = 0
    readnext = 0

    last_value = -1.0E9

    if(datatype == 'ncnr'):
        try:
            fileexist = os.path.exists(expdata)
            if(os.path.isdir(expdata)):
                ev = 2
                return ev, value
            if(fileexist):
                ev = 1
                try:
                    totq = dq * maxpoints
                except:
                    value = 2
                    return ev, value
                # data=open(expdata,'r').readlines()
                data = open(expdata, 'r').read().splitlines()
                nl = len(data)
                for k in xrange(nl):
                    lin = string.split(data[k])
                    try:
                        dum = float(lin[0])
                        dum = float(lin[1])
                        dum = float(lin[2])
                        readnext = 1
                    except:
                        pass
                    if((len(lin) > 2) and readnext == 1):
                        x = locale.atof(lin[0])
                        y = locale.atof(lin[1])
                        z = locale.atof(lin[2])

                        if(z <= 0.0):
                            value = 6
#                            return ev, value

                        if(y < 0.0):
                            value = 7
#                            return ev, value

                        if(x <= last_value):
                            value = 5
                            print('this_value = ', x)
                            print('last_value = ', last_value)
                            print('ERROR')
                            return ev, value
                        else:
                            last_value = x
                        if(squak == 0):
                            eq0 = x
                            squak = 1
                        good += 1
                        lastdum1 = x

                if(totq > (lastdum1 - eq0)):
                    value = 3
                    print('first value in data = ', eq0)
                    print('last value in data = ', lastdum1)
                    print('totq = ', totq)
                    return ev, value
                if(dq < eq0):
                    print('dq = ', dq)
                    print('first_value in data = ', eq0)
                    value = 4
#                    return ev,value
        except:
            return ev, value

    return ev, value


def check_asaxs(variables, **kwargs):

    run_name = variables['run_name'][0]


    return

    expdata = variables['expdata'][0]
    io = variables['io'][0]
    ioe = variables['ioe'][0]
    dq = variables['dq'][0]
    maxpoints = variables['maxpoints'][0]

    error = []

    # if 'no_file_check' not in kwargs:
    error = input_filter.check_name(run_name)
    if(error != []):
        return error

    if(dq <= 0):
        error.append('dq needs to be > 0 : '+str(dq))
        return error
    elif(maxpoints <= 2):
        error.append('maxpoints needs to be greater than 2 : '+str(maxpoints))
        return error

    ev, value = data_check(expdata, dq, maxpoints, 'ncnr')

    if(ev == 0):
        error.append('input data file, '+expdata+', does not exist')
        return error
    elif(ev == 2):
        error.append('input data file, '+expdata+', is a directory!')
        return error
    elif(value == 2):
        error.append(
            'delta q and/or number of points can not be read from GUI (did you enter a number?) : '+str(dq)+' '+str(maxpoints))
        return error
    elif(value == 3):
        error.append('input parameters compared to data in file, '+expdata +
                     ', are not valid : dq*maxpoints > q-range in data file [q-range = last q-value minus first q-value] : dq = '+str(dq)+', maxpoints = '+str(maxpoints))
        return error
    # elif(value == 4):
        #    error.append( 'input parameters compared to data in file, '+expdata[3:]+', are not valid : dq < value of first point; choose a larger value of dq')
        #    return error
    elif(value == 5):
        error.append('in your input file, ' +
                     expdata[3:]+', found that q-values that are not increasing: duplicate points or successive q-values where q[0] < q[1]')
        return error
    # elif(value == 6):
    #    error.append( 'in your input file, '+expdata[3:]+', found error values <= 0 : all values need to be >= 0')
    #    return error
    # elif(value == 7):
    #    error.append( 'in your input file, '+expdata[3:]+', found I(q) values < 0 : all values need to be > 0')
    #    return error

    return error
