# -*- coding: utf-8 -*-
import os
import glob
import numpy
import locale
import string
import math

BOKEH_PLOT = False
MATPLOTLIB_PLOT = True
#MATPLOTLIB_PLOT = False
PLOTALL = False

if MATPLOTLIB_PLOT:
    import matplotlib.pyplot as plt

if BOKEH_PLOT:

    from bokeh.layouts import column
    from bokeh.models import ColumnDataSource, Slider
    from bokeh.plotting import figure
    from bokeh.server.server import Server
    from bokeh.themes import Theme
    from bokeh.sampledata.sea_surface_temperature import sea_surface_temperature

DEBUG = True

if DEBUG:
    import sys
    sys.path.append('./')
    import file_utils as file_utils
else:
    import sassie.analyze.eros.file_utils as file_utils

def update_status(other_self, fraction_done):

    log = other_self.log
    pgui = other_self.run_utils.print_gui

    report_string = 'STATUS\t%f' % fraction_done
    pgui(report_string)
    log.debug(report_string)

    return

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        bar_length  - Optional  : character length of bar (Int)

        from: https://gist.github.com/aubricus/f91fb55dc6ba5557fbab06119420dd6a
    """
    str_format = "{0:." + str(decimals) + "f}"
    percents = str_format.format(100 * (iteration / float(total)))
    filled_length = int(round(bar_length * iteration / float(total)))
    bar = '█' * filled_length + '-' * (bar_length - filled_length)

    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()

    return

def initialize_data(other_self):

    log = other_self.log
    pgui = other_self.run_utils.print_gui
    evars = other_self.evars

    file_utils.read_iq_data(other_self)

    file_utils.read_goal_iq_data(other_self)

    file_utils.validate_data(other_self)

    update_status(other_self, 0.1)

    evars.number_q_points = len(evars.goal_q_data)

    log.info('evars.number_q_points = ' + str(evars.number_q_points))
    pgui('evars.number_q_points = ' + str(evars.number_q_points))

    return

def calculate_average_x2(other_self, weights):
    
    evars = other_self.evars
    mvars = other_self.mvars

    tspec = numpy.zeros(evars.number_q_points, numpy.float32)

    for i in xrange(len(evars.iq_data)):
        for j in xrange(evars.number_q_points):
            tspec[j] += weights[i] * evars.iq_data[i][j]

    tspec = mvars.io * tspec
    evars.average_spec = tspec
    #print 'tspec = ', tspec
    goal = evars.goal_iq_data
    diff = tspec - goal
    diff2 = diff * diff
    sigma = evars.goal_iq_error_data
    sigma2 = sigma * sigma
    ngpoints = len(sigma)

    if(mvars.reduced_x2 == 1):
        #   reduced X2 (assuming n = 0)
        X2 = numpy.sum(diff2 / sigma2) / (ngpoints - 1)
    elif(mvars.reduced_x2 == 0):
        X2 = numpy.sum(diff2 / sigma2)  # chi-square distribution
    elif(mvars.reduced_x2 == 2):
        X2 = numpy.sum(diff2 / goal)  # Pearson's X2 test-statistic
    elif(mvars.reduced_x2 == 3):
        X2 = numpy.sum(numpy.abs(diff)) / numpy.sum(numpy.abs(goal)) # R-factor

    return X2

def calculate_x2(other_self, i):

    evars = other_self.evars
    mvars = other_self.mvars

    tspec = numpy.zeros(evars.number_q_points, numpy.float32)

    for i in xrange(len(evars.iq_data)):
        for j in xrange(evars.number_q_points):
            tspec[j] += evars.weights[i] * evars.iq_data[i][j]

    tspec = mvars.io * tspec
    evars.initial_average_spec = tspec
    
    tspec = mvars.io * evars.iq_data[i]
    goal = evars.goal_iq_data

    diff = tspec - goal
    diff2 = diff * diff

    sigma = evars.goal_iq_error_data
    sigma2 = sigma * sigma
    ngpoints = len(sigma)

    if(mvars.reduced_x2 == 1):
        #   reduced X2 (assuming n = 0)
        X2 = numpy.sum(diff2 / sigma2) / (ngpoints - 1)
    elif(mvars.reduced_x2 == 0):
        X2 = numpy.sum(diff2 / sigma2)  # chi-square distribution
    elif(mvars.reduced_x2 == 2):
        X2 = numpy.sum(diff2 / goal)  # Pearson's X2 test-statistic
    elif(mvars.reduced_x2 == 3):
        X2 = numpy.sum(numpy.abs(diff)) / numpy.sum(numpy.abs(goal)) # R-factor

    return X2

def modify_doc(doc):
    df = sea_surface_temperature.copy()
    source = ColumnDataSource(data=df)

    plot = figure(x_axis_type='datetime', y_range=(0, 25), y_axis_label='Temperature (Celsius)',
                  title="Sea Surface Temperature at 43.18, -70.43")
    plot.line('time', 'temperature', source=source)

    def callback(attr, old, new):
        if new == 0:
            data = df
        else:
            data = df.rolling('{0}D'.format(new)).mean()
        source.data = ColumnDataSource(data=data).data

    slider = Slider(start=0, end=30, value=0, step=1, title="Smoothing by N Days")
    slider.on_change('value', callback)

    doc.add_root(column(slider, plot))

    doc.theme = Theme(filename="theme.yaml")
    
    return

def modify_doc_new(doc):
    ### TODO: can't get other_self initialized this way ... 
    ### see: https://stackoverflow.com/questions/44338066/how-to-reuse-a-bokeh-app-with-different-data-on-an-embedded-bokeh-server

    evars = other_self.evars
    q_data = evars.goal_q_data
    iq_data = evars.goal_iq_data

    data_dict = {"q_data":q_data, "iq_data":iq_data}

    import pandas as pd
    pd_data = pd.DataFrame(data_dict, columns = ['q_data', 'iq_data'])

    source = ColumnDataSource(data=pd_data)

    plot = figure(y_axis_label='q (a/A)',
                  title="Goal Data")
    plot.line('q', 'iq', source=source)

    def callback(attr, old, new):
        if new == 0:
            data = pd_data
        else:
            data = df.rolling('{0}D'.format(new)).mean()
        source.data = ColumnDataSource(data=data).data

    #slider = Slider(start=0, end=30, value=0, step=1, title="Smoothing by N Days")
    #slider.on_change('value', callback)

    #doc.add_root(column(slider, plot))

    doc.theme = Theme(filename="theme.yaml")
    
    return

def create_sorted_weight_data(other_self):

    mvars = other_self.mvars
    evars = other_self.evars

    sorted_indices = numpy.argsort(-evars.weights)
    evars.sorted_weights = evars.weights[sorted_indices]
    evars.sorted_profiles = evars.iq_data[sorted_indices]
    
    print 'sw 1 = ', evars.sorted_weights[0]
    print 'max(evars.weights) = ', max(evars.weights)

    discrepency_and_weights_array = []

    for i in xrange(len(evars.sorted_weights)):
        cluster_weights = evars.sorted_weights[0:i+1]
        cluster_profile = numpy.zeros(evars.number_q_points, numpy.float32)

        for j in xrange(i+1):
            for k in xrange(evars.number_q_points):
                cluster_profile[k] += mvars.io * (cluster_weights[j]/numpy.sum(cluster_weights)) * evars.sorted_profiles[j][k]

        goal = evars.goal_iq_data
        diff = cluster_profile - goal
        diff2 = diff * diff
        sigma = evars.goal_iq_error_data
        sigma2 = sigma * sigma
        ngpoints = len(sigma)

        if(mvars.reduced_x2 == 1):
            #   reduced X2 (assuming n = 0)
            X2 = numpy.sum(diff2 / sigma2) / (ngpoints - 1)
        elif(mvars.reduced_x2 == 0):
            X2 = numpy.sum(diff2 / sigma2)  # chi-square distribution
        elif(mvars.reduced_x2 == 2):
            X2 = numpy.sum(diff2 / goal)  # Pearson's X2 test-statistic
        elif(mvars.reduced_x2 == 3):
            X2 = numpy.sum(numpy.abs(diff)) / numpy.sum(numpy.abs(goal)) # R-factor

        discrepency_and_weights_array.append([numpy.sum(cluster_weights), X2])

    evars.discrepency_and_weights_array = numpy.array(discrepency_and_weights_array)

    return


def plotme_matplotlib(other_self):

    mvars = other_self.mvars
    evars = other_self.evars
    goal_q_data = evars.goal_q_data
    goal_iq_data = evars.goal_iq_data
    goal_iq_error_data = evars.goal_iq_error_data

    my_dpi=96

    fig = plt.figure(figsize=(800/my_dpi, 1400/my_dpi), dpi=my_dpi)
    fig.subplots_adjust(hspace=0.5, wspace=0.4)

    ax1 = fig.add_subplot(7,2,(1,3))
    if PLOTALL:
        for g in xrange(len(evars.iq_data)):
            ax1.plot(evars.q_data, mvars.io * evars.iq_data[g])

    #ax1.plot(evars.goal_q_data, evars.goal_iq_data, 'bo', label='goal')
    ax1.errorbar(evars.goal_q_data, evars.goal_iq_data, yerr=goal_iq_error_data, fmt='.-', label='goal')
    ax1.plot(evars.goal_q_data, evars.initial_average_spec, 'go', label='initial equal weights')
    ax1.plot(evars.goal_q_data, evars.average_spec, 'ro-', linewidth=2, markersize=3, label='final optimized fit')

    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_title("I(q) Profiles")
    ax1.set_xlabel("q (1/$\AA$)")
    ax1.set_ylabel("I(q)")

    ax1r = fig.add_subplot(7,2,(2,4), sharex=ax1)
    ax1r.set_xscale('log')
    ax1r.set_title("I(q) Residuals")
    ax1r.set_xlabel("q (1/$\AA$)")
    #ax1r.set_ylabel("Fractional Residuals ((I(q)_fit - I(q)_goal)/I(q)_goal")
    ax1r.set_ylabel("Fractional Residuals")
    ax1r.plot(evars.goal_q_data, (evars.initial_average_spec - evars.goal_iq_data)/evars.goal_iq_data, 'go', label='initial equal weights')
    ax1r.plot(evars.goal_q_data, (evars.average_spec - evars.goal_iq_data)/evars.goal_iq_data, 'ro', label='final optimized fit')

    ax2 = fig.add_subplot(7,2,5)
    ax2.plot(evars.s_array, 'co')
    ax2.set_xlabel("Accepted Monte Carlo Step")
    ax2.set_ylabel("S")

    ax3 = fig.add_subplot(7,2,6) #, sharex=ax2)
    ax3.plot(evars.x2_array, 'mo')
    ax3.set_xlabel("Accepted Monte Carlo Step")

    ax3s = fig.add_subplot(7,2,8) 
    ax3s.plot(evars.s_array,evars.x2_array, 'yo')
    ax3s.set_xlabel("S")

    ax4 = fig.add_subplot(7,2,7) 
    ax4.plot(evars.g_array, 'ko')
    ax4.set_xlabel("Accepted Monte Carlo Step")
    ax4.set_ylabel("G")

    #ax5 = fig.add_subplot(7,2,(9,14)) 
    ax5 = fig.add_subplot(7,2,(9,12)) 
    ax5.set_yscale('log')
    ax5.plot(evars.initial_weights, 'g-', label='initial weights')
    ax5.plot(evars.weights, 'ro', label='final optimized weights')
    ax5.set_xlabel("I(q) Profile Number")
    ax5.set_ylabel("Fractional Weights")

    create_sorted_weight_data(other_self)

    ax6 = fig.add_subplot(7,2,(13,14)) 
    ax6.plot(evars.discrepency_and_weights_array[:,0], evars.discrepency_and_weights_array[:,1],'ro')
    ax6.set_xlabel("cummulative weight $\sum_k{w_k}$")

    if mvars.reduced_x2 == 3:
        ax3.set_ylabel("R-factor")
        ax3s.set_ylabel("R-factor")
        ax6.set_ylabel("R-factor")
    else:
        ax3.set_ylabel("X2")
        ax3s.set_ylabel("X2")
        ax6.set_ylabel("X2")



    ax1.legend(loc='lower left', fancybox=True, framealpha=0.1)
    ax1r.legend(loc='upper left', fancybox=True, framealpha=0.1)
    #ax5.legend(loc='upper right', fancybox=True, framealpha=0.1)

    plt.show()

    return


def plotme_bokeh(other_self):

    print('executing server; num_procs = 4')
    server = Server({'/': modify_doc}, num_procs=4)
    ### TODO: following line does not work (see method for details)
    
    #server = Server({'/': modify_doc_new}, num_procs=4)
    server.start()
    #server.run_until_shutdown()

    print('Opening Bokeh application on http://localhost:5006/')

    server.io_loop.add_callback(server.show, "/")
    server.io_loop.start()
    #server.io_loop.run_until_shutdown()
    #server.io_loop.stop()
    #import time
    #time.sleep(20)
    #server.stop()

    return


def main(other_self):

    log = other_self.log
    log.debug('in eros_util.main')
    pgui = other_self.run_utils.print_gui
    mvars = other_self.mvars
    evars = other_self.evars

    initialize_data(other_self)

    number_of_profiles = len(evars.iq_data)

    evars.weights = numpy.ones(number_of_profiles, numpy.float32)/number_of_profiles
    sum_weights = numpy.sum(evars.weights)
    print 'sum_weights = ', sum_weights

    ## initialize equal weights with sum == 1.0
    evars.x2 = numpy.zeros(number_of_profiles, numpy.float32)

    ## calculate initial X2
    for i in xrange(number_of_profiles):
        evars.x2[i] = calculate_x2(other_self, i)

    print 'using reduced x2 option = ', mvars.reduced_x2    
    print 'min X2 = ', numpy.min(evars.x2) 
    print 'max X2 = ', numpy.max(evars.x2) 
    print 'average X2 = ', numpy.sum(evars.x2)/len(evars.x2)

    evars.average_x2 = calculate_average_x2(other_self, evars.weights)
    print 'average_x2 (equal weights) = ', evars.average_x2

    print 'minimizing G'

    minimize_g(other_self)

    print 'done minimizing G'
#    print 'evars.weights = ', evars.weights

    if MATPLOTLIB_PLOT:
        plotme_matplotlib(other_self)

    elif BOKEH_PLOT:
        plotme_bokeh(other_self)

    return 

def minimize_g(other_self):
    
    mvars = other_self.mvars
    evars = other_self.evars

    evars.average_x2 = calculate_average_x2(other_self, evars.weights)
    print 'average x2 = ', evars.average_x2

    evars.initial_weights = evars.weights

    number_of_monte_carlo_steps = mvars.number_of_monte_carlo_steps
    weight_step_size_fraction = mvars.weight_step_size_fraction
    theta = mvars.theta
    beta = mvars.beta

    print 'number of Monte Carlo steps = ', number_of_monte_carlo_steps
    print 'weight step size (fraction) = ', weight_step_size_fraction
    print 'theta = ', theta
    print 'beta = ', beta

    accepted = 0
    g_lower_accepted = 0
    boltz_failed = 0
    boltz_accepted = 0

    s_array = [] ; g_array = [] ; x2_array = []
    
    for i in xrange(number_of_monte_carlo_steps):
        #print '.', ; sys.stdout.flush()
        if DEBUG and i > 0:
            print_progress(i+1, number_of_monte_carlo_steps, prefix = 'Progress:', suffix = 'Done', bar_length = 50)
        new_weights = evars.weights + weight_step_size_fraction * (numpy.random.random_sample(len(evars.weights))-0.5)*evars.weights
        new_weights[new_weights < 0.0] = 0.0
        new_weights = new_weights / numpy.sum(new_weights)
        if i==0:
            evars.initial_weights = evars.weights

        S = - numpy.sum(new_weights * numpy.log(new_weights/evars.initial_weights))
        average_x2 = calculate_average_x2(other_self, new_weights)
        G = average_x2 - (theta * S)
        if i == 0:
            last_G = G
            print 'S = ', S, 'x2 = ', average_x2, ' G = ', G
            s_array.append(S) ; g_array.append(G) ; x2_array.append(average_x2)
        else:
            if G < last_G:
                evars.weights = new_weights
                last_G = G
           #     print 'S = ', S, 'x2 = ', average_x2, ' G = ', G
                accepted += 1
                g_lower_accepted += 1
                s_array.append(S) ; g_array.append(G) ; x2_array.append(average_x2)
            else:
                arg = - beta * (G-last_G)                    
                ran = numpy.random.random_sample() 
                if ran < math.exp(arg):
                    evars.weights = new_weights
                    last_G = G
           #         print 'boltz: S = ', S, 'x2 = ', average_x2, ' G = ', G
                    accepted += 1
                    boltz_accepted += 1
                    s_array.append(S) ; g_array.append(G) ; x2_array.append(average_x2)
                else:
                    boltz_failed += 1
    
    ## N = evars.number_q_points
    ## X2 = (1/N) SUM_{i=1}^N { ( (c * I_theory(q_i) - I_goal(q_i) )^2 ) / ( sigma^2(q_i) )  }

    ## vary cluster weights but prevent overfitting
    ##  by S = - SUM{k=1,N} { w_k ln ( w_k / w_k^(0) ) }
    ##  use free-energy function G = X2 - theta * S
    ##  Minimize G with respect to w_k by simulated annealing
    ##  to get optimal set of weights

    ## For large theta, minimizing G leads to small changes in w_k
    ## For small theta (to zero), best agreement with possible over fitting

    evars.s_array = s_array
    evars.g_array = g_array
    evars.x2_array = x2_array

    print 'accepted ', (float(accepted)/float(number_of_monte_carlo_steps))*100.0, ' percent'
    print 'g_lower_accepted ', (float(g_lower_accepted)/float(number_of_monte_carlo_steps))*100.0, ' percent'
    print 'attempted ', boltz_failed + boltz_accepted, ' boltzman moves'
    print 'accepted ', float(boltz_accepted)/(boltz_failed + boltz_accepted)*100.0, ' percent boltzman moves'

    return 

if __name__ == "__main__":

    pass
        

