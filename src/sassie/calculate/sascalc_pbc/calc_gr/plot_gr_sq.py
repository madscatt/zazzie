import numpy as np
from bokeh.plotting import figure, output_file, show
# from bokeh.charts import Line, output_file, show
from bokeh.layouts import row
output_file('argon_85K.html')

# get the data
gr = np.loadtxt('argon_85K_gr.dat')
sq = np.loadtxt('argon_85K_sq.dat')
gr_int = np.loadtxt('gr.int')
gr_ian = np.loadtxt('ian_sq_from_gr.dat')
sim_sq = np.loadtxt('../iq_0to695.dat')
# sim_sq_wr = np.loadtxt('../iq_690_wr.dat')
sim_sq_wor = np.loadtxt('../iq_690_wor.dat')

# rescale sim_sq
i690_sim = np.argmin(np.abs(sim_sq_wor[:, 0] - 2))
s690_sim = sim_sq_wor[i690_sim, 1]

i_sim = np.argmin(np.abs(sim_sq[:, 0] - 2))
s_sim = sim_sq[i_sim, 1]


i_ian = np.argmin(np.abs(gr_ian[:, 0] - 2))
s_ian = gr_ian[i_ian, 1]

sim_sq[:, 1] *= s_ian/s_sim
sim_sq_wor[:, 1] *= s_ian/s690_sim


# plot gr
pg = figure(y_axis_label='g(r)', x_axis_label=u'r (\u212B)')
pg.line(gr[:,0], gr[:,1])

# plot sq
ps = figure(y_axis_label='S(q)', x_axis_label=u'q (1/\u212B)')

ps.line(sq[:,0], sq[:,1], legend='Yarnell et al.', color='black', line_width=2)
ps.line(sim_sq[:,0], sim_sq[:,1], legend='Simulated Box, all', color='orange', line_width=2)
# ps.line(sim_sq_wr[:,0], sim_sq_wr[:,1], legend=r'Simulated Box w/ r2', color='yellow', line_width=2)
ps.line(sim_sq_wor[:,0], sim_sq_wor[:,1], legend='Simulated Box, 690', color='red', line_width=2)


# ps.line(gr_int[:,0], gr_int[:,1], legend='jec simps', color='blue', line_width=2)
# ps.line(gr_int[:,0], gr_int[:,2], legend='jec trapz', color='orange', line_width=2)
ps.line(gr_ian[:,0], gr_ian[:,1], legend='ian simps', color='violet', line_width=2)


show(row(pg, ps))