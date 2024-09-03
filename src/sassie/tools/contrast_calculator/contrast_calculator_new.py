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
import time
import math
import platform
import numpy as np
import sasmol.sasmol as sasmol
import sassie.util.sasconfig as sasconfig
import sassie.util.module_utilities as module_utilities
import sassie.tools.contrast_calculator.contrast_helper as contrast_helper
# import contrast_helper as contrast_helper
from scipy import stats

# -------------------
# for formula parser
from re import findall
# -------------------
import Gnuplot
import Gnuplot.PlotItems
import Gnuplot.funcutils  # for plotting

#    CONTRAST
#
#    11/21/2012  --  broadly adapted from center.py :      jc/ks
#    1/3/2013    --  ready for beta testing; no FASTA
#                    or solvent functionality yet          ks
#    4/1/2013    --  added plotting capability             sk
#    5/1/2013    --  works for any number of protein
#                    dna or rna components with different
#                    amounts of deuteration                sk
#    5/15/2013   --  solvent functionality added
#                    code follows that of Whitten et al.
#                    J. Appl. Cryst. (2008). 41, 222-226   sk
#    6/18/2014   --  now works for any chemical formula.
#                    Inputs are the chemical formula,
#                    mass density, number of exchangeable
#                    hydrogens and fraction of exchangeable
#                    hydrogens that exchange               sk


'''
    CONTRAST_CALCULATOR is the module that calculates SLD, I(0) and contrast v. %D2O 
    for a given molecule or complex based only on sequence or chemical formula.  It will
    read the sequence from a FASTA or PDB file.

'''

if sasconfig.__level__ == "DEBUG":
    DEBUG = True

app = 'contrast_calculator'


class module_variables():

    def __init__(self, parent=None):
        self.app = app


class contrast_calculator_input_variables():

    def __init__(self, parent=None):
        pass


class contrast_calculator():

    def __init__(self, parent=None):
        pass

    def main(self, input_variables, ivariables, solvvariables, chemvariables, txtOutput):

        self.mvars = module_variables()

        self.avars = contrast_calculator_input_variables()

        self.run_utils = module_utilities.run_utils(app, txtOutput)

        self.run_utils.setup_logging(self)

        self.log.debug('in main')

        self.unpack_variables(input_variables)

        self.unpack_ivariables(ivariables)

        self.unpack_solvvariables(solvvariables)

        self.unpack_chemvariables(chemvariables)

        self.run_utils.general_setup(self)

        self.initialization()

        self.contrast()

        self.epilogue()

        return


#   for plotting

    def wait(str=None, prompt='Plot will clear in 10 seconds ...\n'):
        if str is not None:
            print str

        try:
            if (platform.system() == "Linux"):
                import curses
                stdscr = curses.initscr()
                stdscr.addstr('press a key to continue')
                c = stdscr.getch()
                curses.endwin()
#               raw_input(prompt)
        except:
            time.sleep(1)


#    pgui performs this function
#    def print_failure(message,txtOutput):
#
#    txtOutput.put("\n\n>>>> RUN FAILURE <<<<\n")
#    txtOutput.put(">>>> RUN FAILURE <<<<\n")
#    txtOutput.put(">>>> RUN FAILURE <<<<\n\n")
#    txtOutput.put(message)
#
#    return


    def unpack_variables(self, variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_variables')

        mvars.runname = variables['runname'][0]
        mvars.inpath = variables['inpath'][0]
        mvars.outfile = variables['outfile'][0]
        mvars.solute_conc = variables['solute_conc'][0]
        mvars.d2ostep = variables['d2ostep'][0]
        mvars.fexchp = variables['fexchp'][0]
        mvars.fexchn = variables['fexchn'][0]
        mvars.plotflag = variables['plotflag'][0]  # for plotting

        log.debug(vars(mvars))

        return

    def unpack_ivariables(self, ivariables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_ivariables')

        mvars.seqfiles = []
        mvars.numunits = []
        mvars.fracdeut = []
        mvars.moltype = []
        mvars.isFasta = []

        log.debug('len(ivariables): %i' % len(ivariables))
        mvars.numfiles = len(ivariables)
#       print 'numfiles = ', mvars.numfiles

        for i in range(mvars.numfiles):
            mvars.seqfiles.append(ivariables[i][0])
            mvars.numunits.append(ivariables[i][1])
            mvars.fracdeut.append(ivariables[i][2])
            mvars.moltype.append(ivariables[i][3])
            mvars.isFasta.append(ivariables[i][4])

        log.debug(vars(mvars))

        return

    def unpack_solvvariables(self, solvvariables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_solvvariables')

        mvars.solv_comp = []  # chemical formula (see below)
        mvars.solv_conc = []  # concentration (M)

        log.debug('len(solvvariables): %i' % len(solvvariables))
        mvars.numsolv = len(solvvariables)
#       print 'numsolv = ', mvars.numsolv
#       for i in range(mvars.numsolv):
#           print solvvariables[i]

        for i in range(mvars.numsolv):
            if solvvariables[i][0] == "":
                chemical_string = ''
                mvars.numsolv = 0
            else:
                chemical_string = ''
                for key, value in solvvariables[i][0].items():
                    chemical_string += key + str(value)
#           print 'solvent_string = ',chemical_string
            mvars.solv_comp.append(chemical_string)
            mvars.solv_conc.append(solvvariables[i][1])

        log.debug(vars(mvars))

        return

    def unpack_chemvariables(self, chemical_variables):

        log = self.log
        mvars = self.mvars
        log.debug('in unpack_chemvariables')

        mvars.chem_comp = []     # chemical formula (see below)
        mvars.hexchem = []       # number of exchangeable hydrogens
        mvars.chem_density = []  # mass density
        mvars.fexchchem = []     # fraction of exchangeable hydrogens that actually do exchange

        log.debug('len(chemical_variables): %i' % len(chemical_variables))
        mvars.numformulas = len(chemical_variables)
#       print 'numformulas = ', mvars.numformulas
#       for i in range(mvars.numformulas):
#           print chemical_variables[i]

        for i in range(mvars.numformulas):
            if chemical_variables[i][0] == "":
                chemical_string = ''
                mvars.numformulas = 0
            else:
                chemical_string = ''
                for key, value in chemical_variables[i][0].items():
                    chemical_string += key + str(value)
#           print 'chemical_string = ',chemical_string
            mvars.chem_comp.append(chemical_string)
            mvars.hexchem.append(chemical_variables[i][1])
            mvars.fexchchem.append(chemical_variables[i][2])
            mvars.chem_density.append(chemical_variables[i][3])

        log.debug(vars(mvars))

        return

    def initialization(self):
        '''
        method to prepare for contrast calculation
        '''

        log = self.log
        log.debug('in initialization')
        pgui = self.run_utils.print_gui

        mvars = self.mvars
        avars = self.avars

        print 'numfiles, numsolv, numformulas :', mvars.numfiles, mvars.numsolv, mvars.numformulas

        avars.contpath = os.path.join(mvars.runname, 'contrast_calculator')
        direxist = os.path.exists(avars.contpath)
        if (direxist == 0):
            os.system('mkdir -p ' + avars.contpath)

        prefix = mvars.outfile.find('.')
        if prefix == -1:
            coreoutfile = mvars.outfile
        else:
            coreoutfile = mvars.outfile[0:prefix]

        avars.izerofile = coreoutfile + '_izero.txt'
        avars.izfile = open(os.path.join(avars.contpath, avars.izerofile), 'w')

        avars.contrastfile = coreoutfile + '_contrast.txt'
        avars.cotfile = open(os.path.join(
            avars.contpath, avars.contrastfile), 'w')

        avars.scatlendenfile = coreoutfile + '_sld.txt'
        avars.sldfile = open(os.path.join(
            avars.contpath, avars.scatlendenfile), 'w')

        log.debug(vars(mvars))
        log.debug(vars(avars))

    def contrast(self):
        '''
        CONTRAST is the function that calculates SLD, I(0) and contrast v. %D2O for a
            given molecule or complex based only on sequence or chemical formula.  It will
            read the sequence from a FASTA or PDB file.

        INPUT:

            runname:                            project name
            inpath:                             path name for input files (pdb or sequence)
            numfiles:                           number of input files (protein, rna or dna)
            solute_conc:                        concentration of solute
            d2ostep:                            step in fraction D2O used when calculating SLD, contrast and I(0)
            fexchp:                             fraction of exchangeable hydrogen atoms that exchange for protein components
            fexchn:                             fraction lf exchangeable hydrogen atoms that exchange for nucleic acid (rna, dna) components
            seqfiles:                           names of sequence or pdb files
            numunits:                           number of times the sequence from seqfile is used in the protein, rna and/or dna complex
            fracdeut:                           fraction of deuteration for the subunit represented in the seqfile
            moltype:                            type of molecule represented in the seqfile (protein, rna, or dna)
            isFasta:                            indicates whether the file is a FASTA file (0=PDB, 1=FASTA)
            plotflag:                           flag to indicate whether to plot data (0=no plot, 1=plot)
            numsolv:                            number of non-water components in the solvent
            solv_comp:                          chemical formula representing the non-water solvent component
            solv_conc:                          molar concentration of the non-water solvent component
            number_of_chemicals:                number of non-protein (rna or dna) components in the solute
            formula_array:                      chemical formulas for the solute components
            number_exchangeable_hydrogens:      the number of exchangeable hydrogen atoms in the solute component
            fraction_exchangeable_hydrogens:    fraction of exchangeable hydrogens that exchange in the solute compoent
            mass_density:                       mass density of solute component   


        OUTPUT:
            outfile:              	output filename prefix

            outfile_contrast.txt
            outfile_sld.txt
            outfile_izero.txt
            '''

    # mvars: runname, inpath, outfile, solute_conc, d2ostep, fexchp, fexchn, plotflag,
    #        numfiles, seqfiles, numunits, fracdeut, moltype, isFasta, numfiles, solv_comp, solv_conc, numsolv
    #        chem_comp, hexchem, fexchchem, chem_density, numformulas

    # avars: contpath, izfile, cotfile, sldfile, izerofile, contrastfile, scatlendenfile
    #

        log = self.log
        pgui = self.run_utils.print_gui
        log.debug('in contrast')

        mvars = self.mvars
        avars = self.avars

        # ttxt=time.ctime()
        ttxt = time.asctime(time.gmtime(time.time()))
        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n" % (st))
        pgui("DATA FROM RUN: %s \n\n" % (ttxt))

        na = 6.023e23

        vdna = 0.56    # using an average DNA value for both DNA and RNA, since the RNA value is only slightly smaller
        vprot = 0.730

# fraction of exchangeable H that are assumed to exchange
# now obtained from graphical contrast screen as a variable (9/2013)
        fexchprot = mvars.fexchp
        fexchdna = mvars.fexchn

#       print 'fraction of exchangeables'
#       print mvars.fexchp, mvars.fexchn

        working_conc = mvars.solute_conc / 1000  # convert mg/ml to g/ml

# test if this works once program is running in 2.0 format
# if working, use sasproperties dictionaries for all dictionaries
# except atomic_dict since we are still adding elements to it
# and the sasproperties dictionary is not up to date

#       atomic = sasmol.sasproperties.Atomic()
#       protein_dict = atomic.protein_sl()
#       dna_dict = atomic.dna_sl()
#       rna_dict = atomic.rna_sl()
        protein_dict = contrast_helper.protein_sl()
        dna_dict = contrast_helper.dna_sl()
        rna_dict = contrast_helper.rna_sl()
        atomic_dict = contrast_helper.element_sl()

# water properties used to calculate bh2o,bd2o,bxh2o,rh2o,rd2o,rxh2o

        hstats = atomic_dict['H']
#        print hstats
        dstats = atomic_dict['D']
#        print dstats
        ostats = atomic_dict['O']
#        print ostats

        bh = hstats[3]  # in units of 10^-12 cm
        bd = dstats[3]
        bo = ostats[3]
        bxh = hstats[2]
        bxo = ostats[2]

        bh2o = bh * 2 + bo
        bd2o = bd * 2 + bo
        bxh2o = bxh * 2 + bxo
#       print bh, bd, bxh, bh2o, bd2o, bxh2o

        h2o_conc = 1.00e3 / (hstats[0] * 2 + ostats[0])  # molar conc of water
        d2o_conc = 1.00e3 / (dstats[0] * 2 + ostats[0])
        h2o_vol = 1.0e-3 / na / h2o_conc * 1e30          # Ang^3 (10^-24 cm^3)
        d2o_vol = 1.0e-3 / na / d2o_conc * 1e30

#       print h2o_conc, d2o_conc, h2o_vol, d2o_vol


# determine if there are additional components in solvent

        if mvars.numsolv == 0:
            usesolventinfo = 0
        else:
            usesolventinfo = 1

        elementmatrix = []

#        print 'numsolv = ', mvars.numsolv

        # What if nothing was entered initially? This takes care of that possibility.
        if usesolventinfo == 1:
            if mvars.solv_conc[0] == "":
                mvars.numsolv = 0
                usesolventinfo = 0
#        print 'numsolv = ', mvars.numsolv

# parse the solvent composition to get the total number of each element
# formula can't have parenthesis at this stage

        for m in range(mvars.numsolv):
            #            print 'm = ', m
            thissolv_comp = mvars.solv_comp[m]

#            print thissolv_comp
#            print 'length = ', len(thissolv_comp)

            elements = []
            elementarray = []
            numberarray = []

            elements = [k for k in findall(
                r'([A-Z][a-z]*)(\d*)', thissolv_comp) if k]
#            print elements
#            print len(elements)

            for i in range(len(elements)):
                #                print len(elements[i])
                for j in range(len(elements[i])):
                    #                print i,j,elements[i][j]
                    if j == 0:
                        elementarray.append(elements[i][j])
                    if j == 1:
                        if elements[i][j] == "":
                            numberarray.append('1')
                        else:
                            numberarray.append(elements[i][j])
#            print elementarray
#            print numberarray
            elementmatrix.append([elementarray, numberarray])

#        print elementmatrix

# arrays for solvent components

        volsolvarray = []
        bxsolvarray = []
        bhsolvarray = []
        mwsolvarray = []
        concsolvarray = []
        rxsolvarray = []
        rhsolvarray = []

# calculate solvent properties if solvent contains only water

        if usesolventinfo == 0:

            rxh2o = bxh2o * 100 / h2o_vol   # in units of 10^10 cm^-2
            rh2o = bh2o * 100 / h2o_vol     # H2O volume is used for both H2O and D2O to get
            rd2o = bd2o * 100 / h2o_vol     # the accepted value of 6.4 for rd2o
#            rd2o2 = bd2o*100/d2o_vol
#           print rh2o, rd2o, rd2o2, rxh2o

            rsx = rxh2o
            rsd2o = rd2o
            rsh2o = rh2o
#            print rsh2o, rsd2o, rsx

# calculate solute properties if usesolventinfo=1
# H-D exchange between water and small solute molecules in the solvent is
# not considered

        elif usesolventinfo == 1:

            rxh2o = bxh2o * 100 / h2o_vol   # in units of 10^10 cm^-2
            rh2o = bh2o * 100 / h2o_vol     # for printing out to files

            soluteconc = 0.0

            for i in range(mvars.numsolv):

                mwsolv = 0.0
                bhsolv = 0.0
                bxsolv = 0.0
                volsolv = 0.0

                for j in range(len(elementmatrix[i][0])):

                    element = elementmatrix[i][0][j]
                    number = int(elementmatrix[i][1][j])
#                    print 'element, number: ', element,number

                    estats = atomic_dict[element]

                    mwsolv = mwsolv + estats[0] * number
                    volsolv = volsolv + estats[1] * number
                    bhsolv = bhsolv + estats[3] * number
                    bxsolv = bxsolv + estats[2] * number
#                    print mwsolv,volsolv,bhsolv,bxsolv

                concsolvarray.append(float(mvars.solv_conc[i]))
                bxsolvarray.append(bxsolv)
                bhsolvarray.append(bhsolv)
                mwsolvarray.append(mwsolv)
                volsolvarray.append(volsolv)

                rxsolvarray.append(
                    bxsolvarray[i] * 100 / volsolvarray[i])                             # 10^10 cm^-2
                rhsolvarray.append(bhsolvarray[i] * 100 / volsolvarray[i])

                # total conc of solutes
                soluteconc = soluteconc + \
                    (concsolvarray[i] * volsolvarray[i] * na)

#            print concsolvarray
#            print bxsolvarray
#            print bhsolvarray
#            print mwsolvarray
#            print volsolvarray
#            print rxsolvarray
#            print rhsolvarray
#            print 'solute conc: ', soluteconc

# adjust water concentration due to presence of solutes

            # 1e27 is the number of A^3 per liter
            newh2o_conc = h2o_conc * (1e27 - soluteconc) / 1e27

#            print 'new h2o conc: ', newh2o_conc

            bxsolvtot = bxh2o * newh2o_conc
            bhsolvtot = bh2o * newh2o_conc
            bdsolvtot = bd2o * newh2o_conc
            volsolvtot = h2o_vol * newh2o_conc

#            print bhsolvtot, bdsolvtot, bxsolvtot, volsolvtot

            for i in range(mvars.numsolv):

                bxsolvtot = bxsolvtot + concsolvarray[i] * bxsolvarray[i]
                bhsolvtot = bhsolvtot + concsolvarray[i] * bhsolvarray[i]
                bdsolvtot = bdsolvtot + concsolvarray[i] * bhsolvarray[i]
                volsolvtot = volsolvtot + concsolvarray[i] * volsolvarray[i]

#            print bh, bd, bhsolvtot, bdsolvtot, bxsolvtot, volsolvtot

            rhsolvtot = bhsolvtot * 100 / volsolvtot                # 10^10 cm^-2
            rdsolvtot = bdsolvtot * 100 / volsolvtot
            rxsolvtot = bxsolvtot * 100 / volsolvtot

# print rhsolvtot, rdsolvtot, rxsolvtot

            rsx = rxsolvtot
            rsd2o = rdsolvtot
            rsh2o = rhsolvtot

#            print rsh2o, rsd2o, rsx


# Get the properties of the COMPLEX.

# First, determine if there are any chemical forumlas, which indicates that there are components
# of the complex that are not protein, dna or rna.

        formulamatrix = []
        allformulas = []
        vxchemarray = []
        mwchemarray = []
        hchemarray = []
        dchemarray = []
        hexchemarray = []
        fexchchemarray = []
        denschemarray = []
        btotchemarray = []
        bxtotchemarray = []

#        print 'numformulas =', mvars.numformulas


# Each chemical formula is treated as a separate component regardless of the amount of deuteration
# because the mass densities could be different

        if mvars.numformulas > 0:
            #            print 'CHEMICAL FORMULA'

            # parse the chemical formula to get the total number of each element

            for i in range(mvars.numformulas):
                thischem_comp = mvars.chem_comp[i]
                allformulas.append(thischem_comp)

#                print thischem_comp
#                print 'length = ', len(thischem_comp)

                elements = []
                elementarray = []
                numberarray = []

                elements = [k for k in findall(
                    r'([A-Z][a-z]*)(\d*)', thischem_comp) if k]
#                print elements
#                print len(elements)

                for i in range(len(elements)):
                    #                    print len(elements[i])
                    for j in range(len(elements[i])):
                        #                        print i,j,elements[i][j]
                        if j == 0:
                            elementarray.append(elements[i][j])
                        if j == 1:
                            if elements[i][j] == "":
                                numberarray.append('1')
                            else:
                                numberarray.append(elements[i][j])
#                print elementarray
#                print numberarray
                formulamatrix.append([elementarray, numberarray])

# print formulamatrix

# calculate Mw and total x-ray and neutron scattering lengths for each formula
# the calculations assume D is explicitely listed in the chemical formula

            for i in range(mvars.numformulas):

                mwchem = 0.0
                hchem = 0.0
                bh2chem = 0.0
                bxchem = 0.0
                dchem = 0.0

                for j in range(len(formulamatrix[i][0])):

                    elementname = formulamatrix[i][0][j]
                    elementnumber = int(formulamatrix[i][1][j])
#                    print 'element, number: ', elementname,elementnumber

                    estats = atomic_dict[elementname]

                    mwchem = mwchem + estats[0] * elementnumber
                    bh2chem = bh2chem + estats[3] * elementnumber
                    bxchem = bxchem + estats[2] * elementnumber
                    if elementname == 'H':
                        hchem = hchem + elementnumber
                    if elementname == 'D':
                        dchem = dchem + elementnumber
#                print mwchem,bh2chem,bxchem,hchem,dchem

                mwchemarray.append(mwchem)
                hchemarray.append(hchem)
                dchemarray.append(dchem)
                hexchemarray.append(float(mvars.hexchem[i]))
                fexchchemarray.append(float(mvars.fexchchem[i]))
                denschemarray.append(float(mvars.chem_density[i]))
                btotchemarray.append(bh2chem)
                bxtotchemarray.append(bxchem)

            for i in range(len(mwchemarray)):
                vol = mwchemarray[i] / (denschemarray[i] * na)
                vxchemarray.append(vol)

#            print mwchemarray
#            print hchemarray
#            print dchemarray
#            print hexchemarray
#            print fexchchemarray
#            print denschemarray
#            print btotchemarray
#            print bxtotchemarray
#            print vxchemarray


# Get the properties of protein, dna and/or rna in the complex.
# First, determine how many different molecule types there are and the fract deut for each given molecule type.
# Each different fract deut is one component in the total complex.

        proteinmatrix = []
        dnamatrix = []
        rnamatrix = []
        mwprotarray = []
        hprotarray = []
        hexprotarray = []
        btotprotarray = []
        bxtotprotarray = []
        mwdnaarray = []
        hdnaarray = []
        hexdnaarray = []
        btotdnaarray = []
        bxtotdnaarray = []
        mwrnaarray = []
        hrnaarray = []
        hexrnaarray = []
        btotrnaarray = []
        bxtotrnaarray = []

#        print 'numfiles =', mvars.numfiles

        if mvars.numfiles > 0:
            # In case 0 files was entered initially but "Then Click Here" wasn't clicked
            if mvars.seqfiles[0] == "":
                mvars.numfiles = 0

#       print 'numfiles =', mvars.numfiles

        for index in range(mvars.numfiles):

            thismoltype = mvars.moltype[index]
            thisnumunits = int(mvars.numunits[index])
            thisfracdeut = float(mvars.fracdeut[index])
            thisIsFasta = mvars.isFasta[index]
            thisfilename = mvars.seqfiles[index]
#            print thisfilename,thisfracdeut,thisnumunits,thisIsFasta,thismoltype,index

            if thismoltype == "protein":
                proteinmatrix.append(
                    [thisfilename, thisfracdeut, thisnumunits, thisIsFasta])
            elif thismoltype == "dna":
                dnamatrix.append([thisfilename, thisfracdeut,
                                  thisnumunits, thisIsFasta])
            elif thismoltype == "rna":
                rnamatrix.append([thisfilename, thisfracdeut,
                                  thisnumunits, thisIsFasta])
#            else:
#                print("Molecule type is neither protein nor dna nor rna.")      # this issue should be caught in the contrast_calculator_filter

# First calculate the contrast parameters for the protein components.
# Sort the protein values by fract deut and bin accordingly.  Then go
# through each bin and calculate the contrast parameters and store them in
# arrays.  The size of the arrays will depend on the number of components.

        proteinfdval = []
        dnafdval = []
        rnafdval = []
        allfiles = []
        allnumunits = []

        if len(proteinmatrix) > 0:
            #            print 'PROTEIN'
            #            print proteinmatrix
            newproteinmatrix = sorted(proteinmatrix, key=lambda x: x[
                1])                                     # sort by fract deut
#            print newproteinmatrix

            current = newproteinmatrix[0][1]
            proteinfdval.append(newproteinmatrix[0][1])
            numvalues = 1

            # determine how many components there are
            for i in range(1, len(newproteinmatrix)):
                if newproteinmatrix[i][1] != current:
                    numvalues = numvalues + 1
                    proteinfdval.append(newproteinmatrix[i][1])
                    current = newproteinmatrix[i][1]
#            print 'numvalues = ', numvalues                                #number of components
#            print 'proteinfdval = ', proteinfdval                          #value of fract deut

            proteinbins = [[[] for x in range(len(newproteinmatrix))]
                           for y in range(numvalues)]

            # sort into bins according to fract deut
            for i in range(numvalues):
                j = 0
                k = 0
                while j < len(newproteinmatrix):
                    if newproteinmatrix[j][1] == proteinfdval[i]:
                        proteinbins[i][k] = newproteinmatrix[j]
                        j = j + 1
                        k = k + 1
                    else:
                        j = j + 1
                        k = 0
#            print proteinbins

# calculate Mw, number of hydrogens, number of exchangeable hydrogens and
# total x-ray and neutron scattering lengths for each protein component

            for i in range(numvalues):
                mwprottot = 0.0
                hprottot = 0.0
                hexprottot = 0.0
                btotprot = 0.0
                bxtotprot = 0.0

                # eliminate empty spaces in bins
                proteinbins[i] = [e for e in proteinbins[i] if e]
#                print proteinbins[i]
#                print len(proteinbins[i])

                for j in range(len(proteinbins[i])):
                    # read the file
                    filename = proteinbins[i][j][0]
                    nunits = proteinbins[i][j][2]
                    # for printing to output files
                    allfiles.append(filename)
                    allnumunits.append(nunits)
#                    print 'protein filename: ', filename
                    thisfile = os.path.join(mvars.inpath, filename)

                    # not thisIsFasta
                    if proteinbins[i][j][3] != '1':
                        m1 = sasmol.SasMol(0)
                        m1.read_pdb(thisfile)

                        resids = m1.resid()
                        resnames = m1.resname()

                        seq = contrast_helper.sequence(resids, resnames)
                    else:
                        seq = contrast_helper.FASTA_sequence(thisfile)

                    mwprot = 0.0
                    hprot = 0.0
                    hexprot = 0.0
                    bh2oprot = 0.0
                    bxprot = 0.0

                    for m in protein_dict.keys():
                        nres = seq.count(m)
                        aastats = protein_dict[m]
                        mwprot = mwprot + \
                            proteinbins[i][j][2] * nres * \
                            aastats[0]                                      # thisnumunits
                        hprot = hprot + \
                            proteinbins[i][j][2] * nres * aastats[6]
                        hexprot = hexprot + \
                            proteinbins[i][j][2] * nres * aastats[5]
                        bh2oprot = bh2oprot + \
                            proteinbins[i][j][2] * nres * aastats[3]
                        bxprot = bxprot + \
                            proteinbins[i][j][2] * nres * aastats[2]

                    mwprottot = mwprottot + mwprot + \
                        proteinbins[i][j][1] * (hprot - hexprot)
                    hprottot = hprottot + hprot
                    hexprottot = hexprottot + hexprot
                    btotprot = btotprot + bh2oprot + \
                        (hprot - hexprot) * \
                        proteinbins[i][j][1] * \
                        (bd - bh)                    # thisfracdeut
                    bxtotprot = bxtotprot + bxprot

                mwprotarray.append(mwprottot)
                hprotarray.append(hprottot)
                hexprotarray.append(hexprottot)
                btotprotarray.append(btotprot)
                bxtotprotarray.append(bxtotprot)

#            print mwprotarray
#            print hexprotarray
#            print btotprotarray
#            print bxtotprotarray

# Now repeat for the DNA contrast parameters

        if len(dnamatrix) > 0:
            #            print 'DNA'
            #            print dnamatrix
            newdnamatrix = sorted(dnamatrix, key=lambda x: x[
                1])                                           # sort by fract deut
#            print newdnamatrix

            current = newdnamatrix[0][1]
            dnafdval.append(newdnamatrix[0][1])
            numdnavalues = 1

            # determine how many components there are
            for i in range(1, len(newdnamatrix)):
                if newdnamatrix[i][1] != current:
                    numdnavalues = numdnavalues + 1
                    dnafdval.append(newdnamatrix[i][1])
                    current = newdnamatrix[i][1]
#            print 'numdnavalues = ', numdnavalues                          #number of components
#            print 'dnafdval = ', dnafdval                                  #value of fract deut

            dnabins = [[[] for x in range(len(newdnamatrix))]
                       for y in range(numdnavalues)]

            # sort into bins according to fract deut
            for i in range(numdnavalues):
                j = 0
                k = 0
                while j < len(newdnamatrix):
                    if newdnamatrix[j][1] == dnafdval[i]:
                        dnabins[i][k] = newdnamatrix[j]
                        j = j + 1
                        k = k + 1
                    else:
                        j = j + 1
                        k = 0
#            print dnabins

# calculate Mw, number of hydrogens, number of exchangeable hydrogens and
# total x-ray and neutron scattering lengths for each DNA component

            for i in range(numdnavalues):
                mwdnatot = 0.0
                hdnatot = 0.0
                hexdnatot = 0.0
                btotdna = 0.0
                bxtotdna = 0.0

                # eliminate empty spaces in bins
                dnabins[i] = [e for e in dnabins[i] if e]
#                print dnabins[i]
#                print len(dnabins[i])
                for j in range(len(dnabins[i])):
                    # read the file
                    filename = dnabins[i][j][0]
                    nunits = dnabins[i][j][2]
                    allfiles.append(filename)
                    allnumunits.append(nunits)
#                    print 'DNA filename', filename
                    thisfile = os.path.join(mvars.inpath, filename)

                    if dnabins[i][j][3] != '1':                             # not thisIsFasta
                        m1 = sasmol.SasMol(0)
                        m1.read_pdb(thisfile)

                        resids = m1.resid()
                        resnames = m1.resname()

                        seq = contrast_helper.sequence(resids, resnames)
                    else:
                        seq = contrast_helper.FASTA_sequence(thisfile)

                    mwdna = 0.0
                    hdna = 0.0
                    hexdna = 0.0
                    bh2odna = 0.0
                    bxdna = 0.0

                    for m in dna_dict.keys():
                        nres = seq.count(m)
                        aastats = dna_dict[m]
                        mwdna = mwdna + dnabins[i][j][2] * \
                            nres * \
                            aastats[0]                               # thisnumunits
                        hdna = hdna + dnabins[i][j][2] * nres * aastats[6]
                        hexdna = hexdna + dnabins[i][j][2] * nres * aastats[5]
                        bh2odna = bh2odna + \
                            dnabins[i][j][2] * nres * aastats[3]
                        bxdna = bxdna + dnabins[i][j][2] * nres * aastats[2]

                    mwdnatot = mwdnatot + mwdna + \
                        dnabins[i][j][1] * (hdna - hexdna)
                    hdnatot = hdnatot + hdna
                    hexdnatot = hexdnatot + hexdna
                    btotdna = btotdna + bh2odna + \
                        (hdna - hexdna) * \
                        dnabins[i][j][1] * \
                        (bd - bh)                        # thisfracdeut
                    bxtotdna = bxtotdna + bxdna

                mwdnaarray.append(mwdnatot)
                hdnaarray.append(hdnatot)
                hexdnaarray.append(hexdnatot)
                btotdnaarray.append(btotdna)
                bxtotdnaarray.append(bxtotdna)

#            print mwdnaarray
#            print hdnaarray
#            print hexdnaarray
#            print btotdnaarray
#            print bxtotdnaarray

# Now repeat for RNA contrast parameters

        if len(rnamatrix) > 0:
            #            print 'RNA'
            #            print rnamatrix
            newrnamatrix = sorted(rnamatrix, key=lambda x: x[
                1])                                       # sort by fract deut
#            print newrnamatrix

            current = newrnamatrix[0][1]
            rnafdval.append(newrnamatrix[0][1])
            numrnavalues = 1

            # determine how many components there are
            for i in range(1, len(newrnamatrix)):
                if newrnamatrix[i][1] != current:
                    numrnavalues = numrnavalues + 1
                    rnafdval.append(newrnamatrix[i][1])
                    current = newrnamatrix[i][1]
#            print 'numrnavalues = ', numrnavalues                          #number of components
#            print 'rnafdval = ', rnafdval                                  #value of fract deut

            rnabins = [[[] for x in range(len(newrnamatrix))]
                       for y in range(numrnavalues)]

            # sort into bins according to fract deut
            for i in range(numrnavalues):
                j = 0
                k = 0
                while j < len(newrnamatrix):
                    if newrnamatrix[j][1] == rnafdval[i]:
                        rnabins[i][k] = newrnamatrix[j]
                        j = j + 1
                        k = k + 1
                    else:
                        j = j + 1
                        k = 0
#            print rnabins

# calculate Mw, number of hydrogens, number of exchangeable hydrogens and
# total x-ray and neutron scattering lengths for each RNA component

            for i in range(numrnavalues):
                mwrnatot = 0.0
                hrnatot = 0.0
                hexrnatot = 0.0
                btotrna = 0.0
                bxtotrna = 0.0

                # eliminate empty spaces in bins
                rnabins[i] = [e for e in rnabins[i] if e]
#                print rnabins[i]
#                print len(rnabins[i])
                for j in range(len(rnabins[i])):
                    # read the file
                    filename = rnabins[i][j][0]
                    nunits = rnabins[i][j][2]
                    allfiles.append(filename)
                    allnumunits.append(nunits)
#                    print 'RNA filename', filename
                    thisfile = os.path.join(mvars.inpath, filename)

                    if rnabins[i][j][3] != '1':                             # not thisIsFasta
                        m1 = sasmol.SasMol(0)
                        m1.read_pdb(thisfile)

                        resids = m1.resid()
                        resnames = m1.resname()

                        seq = contrast_helper.sequence(resids, resnames)
                    else:
                        seq = contrast_helper.FASTA_sequence(thisfile)

                    mwrna = 0.0
                    hrna = 0.0
                    hexrna = 0.0
                    bh2orna = 0.0
                    bxrna = 0.0

                    for m in rna_dict.keys():
                        nres = seq.count(m)
                        aastats = rna_dict[m]
                        mwrna = mwrna + rnabins[i][j][2] * \
                            nres * \
                            aastats[0]                               # thisnumunits
                        hrna = hrna + rnabins[i][j][2] * nres * aastats[6]
                        hexrna = hexrna + rnabins[i][j][2] * nres * aastats[5]
                        bh2orna = bh2orna + \
                            rnabins[i][j][2] * nres * aastats[3]
                        bxrna = bxrna + rnabins[i][j][2] * nres * aastats[2]

                    mwrnatot = mwrnatot + mwrna + \
                        rnabins[i][j][1] * (hrna - hexrna)
                    hrnatot = hrnatot + hrna
                    hexrnatot = hexrnatot + hexrna
                    btotrna = btotrna + bh2orna + \
                        (hrna - hexrna) * \
                        rnabins[i][j][1] * \
                        (bd - bh)                        # thisfracdeut
                    bxtotrna = bxtotrna + bxrna

                mwrnaarray.append(mwrnatot)
                hrnaarray.append(hrnatot)
                hexrnaarray.append(hexrnatot)
                btotrnaarray.append(btotrna)
                bxtotrnaarray.append(bxtotrna)

#            print mwrnaarray
#            print hrnaarray
#            print hexrnaarray
#            print btotrnaarray
#            print bxtotrnaarray


# print date and list of filenames and chemical formulas used to output files

        avars.izfile.write('#Date: ' + ttxt + '\n\n')
        avars.cotfile.write('#Date: ' + ttxt + '\n\n')
        avars.sldfile.write('#Date: ' + ttxt + '\n\n')

        for i in range(len(allfiles)):
            if i == 0:
                files = str(allfiles[i]) + ' (' + str(allnumunits[i]) + ')'
            else:
                files = files + ', ' + \
                    str(allfiles[i]) + ' (' + str(allnumunits[i]) + ')'
        if len(allfiles) > 0:
            #            print 'Files used: ', files
            avars.izfile.write('#Files used: ' + files + '\n\n')
            avars.cotfile.write('#Files used: ' + files + '\n\n')
            avars.sldfile.write('#Files used: ' + files + '\n\n')

        for i in range(len(allformulas)):
            if i == 0:
                formulas = str(i + 1) + '. ' + str(allformulas[i])
            else:
                formulas = formulas + ', ' + \
                    str(i + 1) + '. ' + str(allformulas[i])
        if len(allformulas) > 0:
            #            print 'Chemical formulas used: ', formulas
            avars.izfile.write('#Chemical formulas used: ' + formulas + '\n\n')
            avars.cotfile.write(
                '#Chemical formulas used: ' + formulas + '\n\n')
            avars.sldfile.write(
                '#Chemical formulas used: ' + formulas + '\n\n')


# print description of the solvent components and their parameters
# Solvent component, molar conc, volume, mass, bhtot, bxtot, rxtot, rhtot

        avars.izfile.write('#Solvent Components:\n')
        avars.cotfile.write('#Solvent Components:\n')
        avars.sldfile.write('#Solvent Components:\n')

        avars.izfile.write(
            '#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
        avars.cotfile.write(
            '#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
        avars.sldfile.write(
            '#Component, molar conc, volume (A^3), Mw (kDA), x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')

# print water values

        solvstring = 'H2O'
        concstring = '%6.2f' % (h2o_conc)
        volstring = '%5.1f' % (h2o_vol)
        massstring = ' 0.018'
        bxstring = '%7.3f' % (bxh2o)
        bhstring = '%7.3f' % (bh2o)
        rxstring = '%7.3f' % (rxh2o)
        rhstring = '%7.3f' % (rh2o)

        avars.izfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' + massstring +
                           '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')
        avars.cotfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' +
                            massstring + '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')
        avars.sldfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' +
                            massstring + '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')


# print solute values

        if usesolventinfo == 1:

            for i in range(mvars.numsolv):
                solvstring = mvars.solv_comp[i]
                concstring = '%6.2f' % (concsolvarray[i])
                volstring = '%5.1f' % (volsolvarray[i])
                massstring = '%6.3f' % (mwsolvarray[i] / 1000.0)
                bxstring = '%7.3f' % (bxsolvarray[i])
                bhstring = '%7.3f' % (bhsolvarray[i])
                rxstring = '%7.3f' % (rxsolvarray[i])
                rhstring = '%7.3f' % (rhsolvarray[i])

                avars.izfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' + massstring +
                                   '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')
                avars.cotfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' +
                                    massstring + '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')
                avars.sldfile.write('#\t' + solvstring + '\t' + concstring + '\t' + volstring + '\t' +
                                    massstring + '\t' + bxstring + '\t' + bhstring + '\t' + rxstring + '\t' + rhstring + '\n')

# print total SL and SLD values if there are non-water solvent components

            avars.izfile.write(
                '#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
            avars.cotfile.write(
                '#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')
            avars.sldfile.write(
                '#Totals: x-ray SL, neutron SL (10^-12 cm), x-ray SLD, neutron SLD (10^10 cm^-2)\n')

            bxstring = '%7.3f' % (bxsolvtot)
            bhstring = '%7.3f' % (bhsolvtot)
            rxstring = '%7.3f' % (rsx)
            rhstring = '%7.3f' % (rsh2o)

            avars.izfile.write('#\t' + bxstring + '\t' + bhstring +
                               '\t' + rxstring + '\t' + rhstring + '\n')
            avars.cotfile.write('#\t' + bxstring + '\t' + bhstring +
                                '\t' + rxstring + '\t' + rhstring + '\n')
            avars.sldfile.write('#\t' + bxstring + '\t' + bhstring +
                                '\t' + rxstring + '\t' + rhstring + '\n')


# calculate x-ray properties of the system:  combine all protein
# components, all NA (DNA+RNA) and molecule components to find the
# fraction of proteins, nucleic acids and other molecules in the complex.

        avars.izfile.write('\n#Complex concentration:  ' +
                           str(mvars.solute_conc) + ' mg/ml\n')

# Mw and x-ray scattering length

#        print 'X-RAY'

        mwxprot = 0
        bxprot = 0
        mwxdna = 0
        bxdna = 0
        mwxrna = 0
        bxrna = 0
        mwxchem = 0

        fracx = []
        volx = []
        pvolx = []
        bx = []
        rx = []
        rxbar = []
        chemxsldtitle = []
        chemxcontrasttitle = []
        chemxmwtitle = []

        for i in range(len(mwprotarray)):
            mwxprot = mwxprot + mwprotarray[i] / 1.0e3
            bxprot = bxprot + bxtotprotarray[i]
#        print mwxprot, bxprot

        for i in range(len(mwdnaarray)):
            mwxdna = mwxdna + mwdnaarray[i] / 1.0e3
            bxdna = bxdna + bxtotdnaarray[i]
#        print mwxdna, bxdna

        for i in range(len(mwrnaarray)):
            mwxrna = mwxrna + mwrnaarray[i] / 1.0e3
            bxrna = bxrna + bxtotrnaarray[i]
#        print mwxrna, bxrna

        bxtotprot = bxprot
        bxtotdna = bxdna + bxrna
        bx.append(bxtotprot)
        bx.append(bxtotdna)
        vxtotdna = vdna * (mwxdna + mwxrna) * 1.0e3 / na
        vxtotprot = vprot * mwxprot * 1.0e3 / na
        volx.append(vxtotprot)
        volx.append(vxtotdna)
        pvolx.append(vprot)
        pvolx.append(vdna)

#        print 'mwchemarray length = ', len(mwchemarray)

        for i in range(len(mwchemarray)):
            mwxchem = mwxchem + mwchemarray[i] / 1.0e3
#        print mwxprot, mwxdna, mwxrna, mwxchem

        mwxcomp = mwxprot + mwxdna + mwxrna + mwxchem
#        print mwxcomp

        fxprot = mwxprot / mwxcomp
        fxdna = (mwxdna + mwxrna) / mwxcomp
        fracx.append(fxprot)
        fracx.append(fxdna)

        for i in range(len(mwchemarray)):
            fraction = mwchemarray[i] / (1.0e3 * mwxcomp)
            volume = vxchemarray[i]
            partial_vol = 1.0 / denschemarray[i]
            fracx.append(fraction)
            volx.append(volume)
            pvolx.append(partial_vol)
            bx.append(bxtotchemarray[i])
            chemxsldtitle.append('Molecule ' + str(i + 1))
            chemxcontrasttitle.append('Molecule ' + str(i + 1))
            chemxmwtitle.append('Molecule ' + str(i + 1) + ' Mw')

#        print fracx
#        print volx, pvolx
#        print bx

# x-ray SLD, contrast and I(0)

        rxdna = 0.0
        rxdnabar = 0.0
        rxprot = 0.0
        rxprotbar = 0.0
        izeq = 0.0
        contrx = 0.0
        sldx = 0.0

        bxtotprot = bxtotprot * 1.0e-12
        bxtotdna = bxtotdna * 1.0e-12

        if (mwxdna + mwxrna) > 0:
            rxdna = (bxtotdna / vxtotdna) / 1.0e10
            rxdnabar = (rxdna - rsx)
#        print rxdna,rxdnabar

        if mwxprot > 0:
            rxprot = (bxtotprot / vxtotprot) / 1.0e10
            rxprotbar = (rxprot - rsx)
#        print rxprot,rxprotbar

        rx.append(rxprot)
        rx.append(rxdna)
        rxbar.append(rxprotbar)
        rxbar.append(rxdnabar)

        if mwxchem > 0:
            for i in range(2, len(mwchemarray) + 2):
                b = bx[i] * 1.0e-12
                v = volx[i]
                rhox = (b / v) / 1.0e10
                rx.append(rhox)
                rhobarx = (rhox - rsx)
                rxbar.append(rhobarx)

#        print rx, rsx
#        print rxbar

        for i in range(len(mwchemarray) + 2):
            izeq = izeq + fracx[i] * pvolx[i] * rxbar[i] * 1.0e10
            sldx = sldx + fracx[i] * rx[i]
            contrx = contrx + fracx[i] * rxbar[i]
#        print sldx, contrx

        izerox = working_conc * (mwxcomp) * 1.0e3 / na * (izeq)**2
#        print izerox


#    write information to files

        protxsldtitle = 'Protein'
        dnaxsldtitle = 'NA'
        protxcontrasttitle = 'Protein'
        dnaxcontrasttitle = 'NA'
        protxmwtitle = 'Protein Mw'
        dnaxmwtitle = 'NA Mw'

        if mwxprot > 0:
            xsldtitle = protxsldtitle
            xcontrasttitle = protxcontrasttitle
            xmwtitle = protxmwtitle
            if (mwxdna + mwxrna) > 0:
                xsldtitle = xsldtitle + ', ' + dnaxsldtitle
                xcontrasttitle = xcontrasttitle + ', ' + dnaxcontrasttitle
                xmwtitle = xmwtitle + ', ' + dnaxmwtitle
            if mwxchem > 0:
                for i in range(len(mwchemarray)):
                    xsldtitle = xsldtitle + ', ' + chemxsldtitle[i]
                    xcontrasttitle = xcontrasttitle + \
                        ', ' + chemxcontrasttitle[i]
                    xmwtitle = xmwtitle + ', ' + chemxmwtitle[i]
        elif (mwxdna + mwxrna) > 0:
            xsldtitle = dnaxsldtitle
            xcontrasttitle = dnaxcontrasttitle
            xmwtitle = dnaxmwtitle
            if mwxchem > 0:
                for i in range(len(mwchemarray)):
                    xldtitle = xsldtitle + ', ' + chemxsldtitle[i]
                    xcontrasttitle = xcontrasttitle + \
                        ', ' + chemxcontrasttitle[i]
                    xmwtitle = xmwtitle + ', ' + chemxmwtitle[i]
        elif mwxchem > 0:
            xsldtitle = chemxsldtitle[0]
            xcontrasttitle = chemxcontrasttitle[0]
            xmwtitle = chemxmwtitle[0]
            for i in range(1, len(mwchemarray)):
                xsldtitle = xsldtitle + ', ' + chemxsldtitle[i]
                xcontrasttitle = xcontrasttitle + ', ' + chemxcontrasttitle[i]
                xmwtitle = xmwtitle + ', ' + chemxmwtitle[i]

#        print xsldtitle
#        print xcontrasttitle
#        print xmwtitle

        avars.izfile.write('\n#XRAY I(0):\n')
        avars.sldfile.write('\n#XRAY SLD (10^10 cm^-2):\n')
        avars.cotfile.write('\n#XRAY Contrast (10^10 cm^-2):\n')

        avars.izfile.write(
            '#' + xmwtitle + ', Complex Mw (kDa), I(0) (cm^-1)\n')
        avars.sldfile.write('#' + xsldtitle + ', Complex, Solvent\n')
        avars.cotfile.write('#' + xcontrasttitle + ', Complex\n')

        xsldtable = []
        xcontrasttable = []
        xmwtable = []

        mwxna = mwxdna + mwxrna

        if mwxprot > 0:
            xsldtable.append(rx[0])
            xcontrasttable.append(rxbar[0])
            xmwtable.append(mwxprot)

            if mwxna > 0:
                xsldtable.append(rx[1])
                xcontrasttable.append(rxbar[1])
                xmwtable.append(mwxna)
            if mwxchem > 0:
                for j in range(2, len(mwchemarray) + 2):
                    xsldtable.append(rx[j])
                    xcontrasttable.append(rxbar[j])
                    xmwtable.append(mwchemarray[j - 2] / 1.0e3)
        elif mwxna > 0:
            xsldtable.append(rx[1])
            xcontrasttable.append(rxbar[1])
            xmwtable.append(mwxna)
            if mwxchem > 0:
                for j in range(2, len(mwchemarray) + 2):
                    xsldtable.append(rx[j])
                    xcontrasttable.append(rxbar[j])
                    xmwtable.append(mwchemarray[j - 2] / 1.0e3)
        elif mwxchem > 0:
            for j in range(2, len(mwchemarray) + 2):
                xsldtable.append(rx[j])
                xcontrasttable.append(rxbar[j])
                xmwtable.append(mwchemarray[j - 2] / 1.0e3)


#        print xmwtable
#        print xcontrasttable
#        print xsldtable

        xmwvalstring = '\t'.join(['%8.3f' % (x) for x in xmwtable])
        xmwcompstring = '%8.3f' % (mwxcomp)
        xizerostring = '%7.3f' % (izerox)
        avars.izfile.write('#' + xmwvalstring + '\t' +
                           xmwcompstring + '\t' + xizerostring + '\n')
        xsldvalstring = '\t'.join(['%7.3f' % (x) for x in xsldtable])
        xsldcompstring = '%7.3f' % (sldx)
        xrsstring = '%7.3f' % (rsx)
        avars.sldfile.write('#' + xsldvalstring + '\t' +
                            xsldcompstring + '\t' + xrsstring + '\n')
        xcontrastvalstring = '\t'.join(['%7.3f' % (x) for x in xcontrasttable])
        xcontrastcompstring = '%7.3f' % (contrx)
        avars.cotfile.write('#' + xcontrastvalstring +
                            '\t' + xcontrastcompstring + '\n')


# calculate neutron properties of the system:  number of components for
# each molecule type depends on proteinfdval, dnafdval and rnafdval.
# Combine DNA and RNA results for same fract deut.  (It is unlikely that
# there will be both DNA and RNA in a complex, but both will be picked up
# if they exist.)

        avars.izfile.write('\n#NEUTRONS:\n')
        avars.sldfile.write('\n#NEUTRON SLDs:\n')
        avars.cotfile.write('\n#NEUTRON Contrast:\n')

# determine the number of fd2o values

#        print 'NEUTRON'

        fd2o = []
        rs = []
        d2ovals = (100 / float(mvars.d2ostep)) + 1
#        print d2ovals

        for i in range(int(d2ovals)):
            newfd2o = (i * (float(mvars.d2ostep))) / 100
            fd2o.append(newfd2o)

        fd2otitle = 'frac D2O'
#        print fd2otitle
#        print fd2o
#        print len(fd2o)

# determine the solvent from the initial values calculated for the solvent
# + any non-water components

        for i in range(len(fd2o)):
            newrs = (rsh2o + (rsd2o - rsh2o) * fd2o[i]) * 1.0e10
            rs.append(newrs)
#        print rs


# calculate chemical formula parameters
# column headings for SLD and contrast will depend on the number of
# chemical formulas

        chemsld = [[[] for x in range(len(mwchemarray))]
                   for y in range(len(fd2o))]
        chemcontrast = [[[] for x in range(len(mwchemarray))]
                        for y in range(len(fd2o))]
        chemmw = [[[] for x in range(len(mwchemarray))]
                  for y in range(len(fd2o))]
        chemsldtitle = []
        chemcontrasttitle = []
        chemmwtitle = []
        chemmatchtitle = []
        chemmatchpoint = []

        for i in range(len(mwchemarray)):

            #            print 'CHEMICAL FORMULA'

            working_bchem = 0.0
            working_mwchem = 0.0
            working_vtotchem = 0.0
            working_rchem = 0.0
            working_rchembar = 0.0
            temparray = []

            chemsldtitle.append('Molecule ' + str(i + 1))
            chemcontrasttitle.append('Molecule ' + str(i + 1))
            chemmwtitle.append('Molecule ' + str(i + 1) + ' Mw')
            chemmatchtitle.append('Molecule ' + str(i + 1) + ' Match Point')

            for j in range(len(fd2o)):

                working_bchem = (btotchemarray[
                    i] + hexchemarray[i] * fexchchemarray[i] * fd2o[j] * (bd - bh)) * 1.0e-12
                working_mwchem = mwchemarray[
                    i] + fd2o[j] * hexchemarray[i] * fexchchemarray[i]
                working_vtotchem = (working_mwchem) / (denschemarray[i] * na)
                working_rchem = working_bchem / working_vtotchem
                working_rchembar = working_rchem - rs[j]
                chemmw[j][i] = working_mwchem / 1.0e3
                chemsld[j][i] = working_rchem / 1.0e10
                chemcontrast[j][i] = working_rchembar / 1.0e10
                temparray.append(working_rchembar)

# calculate protein match point and print to screen and files

            x = np.array(fd2o)
            y = np.array(temparray)
            slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
                x, y)
            x_intercept = -intercept / slope
#            print slope, intercept, x_intercept
            matchpoint = x_intercept * 100
            chemmatchpoint.append(matchpoint)
            matchstring = '%.2f' % (matchpoint)
            gui_string = chemmatchtitle[i] + ': ' + matchstring + " %D2O\n\n"
            pgui(gui_string)
            avars.izfile.write(
                '#' + chemmatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.sldfile.write(
                '#' + chemmatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.cotfile.write(
                '#' + chemmatchtitle[i] + ': ' + matchstring + " %D2O\n")

#        print chemmwtitle
#        print chemmw
#        print chemsldtitle
#        print chemsld
#        print chemcontrasttitle
#        print chemcontrast
#        print chemmatchtitle
#        print chemmatchpoint


# column headings for SLD and contrast will depend on fract deut
# calculate protein parameters

        protsld = [[[] for x in range(len(mwprotarray))]
                   for y in range(len(fd2o))]
        protcontrast = [[[] for x in range(len(mwprotarray))]
                        for y in range(len(fd2o))]
        protmw = [[[] for x in range(len(mwprotarray))]
                  for y in range(len(fd2o))]
        protsldtitle = []
        protcontrasttitle = []
        protmwtitle = []
        protmatchtitle = []
        protmatchpoint = []

        for i in range(len(mwprotarray)):

            #            print 'PROTEIN'

            working_bprot = 0.0
            working_mwprot = 0.0
            working_vtotprot = 0.0
            working_rprot = 0.0
            working_rprotbar = 0.0
            temparray = []

            percentdeut = 100 * proteinfdval[i]
            if percentdeut == 0:
                protsldtitle.append('Protein')
                protcontrasttitle.append('Protein')
                protmwtitle.append('Protein Mw')
                protmatchtitle.append('Protein Match Point')
            else:
                protsldtitle.append(str(percentdeut) + ' %D Protein')
                protcontrasttitle.append(str(percentdeut) + ' %D Protein')
                protmwtitle.append(str(percentdeut) + ' %D Protein Mw')
                protmatchtitle.append(
                    str(percentdeut) + ' %D Protein Match Point')

            for j in range(len(fd2o)):

                working_bprot = (
                    btotprotarray[i] + hexprotarray[i] * fexchprot * fd2o[j] * (bd - bh)) * 1.0e-12
                working_mwprot = mwprotarray[i] + \
                    fd2o[j] * hexprotarray[i] * fexchprot
                working_vtotprot = (vprot * working_mwprot) / na
                working_rprot = working_bprot / working_vtotprot
                working_rprotbar = working_rprot - rs[j]
                protmw[j][i] = working_mwprot / 1.0e3
                protsld[j][i] = working_rprot / 1.0e10
                protcontrast[j][i] = working_rprotbar / 1.0e10
                temparray.append(working_rprotbar)

# calculate protein match point and print to screen and files

            x = np.array(fd2o)
            y = np.array(temparray)
            slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
                x, y)
            x_intercept = -intercept / slope
#            print slope, intercept, x_intercept
            matchpoint = x_intercept * 100
            protmatchpoint.append(matchpoint)
            matchstring = '%.2f' % (matchpoint)
            gui_string = protmatchtitle[i] + ': ' + matchstring + " %D2O\n\n"
            pgui(gui_string)
            avars.izfile.write(
                '#' + protmatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.sldfile.write(
                '#' + protmatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.cotfile.write(
                '#' + protmatchtitle[i] + ': ' + matchstring + " %D2O\n")

#        print protmwtitle
#        print protmw
#        print protsldtitle
#        print protsld
#        print protcontrasttitle
#        print protcontrast
#        print protmatchtitle
#        print protmatchpoint


# calculate dna params

        dnasld = [[[] for x in range(len(mwdnaarray))]
                  for y in range(len(fd2o))]
        dnacontrast = [[[] for x in range(len(mwdnaarray))]
                       for y in range(len(fd2o))]
        dnamw = [[[] for x in range(len(mwdnaarray))]
                 for y in range(len(fd2o))]
        dnasldtitle = []
        dnacontrasttitle = []
        dnamwtitle = []
        dnamatchtitle = []
        dnamatchpoint = []

        for i in range(len(mwdnaarray)):

            #            print 'DNA'

            working_bdna = 0.0
            working_mwdna = 0.0
            working_vtotdna = 0.0
            working_rdna = 0.0
            working_rdnabar = 0.0
            tempdnaarray = []

            dnapercentdeut = 100 * dnafdval[i]
            if dnapercentdeut == 0:
                dnasldtitle.append('DNA')
                dnacontrasttitle.append('DNA')
                dnamwtitle.append('DNA Mw')
                dnamatchtitle.append('DNA Match Point')
            else:
                dnasldtitle.append(str(dnapercentdeut) + ' %D DNA')
                dnacontrasttitle.append(str(dnapercentdeut) + ' %D DNA')
                dnamwtitle.append(str(dnapercentdeut) + ' %D DNA Mw')
                dnamatchtitle.append(
                    str(dnapercentdeut) + ' %D DNA Match Point')

            for j in range(len(fd2o)):

                working_bdna = (
                    btotdnaarray[i] + hexdnaarray[i] * fexchdna * fd2o[j] * (bd - bh)) * 1.0e-12
                working_mwdna = mwdnaarray[i] + \
                    fd2o[j] * hexdnaarray[i] * fexchdna
                working_vtotdna = (vdna * working_mwdna) / na
                working_rdna = working_bdna / working_vtotdna
                working_rdnabar = working_rdna - rs[j]
                dnamw[j][i] = working_mwdna / 1.0e3
                dnasld[j][i] = working_rdna / 1.0e10
                dnacontrast[j][i] = working_rdnabar / 1.0e10
                tempdnaarray.append(working_rdnabar)

# calculate DNA match point

            x = np.array(fd2o)
            y = np.array(tempdnaarray)
            slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
                x, y)
            x_intercept = -intercept / slope
#        print slope, intercept, x_intercept
            matchpoint = x_intercept * 100
            dnamatchpoint.append(matchpoint)
            matchstring = '%.2f' % (matchpoint)
            gui_string = dnamatchtitle[i] + ': ' + matchstring + " %D2O\n\n"
            pgui(gui_string)
            avars.izfile.write(
                '#' + dnamatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.sldfile.write(
                '#' + dnamatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.cotfile.write(
                '#' + dnamatchtitle[i] + ': ' + matchstring + " %D2O\n")

#        print dnamwtitle
#        print dnamw
#        print dnasldtitle
#        print dnasld
#        print dnacontrasttitle
#        print dnacontrast
#        print dnamatchtitle
#        print dnamatchpoint

# calculate rna params

        rnasld = [[[] for x in range(len(mwrnaarray))]
                  for y in range(len(fd2o))]
        rnacontrast = [[[] for x in range(len(mwrnaarray))]
                       for y in range(len(fd2o))]
        rnamw = [[[] for x in range(len(mwrnaarray))]
                 for y in range(len(fd2o))]
        rnasldtitle = []
        rnacontrasttitle = []
        rnamwtitle = []
        rnamatchtitle = []
        rnamatchpoint = []

        for i in range(len(mwrnaarray)):

            #            print 'RNA'

            working_brna = 0.0
            working_mwrna = 0.0
            working_vtotrna = 0.0
            working_rrna = 0.0
            working_rrnabar = 0.0
            temprnaarray = []

            rnapercentdeut = 100 * rnafdval[i]
            if rnapercentdeut == 0:
                rnasldtitle.append('RNA')
                rnacontrasttitle.append('RNA')
                rnamwtitle.append('RNA Mw')
                rnamatchtitle.append('RNA Match Point')
            else:
                rnasldtitle.append(str(rnapercentdeut) + ' %D RNA')
                rnacontrasttitle.append(str(rnapercentdeut) + ' %D RNA')
                rnamwtitle.append(str(rnapercentdeut) + ' %D RNA Mw')
                rnamatchtitle.append(
                    str(rnapercentdeut) + ' %D RNA Match Point')

            for j in range(len(fd2o)):

                working_brna = (
                    btotrnaarray[i] + hexrnaarray[i] * fexchdna * fd2o[j] * (bd - bh)) * 1.0e-12
                working_mwrna = mwrnaarray[i] + \
                    fd2o[j] * hexrnaarray[i] * fexchdna
                working_vtotrna = (vdna * working_mwrna) / na
                working_rrna = working_brna / working_vtotrna
                working_rrnabar = working_rrna - rs[j]
                rnamw[j][i] = working_mwrna / 1.0e3
                rnasld[j][i] = working_rrna / 1.0e10
                rnacontrast[j][i] = working_rrnabar / 1.0e10
                temprnaarray.append(working_rrnabar)

# calculate RNA match point

            x = np.array(fd2o)
            y = np.array(temprnaarray)
            slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
                x, y)
            x_intercept = -intercept / slope
#            print slope, intercept, x_intercept
            matchpoint = x_intercept * 100
            rnamatchpoint.append(matchpoint)
            matchstring = '%.2f' % (matchpoint)
            gui_string = rnamatchtitle[i] + ': ' + matchstring + " %D2O\n\n"
            pgui(gui_string)
            avars.izfile.write(
                '#' + rnamatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.sldfile.write(
                '#' + rnamatchtitle[i] + ': ' + matchstring + " %D2O\n")
            avars.cotfile.write(
                '#' + rnamatchtitle[i] + ': ' + matchstring + " %D2O\n")

#        print rnamwtitle
#        print rnamw
#        print rnasldtitle
#        print rnasld
#        print rnacontrasttitle
#        print rnacontrast
#        print rnamatchtitle
#        print rnamatchpoint


# calculate complex parameters
# first calculate mwcomp and then use the result to calculate fdna, fprot,
# fchem, sldcomp, contrastcomp, izero as func of fd2o

        mwtotprot = []
        mwtotdna = []
        mwtotrna = []
        mwtotchem = []
        mwcomp = []
        sldcomp = []
        contrastcomp = []

#    calculate total Mw of the complex as a func of fd2o

#        print 'COMPLEX'

        for j in range(len(fd2o)):

            working_mwtotprot = 0.0
            working_mwtotdna = 0.0
            working_mwtotrna = 0.0
            working_mwtotchem = 0.0
            working_mwcomp = 0.0

            for i in range(len(mwprotarray)):
                working_mwtotprot = working_mwtotprot + protmw[j][i]
            mwtotprot.append(working_mwtotprot)

            for i in range(len(mwdnaarray)):
                working_mwtotdna = working_mwtotdna + dnamw[j][i]
            mwtotdna.append(working_mwtotdna)

            for i in range(len(mwrnaarray)):
                working_mwtotrna = working_mwtotrna + rnamw[j][i]
            mwtotrna.append(working_mwtotrna)

            for i in range(len(mwchemarray)):
                working_mwtotchem = working_mwtotchem + chemmw[j][i]
            mwtotchem.append(working_mwtotchem)

            working_mwcomp = mwtotprot[j] + \
                mwtotdna[j] + mwtotrna[j] + mwtotchem[j]
            mwcomp.append(working_mwcomp)

#        print mwtotprot
#        print mwtotdna
#        print mwtotrna
#        print mwtotchem
#        print mwcomp

# calculate Mw fraction vs Mw total for each protein, DNA, RNA and
# molecule component

        fprot = [[[] for x in range(len(mwprotarray))]
                 for y in range(len(fd2o))]
        fdna = [[[] for x in range(len(mwdnaarray))] for y in range(len(fd2o))]
        frna = [[[] for x in range(len(mwrnaarray))] for y in range(len(fd2o))]
        fchem = [[[] for x in range(len(mwchemarray))]
                 for y in range(len(fd2o))]

        for j in range(len(fd2o)):

            for i in range(len(mwprotarray)):
                fprot[j][i] = protmw[j][i] / mwcomp[j]

            for i in range(len(mwdnaarray)):
                fdna[j][i] = dnamw[j][i] / mwcomp[j]

            for i in range(len(mwrnaarray)):
                frna[j][i] = rnamw[j][i] / mwcomp[j]

            for i in range(len(mwchemarray)):
                fchem[j][i] = chemmw[j][i] / mwcomp[j]

#        print fprot
#        print fdna
#        print frna
#        print fchem


# first calculate parameters for protein, DNA, RNA and molecule component
# and then add them up for each fd2o

        sldtotprot = []
        sldtotdna = []
        sldtotrna = []
        sldtotchem = []
        sldcomp = []
        izerocomp = []
        sqrtizerocomp = []
        contrastcomp = []
        contrasttotprot = []
        contrasttotdna = []
        contrasttotrna = []
        contrasttotchem = []

        for j in range(len(fd2o)):

            working_rtotprot = 0.0
            working_rbartotprot = 0.0
            working_rtotdna = 0.0
            working_rbartotdna = 0.0
            working_rtotrna = 0.0
            working_rbartotrna = 0.0
            working_rtotchem = 0.0
            working_rbartotchem = 0.0
            working_vchem = 0.0
            working_rcomp = 0.0
            working_rbarcomp = 0.0
            working_izero = 0.0
            working_sqrti = 0.0

            rs[j] = rs[j] / 1.0e10                          # for printing out

            for i in range(len(mwprotarray)):
                working_rtotprot = working_rtotprot + \
                    fprot[j][i] * protsld[j][i]
                working_rbartotprot = working_rbartotprot + \
                    fprot[j][i] * protcontrast[j][i]
            sldtotprot.append(working_rtotprot)
            contrasttotprot.append(working_rbartotprot)

            for i in range(len(mwdnaarray)):
                working_rtotdna = working_rtotdna + fdna[j][i] * dnasld[j][i]
                working_rbartotdna = working_rbartotdna + \
                    fdna[j][i] * dnacontrast[j][i]
            sldtotdna.append(working_rtotdna)
            contrasttotdna.append(working_rbartotdna)

            for i in range(len(mwrnaarray)):
                working_rtotrna = working_rtotrna + frna[j][i] * rnasld[j][i]
                working_rbartotrna = working_rbartotrna + \
                    frna[j][i] * rnacontrast[j][i]
            sldtotrna.append(working_rtotrna)
            contrasttotrna.append(working_rbartotrna)

            for i in range(len(mwchemarray)):
                working_rtotchem = working_rtotchem + \
                    fchem[j][i] * chemsld[j][i]
                working_rbartotchem = working_rbartotchem + \
                    fchem[j][i] * chemcontrast[j][i]
                working_vchem = 1.0 / denschemarray[i]
            sldtotchem.append(working_rtotchem)
            contrasttotchem.append(working_rbartotchem)

            working_rcomp = sldtotprot[j] + \
                sldtotdna[j] + sldtotrna[j] + sldtotchem[j]
            working_rbarcomp = contrasttotprot[
                j] + contrasttotdna[j] + contrasttotrna[j] + contrasttotchem[j]
            sldcomp.append(working_rcomp)
            contrastcomp.append(working_rbarcomp)

            working_izero = working_conc * mwcomp[j] * 1.0e3 / na * (vprot * contrasttotprot[j] * 1.0e10 + vdna * (
                contrasttotdna[j] + contrasttotrna[j]) * 1.0e10 + working_vchem * contrasttotchem[j] * 1.0e10)**2

            working_sqrti = math.sqrt(working_izero)

            if contrastcomp[j] < 0.0:
                working_sqrti = -working_sqrti

            izerocomp.append(working_izero)
            sqrtizerocomp.append(working_sqrti)


#        print 'SLD'
#        print sldtotprot
#        print sldtotdna
#        print sldtotrna
#        print sldtotchem
#        print sldcomp
#        print 'Contrast'
#        print contrasttotprot
#        print contrasttotdna
#        print contrasttotrna
#        print contrasttotchem
#        print contrastcomp
#        print 'izero'
#        print izerocomp
#        print sqrtizerocomp

# calculate Complex match point

        x = np.array(fd2o)
        y = np.array(contrastcomp)
        slope, intercept, r_value, p_value, slope_std_error = stats.linregress(
            x, y)
        x_intercept = -intercept / slope
#        print slope, intercept, x_intercept
        matchpoint = x_intercept * 100
        matchstring = '%.2f' % (matchpoint)
        gui_string = 'Complex Match Point: ' + matchstring + " %D2O\n\n"
        pgui(gui_string)
        avars.izfile.write('#Complex Match Point: ' + matchstring + " %D2O\n")
        avars.sldfile.write('#Complex Match Point: ' + matchstring + " %D2O\n")
        avars.cotfile.write('#Complex Match Point: ' + matchstring + " %D2O\n")


# print out fraction of exchangeable protein, nucleic acid and molecule
# hydrogens that actually do exchange to output files

        fexchp_string = '%.2f' % (mvars.fexchp)
        fexchn_string = '%.2f' % (mvars.fexchn)
        avars.izfile.write('\n#Fraction of exchanged protein hydrogens: ' +
                           fexchp_string + '\n')
        avars.izfile.write(
            '#Fraction of exchanged nucleic acid hydrogens: ' + fexchn_string + '\n')
        avars.sldfile.write(
            '\n#Fraction of exchanged protein hydrogens: ' + fexchp_string + '\n')
        avars.sldfile.write(
            '#Fraction of exchanged nucleic acid hydrogens: ' + fexchn_string + '\n')
        avars.cotfile.write(
            '\n#Fraction of exchanged protein hydrogens: ' + fexchp_string + '\n')
        avars.cotfile.write(
            '#Fraction of exchanged nucleic acid hydrogens: ' + fexchn_string + '\n')

        for i in range(len(mwchemarray)):
            fexchc = fexchchemarray[i]
            fexchc_string = '%.2f' % (fexchc)
            avars.izfile.write('#Fraction of exchanged molecule ' +
                               str(i + 1) + ' hydrogens: ' + fexchc_string + '\n')
            avars.sldfile.write('#Fraction of exchanged molecule ' +
                                str(i + 1) + ' hydrogens: ' + fexchc_string + '\n')
            avars.cotfile.write('#Fraction of exchanged molecule ' +
                                str(i + 1) + ' hydrogens: ' + fexchc_string + '\n')


# combine protein, dna, rna and molecule matrices to create contrast, sld
# and izero tables to print to files

        if len(mwprotarray) > 0:
            sldtitle = protsldtitle[0]
            contrasttitle = protcontrasttitle[0]
            mwtitle = protmwtitle[0]
            for i in range(1, len(mwprotarray)):
                sldtitle = sldtitle + ', ' + protsldtitle[i]
                contrasttitle = contrasttitle + ', ' + protcontrasttitle[i]
                mwtitle = mwtitle + ', ' + protmwtitle[i]
            if len(mwdnaarray) > 0:
                for i in range(len(mwdnaarray)):
                    sldtitle = sldtitle + ', ' + dnasldtitle[i]
                    contrasttitle = contrasttitle + ', ' + dnacontrasttitle[i]
                    mwtitle = mwtitle + ', ' + dnamwtitle[i]
            if len(mwrnaarray) > 0:
                for i in range(len(mwrnaarray)):
                    sldtitle = sldtitle + ', ' + rnasldtitle[i]
                    contrasttitle = contrasttitle + ', ' + rnacontrasttitle[i]
                    mwtitle = mwtitle + ', ' + rnamwtitle[i]
            if len(mwchemarray) > 0:
                for i in range(len(mwchemarray)):
                    sldtitle = sldtitle + ', ' + chemsldtitle[i]
                    contrasttitle = contrasttitle + ', ' + chemcontrasttitle[i]
                    mwtitle = mwtitle + ', ' + chemmwtitle[i]

        elif len(mwdnaarray) > 0:
            sldtitle = dnasldtitle[0]
            contrasttitle = dnacontrasttitle[0]
            mwtitle = dnamwtitle[0]
            for i in range(1, len(mwdnaarray)):
                sldtitle = sldtitle + ', ' + dnasldtitle[i]
                contrasttitle = contrasttitle + ', ' + dnacontrasttitle[i]
                mwtitle = mwtitle + ', ' + dnamwtitle[i]
            if len(mwrnaarray) > 0:
                for i in range(len(mwrnaarray)):
                    sldtitle = sldtitle + ', ' + rnasldtitle[i]
                    contrasttitle = contrasttitle + ', ' + rnacontrasttitle[i]
                    mwtitle = mwtitle + ', ' + rnamwtitle[i]
            if len(mwchemarray) > 0:
                for i in range(len(mwchemarray)):
                    sldtitle = sldtitle + ', ' + chemsldtitle[i]
                    contrasttitle = contrasttitle + ', ' + chemcontrasttitle[i]
                    mwtitle = mwtitle + ', ' + chemmwtitle[i]

        elif len(mwrnaarray) > 0:
            sldtitle = rnasldtitle[0]
            contrasttitle = rnacontrasttitle[0]
            mwtitle = rnamwtitle[0]
            for i in range(1, len(mwrnaarray)):
                sldtitle = sldtitle + ', ' + rnasldtitle[i]
                contrasttitle = contrasttitle + ', ' + rnacontrasttitle[i]
                mwtitle = mwtitle + ', ' + rnamwtitle[i]
            if len(mwchemarray) > 0:
                for i in range(len(mwchemarray)):
                    sldtitle = sldtitle + ', ' + chemsldtitle[i]
                    contrasttitle = contrasttitle + ', ' + chemcontrasttitle[i]
                    mwtitle = mwtitle + ', ' + chemmwtitle[i]

        elif len(mwchemarray) > 0:
            sldtitle = chemsldtitle[0]
            contrasttitle = chemcontrasttitle[0]
            mwtitle = chemmwtitle[0]
            for i in range(1, len(mwchemarray)):
                sldtitle = sldtitle + ', ' + chemsldtitle[i]
                contrasttitle = contrasttitle + ', ' + chemcontrasttitle[i]
                mwtitle = mwtitle + ', ' + chemmwtitle[i]

#        print sldtitle
#        print contrasttitle
#        print mwtitle

        avars.sldfile.write('\n#NEUTRON SLDs (10^10 cm^-2):\n')
        avars.cotfile.write('\n#NEUTRON Contrast (10^10 cm^-2):\n')
        avars.izfile.write('\n# frac D2O, ' + mwtitle +
                           ', Complex Mw (kDa), I(0) (cm^-1), sqrtI(0)\n')
        avars.sldfile.write('# frac D2O, ' + sldtitle + ', Complex, Solvent\n')
        avars.cotfile.write('# frac D2O, ' + contrasttitle + ', Complex\n')

        if len(mwprotarray) > 0:
            sldtable = protsld
            contrasttable = protcontrast
            mwtable = protmw
            for j in range(len(fd2o)):
                if len(mwdnaarray) > 0:
                    sldtable[j].extend(dnasld[j])
                    contrasttable[j].extend(dnacontrast[j])
                    mwtable[j].extend(dnamw[j])
                if len(mwrnaarray) > 0:
                    sldtable[j].extend(rnasld[j])
                    contrasttable[j].extend(rnacontrast[j])
                    mwtable[j].extend(rnamw[j])
                if len(mwchemarray) > 0:
                    sldtable[j].extend(chemsld[j])
                    contrasttable[j].extend(chemcontrast[j])
                    mwtable[j].extend(chemmw[j])
        elif len(mwdnaarray) > 0:
            sldtable = dnasld
            contrasttable = dnacontrast
            mwtable = dnamw
            for j in range(len(fd2o)):
                if len(mwrnaarray) > 0:
                    sldtable[j].extend(rnasld[j])
                    contrasttable[j].extend(rnacontrast[j])
                    mwtable[j].extend(rnamw[j])
                if len(mwchemarray) > 0:
                    sldtable[j].extend(chemsld[j])
                    contrasttable[j].extend(chemcontrast[j])
                    mwtable[j].extend(chemmw[j])
        elif len(mwrnaarray) > 0:
            sldtable = rnasld
            contrasttable = rnacontrast
            mwtable = rnamw
            for j in range(len(fd2o)):
                if len(mwchemarray) > 0:
                    sldtable[j].extend(chemsld[j])
                    contrasttable[j].extend(chemcontrast[j])
                    mwtable[j].extend(chemmw[j])
        elif len(mwchemarray) > 0:
            sldtable = chemsld
            contrasttable = chemcontrast
            mwtable = chemmw

        for j in range(len(fd2o)):
            fd2ostring = '%5.2f' % (fd2o[j])
            mwvalstring = '\t'.join(['%8.3f' % (x) for x in mwtable[j]])
            mwcompstring = '%8.3f' % (mwcomp[j])
            izerostring = '%7.3f' % (izerocomp[j])
            sqrtizerostring = '%7.3f' % (sqrtizerocomp[j])
            avars.izfile.write(fd2ostring + '\t' + mwvalstring + '\t' + mwcompstring +
                               '\t' + izerostring + '\t' + sqrtizerostring + '\n')
            sldvalstring = '\t'.join(['%7.3f' % (x) for x in sldtable[j]])
            sldcompstring = '%7.3f' % (sldcomp[j])
            rsstring = '%7.3f' % (rs[j])
            avars.sldfile.write(fd2ostring + '\t' + sldvalstring + '\t' +
                                sldcompstring + '\t' + rsstring + '\n')
            contrastvalstring = '\t'.join(
                ['%7.3f' % (x) for x in contrasttable[j]])
            contrastcompstring = '%7.3f' % (contrastcomp[j])
            avars.cotfile.write(fd2ostring + '\t' + contrastvalstring +
                                '\t' + contrastcompstring + '\n')

        gui_string = "\nFiles " + avars.izerofile + ", " + avars.scatlendenfile + \
            " and " + avars.contrastfile + " written to ./" + avars.contpath + "."
        pgui(gui_string)

# plot only if there are two components or less

        if (len(mwprotarray) + len(mwdnaarray) + len(mwrnaarray) + len(mwchemarray)) > 2:
            mvars.plotflag = 0

        if (mvars.plotflag == 1):
            graph = Gnuplot.Gnuplot(debug=1)
            graph.clear()
            graph('set title "SqrtI(0) vs D2O Fraction"')
            graph('set zeroaxis')
            graph('set xtics 0,.05,1')
            graph.xlabel('D2O Fraction')
            graph.ylabel('sqrtI(0)')
            graph2 = Gnuplot.Gnuplot(debug=1)
            graph2.clear()
            graph2('set title "Scattering Length Density vs D2O Fraction"')
            graph2('set zeroaxis')
            graph2('set xtics 0,.05,1')
            graph2.xlabel('D2O Fraction')
            graph2.ylabel('SLD x 10^10 cm^-2')
            graph3 = Gnuplot.Gnuplot(debug=1)
            graph3.clear()
            graph3('set title "Contrast vs D2O Fraction"')
            graph3('set zeroaxis')
            graph3('set xtics 0,.05,1')
            graph3.xlabel('D2O Fraction')
            graph3.ylabel('Contrast x 10^10 cm^-2')
            graph4 = Gnuplot.Gnuplot(debug=1)
            graph4.clear()
            graph4('set title "I(0) vs D2O Fraction"')
            graph4('set xtics 0,.05,1')
            graph4.xlabel('D2O Fraction')
            graph4.ylabel('I(0) cm^-1')

# arrays for plots
            izerocoords = []
            sizerocoords = []
            rscoords = []
            rcomp1coords = []
            rcomp2coords = []
            ccomp1coords = []
            ccomp2coords = []
            rcomplexcoords = []
            ccomplexcoords = []

            for j in range(len(fd2o)):
                izerocoords.append([fd2o[j], izerocomp[j]])
                sizerocoords.append([fd2o[j], sqrtizerocomp[j]])
                rscoords.append([fd2o[j], rs[j]])
                rcomplexcoords.append([fd2o[j], sldcomp[j]])
                ccomplexcoords.append([fd2o[j], contrastcomp[j]])

                if len(mwprotarray) == 2:
                    rcomp1coords.append([fd2o[j], protsld[j][0]])
                    ccomp1coords.append([fd2o[j], protcontrast[j][0]])
                    comp1contrasttitle = protcontrasttitle[0]
                    comp1sldtitle = protsldtitle[0]
                    rcomp2coords.append([fd2o[j], protsld[j][1]])
                    ccomp2coords.append([fd2o[j], protcontrast[j][1]])
                    comp2contrasttitle = protcontrasttitle[1]
                    comp2sldtitle = protsldtitle[1]
                    continue
                elif len(mwprotarray) == 1:
                    rcomp1coords.append([fd2o[j], protsld[j][0]])
                    ccomp1coords.append([fd2o[j], protcontrast[j][0]])
                    comp1contrasttitle = protcontrasttitle[0]
                    comp1sldtitle = protsldtitle[0]
                    if len(mwdnaarray) == 1:
                        rcomp2coords.append([fd2o[j], dnasld[j][0]])
                        ccomp2coords.append([fd2o[j], dnacontrast[j][0]])
                        comp2contrasttitle = dnacontrasttitle[0]
                        comp2sldtitle = dnasldtitle[0]
                    elif len(mwrnaarray) == 1:
                        rcomp2coords.append([fd2o[j], rnasld[j][0]])
                        ccomp2coords.append([fd2o[j], rnacontrast[j][0]])
                        comp2contrasttitle = rnacontrasttitle[0]
                        comp2sldtitle = rnasldtitle[0]
                    elif len(mwchemarray) == 1:
                        rcomp2coords.append([fd2o[j], chemsld[j][0]])
                        ccomp2coords.append([fd2o[j], chemcontrast[j][0]])
                        comp2contrasttitle = chemcontrasttitle[0]
                        comp2sldtitle = chemsldtitle[0]
                    else:
                        rcomp2coords.append([fd2o[j], 0.0])
                        ccomp2coords.append([fd2o[j], 0.0])
                        comp2contrasttitle = 'None'
                        comp2sldtitle = 'None'
                    continue
                elif len(mwdnaarray) == 2:
                    rcomp1coords.append([fd2o[j], dnasld[j][0]])
                    ccomp1coords.append([fd2o[j], dnacontrast[j][0]])
                    comp1contrasttitle = dnacontrasttitle[0]
                    comp1sldtitle = dnasldtitle[0]
                    rcomp2coords.append([fd2o[j], dnasld[j][1]])
                    ccomp2coords.append([fd2o[j], dnacontrast[j][1]])
                    comp2contrasttitle = dnacontrasttitle[1]
                    comp2sldtitle = dnasldtitle[1]
                    continue
                elif len(mwdnaarray) == 1:
                    rcomp1coords.append([fd2o[j], dnasld[j][0]])
                    ccomp1coords.append([fd2o[j], dnacontrast[j][0]])
                    comp1contrasttitle = dnacontrasttitle[0]
                    comp1sldtitle = dnasldtitle[0]
                    if len(mwrnaarray) == 1:
                        rcomp2coords.append([fd2o[j], rnasld[j][0]])
                        ccomp2coords.append([fd2o[j], rnacontrast[j][0]])
                        comp2contrasttitle = rnacontrasttitle[0]
                        comp2sldtitle = rnasldtitle[0]
                    elif len(mwchemarray) == 1:
                        rcomp2coords.append([fd2o[j], chemsld[j][0]])
                        ccomp2coords.append([fd2o[j], chemcontrast[j][0]])
                        comp2contrasttitle = chemcontrasttitle[0]
                        comp2sldtitle = chemsldtitle[0]
                    else:
                        rcomp2coords.append([fd2o[j], 0.0])
                        ccomp2coords.append([fd2o[j], 0.0])
                        comp2contrasttitle = 'None'
                        comp2sldtitle = 'None'
                    continue
                elif len(mwrnaarray) == 2:
                    rcomp1coords.append([fd2o[j], rnasld[j][0]])
                    ccomp1coords.append([fd2o[j], rnacontrast[j][0]])
                    comp1contrasttitle = rnacontrasttitle[0]
                    comp1sldtitle = rnasldtitle[0]
                    rcomp2coords.append([fd2o[j], rnasld[j][1]])
                    ccomp2coords.append([fd2o[j], rnacontrast[j][1]])
                    comp2contrasttitle = rnacontrasttitle[1]
                    comp2sldtitle = rnasldtitle[1]
                elif len(mwrnaarray) == 1:
                    rcomp1coords.append([fd2o[j], rnasld[j][0]])
                    ccomp1coords.append([fd2o[j], rnacontrast[j][0]])
                    comp1contrasttitle = rnacontrasttitle[0]
                    comp1sldtitle = rnasldtitle[0]
                    if len(mwchemarray) == 1:
                        rcomp2coords.append([fd2o[j], chemsld[j][0]])
                        ccomp2coords.append([fd2o[j], chemcontrast[j][0]])
                        comp2contrasttitle = chemcontrasttitle[0]
                        comp2sldtitle = chemsldtitle[0]
                    else:
                        rcomp2coords.append([fd2o[j], 0.0])
                        ccomp2coords.append([fd2o[j], 0.0])
                        comp2contrasttitle = 'None'
                        comp2sldtitle = 'None'
                    continue
                elif len(mwchemarray) == 2:
                    rcomp1coords.append([fd2o[j], chemsld[j][0]])
                    ccomp1coords.append([fd2o[j], chemcontrast[j][0]])
                    comp1contrasttitle = chemcontrasttitle[0]
                    comp1sldtitle = chemsldtitle[0]
                    rcomp2coords.append([fd2o[j], chemsld[j][1]])
                    ccomp2coords.append([fd2o[j], chemcontrast[j][1]])
                    comp2contrasttitle = chemcontrasttitle[1]
                    comp2sldtitle = chemsldtitle[1]
                    continue
                elif len(mwchemarray) == 1:
                    rcomp1coords.append([fd2o[j], chemsld[j][0]])
                    ccomp1coords.append([fd2o[j], chemcontrast[j][0]])
                    comp1contrasttitle = chemcontrasttitle[0]
                    comp1sldtitle = chemsldtitle[0]
                    rcomp2coords.append([fd2o[j], 0.0])
                    ccomp2coords.append([fd2o[j], 0.0])
                    comp2contrasttitle = 'None'
                    comp2sldtitle = 'None'

            try:
                graph.plot(Gnuplot.Data(
                    sizerocoords, using='1:2 w lines lw 2', title=''))
            except:
                message = 'error trying to plot sqrtI(0) vs D2O Fraction data'
                message += ' :  stopping here'
                pgui(message)

            try:
                graph2.plot(Gnuplot.Data(rscoords, using='1:2 w lines lw 2', title='Solvent'),
                            Gnuplot.Data(
                                rcomp1coords, using='1:2 w lines lw 2', title=comp1sldtitle),
                            Gnuplot.Data(
                                rcomp2coords, using='1:2 w lines lw 2', title=comp2sldtitle),
                            Gnuplot.Data(rcomplexcoords, using='1:2 w lines lw 2', title='Complex'))
            except:
                message = 'error trying plot Scattering Length Density vs D2O Fraction data'
                message += ' :  stopping here'
                pgui(message)

            try:
                graph3.plot(Gnuplot.Data(ccomp1coords, using='1:2 w lines lw 2', title=comp1contrasttitle),
                            Gnuplot.Data(
                                ccomp2coords, using='1:2 w lines lw 2', title=comp2contrasttitle),
                            Gnuplot.Data(ccomplexcoords, using='1:2 w lines lw 2', title='Complex'))
            except:
                message = 'error trying plot Contrast vs D2O Fraction data'
                message += ' :  stopping here'
                pgui(message)

            try:
                graph4.plot(Gnuplot.Data(
                    izerocoords, using='1:2 w lp pt 5 lw 2', title=''))
            except:
                message = 'error trying plot I(0) vs D2O Fraction data'
                message += ' :  stopping here'
                pgui(message)

        else:
            pgui(
                "\n\nResults are not plotted for complexes with more than 2 components.")

        return

    def epilogue(self):
        '''
        method to print out results and to move results
        to appropriate places.
        '''

        log = self.log
        pgui = self.run_utils.print_gui

        log.debug('in epilogue')

        self.run_utils.clean_up(log)

        st = ''.join(['=' for x in xrange(60)])
        pgui("\n%s \n\n" % (st))
        fraction_done = 1
        report_string = 'STATUS\t' + str(fraction_done)
        pgui(report_string)

        pgui('\nCONTRAST_CALCULATOR IS DONE')

        time.sleep(1.0)

        return
