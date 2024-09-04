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

#	SASUTIL
#
#	12/10/2009	--	initial coding				:	jc
#	11/24/2011	-- 	moved to seperate file			:	jc	
#	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	SasUtil holds general methods for file naming, os differences,
        chemical formula parsing, etc.
'''
import sasmol.properties as properties
import string
   

NAME, NUM, LPAREN, RPAREN, EOS = list(range(5))
import re
_lexer = re.compile(r"[A-Z][a-z]*|\d+|[()]|<EOS>").match
del re



def get_full_filename(path,filename,**kwargs):
	
	web_flag = False

	if 'web_path' in kwargs:				# dict kwarg['web_path'] = web_path
		web_flag = True
		web_path = kwargs['web_path']
		print(('web_path = ',web_path))		

	return

# symbol, name, atomic number, molecular weight
_data = r"""'Ac', 'Actinium', 89, 227
'Ag', 'Silver', 47, 107.868
'Al', 'Aluminum', 13, 26.98154
'Am', 'Americium', 95, 243
'Ar', 'Argon', 18, 39.948
'As', 'Arsenic', 33, 74.9216
'At', 'Astatine', 85, 210
'Au', 'Gold', 79, 196.9665
'B', 'Boron', 5, 10.81
'Ba', 'Barium', 56, 137.33
'Be', 'Beryllium', 4, 9.01218
'Bi', 'Bismuth', 83, 208.9804
'Bk', 'Berkelium', 97, 247
'Br', 'Bromine', 35, 79.904
'C', 'Carbon', 6, 12.011
'Ca', 'Calcium', 20, 40.08
'Cd', 'Cadmium', 48, 112.41
'Ce', 'Cerium', 58, 140.12
'Cf', 'Californium', 98, 251
'Cl', 'Chlorine', 17, 35.453
'Cm', 'Curium', 96, 247
'Co', 'Cobalt', 27, 58.9332
'Cr', 'Chromium', 24, 51.996
'Cs', 'Cesium', 55, 132.9054
'Cu', 'Copper', 29, 63.546
'Dy', 'Dysprosium', 66, 162.50
'Er', 'Erbium', 68, 167.26
'Es', 'Einsteinium', 99, 254
'Eu', 'Europium', 63, 151.96
'F', 'Fluorine', 9, 18.998403
'Fe', 'Iron', 26, 55.847
'Fm', 'Fermium', 100, 257
'Fr', 'Francium', 87, 223
'Ga', 'Gallium', 31, 69.735
'Gd', 'Gadolinium', 64, 157.25
'Ge', 'Germanium', 32, 72.59
'H', 'Hydrogen', 1, 1.0079
'He', 'Helium', 2, 4.0026
'Hf', 'Hafnium', 72, 178.49
'Hg', 'Mercury', 80, 200.59
'Ho', 'Holmium', 67, 164.9304
'I', 'Iodine', 53, 126.9045
'In', 'Indium', 49, 114.82
'Ir', 'Iridium', 77, 192.22
'K', 'Potassium', 19, 39.0983
'Kr', 'Krypton', 36, 83.80
'La', 'Lanthanum', 57, 138.9055
'Li', 'Lithium', 3, 6.94
'Lr', 'Lawrencium', 103, 260
'Lu', 'Lutetium', 71, 174.96
'Md', 'Mendelevium', 101, 258
'Mg', 'Magnesium', 12, 24.305
'Mn', 'Manganese', 25, 54.9380
'Mo', 'Molybdenum', 42, 95.94
'N', 'Nitrogen', 7, 14.0067
'Na', 'Sodium', 11, 22.98977
'Nb', 'Niobium', 41, 92.9064
'Nd', 'Neodymium', 60, 144.24
'Ne', 'Neon', 10, 20.17
'Ni', 'Nickel', 28, 58.71
'No', 'Nobelium', 102, 259
'Np', 'Neptunium', 93, 237.0482
'O', 'Oxygen', 8, 15.9994
'Os', 'Osmium', 76, 190.2
'P', 'Phosphorous', 15, 30.97376
'Pa', 'Proactinium', 91, 231.0359
'Pb', 'Lead', 82, 207.2
'Pd', 'Palladium', 46, 106.4
'Pm', 'Promethium', 61, 145
'Po', 'Polonium', 84, 209
'Pr', 'Praseodymium', 59, 140.9077
'Pt', 'Platinum', 78, 195.09
'Pu', 'Plutonium', 94, 244
'Ra', 'Radium', 88, 226.0254
'Rb', 'Rubidium', 37, 85.467
'Re', 'Rhenium', 75, 186.207
'Rh', 'Rhodium', 45, 102.9055
'Rn', 'Radon', 86, 222
'Ru', 'Ruthenium', 44, 101.07
'S', 'Sulfur', 16, 32.06
'Sb', 'Antimony', 51, 121.75
'Sc', 'Scandium', 21, 44.9559
'Se', 'Selenium', 34, 78.96
'Si', 'Silicon', 14, 28.0855
'Sm', 'Samarium', 62, 150.4
'Sn', 'Tin', 50, 118.69
'Sr', 'Strontium', 38, 87.62
'Ta', 'Tantalum', 73, 180.947
'Tb', 'Terbium', 65, 158.9254
'Tc', 'Technetium', 43, 98.9062
'Te', 'Tellurium', 52, 127.60
'Th', 'Thorium', 90, 232.0381
'Ti', 'Titanium', 22, 47.90
'Tl', 'Thallium', 81, 204.37
'Tm', 'Thulium', 69, 168.9342
'U', 'Uranium', 92, 238.029
'Unh', 'Unnilhexium', 106, 263
'Unp', 'Unnilpentium', 105, 260
'Unq', 'Unnilquadium', 104, 260
'Uns', 'Unnilseptium', 107, 262
'V', 'Vanadium', 23, 50.9415
'W', 'Tungsten', 74, 183.85
'Xe', 'Xenon', 54, 131.30
'Y', 'Yttrium', 39, 88.9059
'Yb', 'Ytterbium', 70, 173.04
'Zn', 'Zinc', 30, 65.38
'Zr', 'Zirconium', 40, 91.22
'D', 'Deuterium', 1, 2.158"""


class Element:
        def __init__(self, symbol,name,atomicnumber,molweight):
            self.sym = symbol
            self.name = name
            self.ano = atomicnumber
            self.mw = molweight

        def getweight(self):
            return self.mw

        def addsyms(self, weight, result):
            result[self.sym] = result.get(self.sym, 0) + weight

def build_dict(s):
    answer = {}
    for line in s.split("\n"):
        symbol, name, num, weight = eval(line)
        answer[symbol] = Element(symbol, name, num, weight)
    return answer

class ElementSequence:
        def __init__(self, *seq):
            self.seq = list(seq)
            self.count = 1

        def append(self, thing):
            self.seq.append(thing)

        def getweight(self):
            sum = 0.0
            for thing in self.seq:
                sum = sum + thing.getweight()
            return sum * self.count

        def set_count(self, n):
            self.count = n

        def __len__(self):
            return len(self.seq)

        def addsyms(self, weight, result):
            totalweight = weight * self.count
            for thing in self.seq:
                thing.addsyms(totalweight, result)

        def displaysyms(self,sym2elt):
            result = {}
            self.addsyms(1, result)
            items = list(result.items())
            items.sort()
            for sym, count in items:
                print((sym, " : ",count))


class Tokenizer:
        def __init__(self, input):
            self.input = input + "<EOS>"
            self.i = 0

        def gettoken(self):
            global ttype, tvalue
            self.lasti = self.i
            m = _lexer(self.input, self.i)
            if m is None:
                self.error("unexpected character")
            self.i = m.end()
            tvalue = m.group()
            if tvalue == "(":
                ttype = LPAREN
            elif tvalue == ")":
                ttype = RPAREN
            elif tvalue == "<EOS>":
                ttype = EOS
            elif "0" <= tvalue[0] <= "9":
                ttype = NUM
                tvalue = int(tvalue)
            else:
                ttype = NAME

        def error(self, msg):
            emsg = msg + ":\n"
            emsg = emsg + self.input[:-5] + "\n"  # strip <EOS>
            emsg = emsg + " " * self.lasti + "^\n"
            raise ValueError(emsg)

def parse(s,sym2elt):
        global t, ttype, tvalue
        t = Tokenizer(s)
        t.gettoken()
        seq = parse_sequence(sym2elt)
        if ttype != EOS:
            t.error("expected end of input")
        return seq

def parse_sequence(sym2elt):
        global t, ttype, tvalue
        seq = ElementSequence()
        while ttype in (LPAREN, NAME):
        # parenthesized expression or straight name
            if ttype == LPAREN:
                t.gettoken()
                thisguy = parse_sequence(sym2elt)
                if ttype != RPAREN:
                    t.error("expected right paren")
                t.gettoken()
            else:
                assert ttype == NAME
                if tvalue in sym2elt:
                    thisguy = ElementSequence(sym2elt[tvalue])
                else:
                    t.error("'" + tvalue + "' is not an element symbol")
                t.gettoken()
        # followed by optional count
            if ttype == NUM:
                thisguy.set_count(tvalue)
                t.gettoken()
            seq.append(thisguy)
        if len(seq) == 0:
            t.error("empty sequence")
        return seq



def get_chemical_formula(formula_string):

    Atomic = properties.Atomic()
    #standard_atomic_weights = Atomic.amu(keep_lower_case=True) 
    amu = Atomic.amu(keep_lower_case=True) 
    sym2elt = build_dict(_data)

    #sym2elt = amu

    formula_dictionary = {}

    #print amu

    error = []

    try:
        seq = parse(formula_string.strip(" "),sym2elt)
        #seq.displaysyms(sym2elt)
        seq.addsyms(1,formula_dictionary)
        items = list(formula_dictionary.items())
        items.sort()

        #for sym, count in items:
        #    print sym," :: ",count

    except (ValueError, detail):
        print((str(detail)))
        error.append(detail)

    return error,formula_dictionary


if __name__ == "__main__":

    path = './'

    filename = 'run_0'

    get_full_filename(path,filename,web_path="/Users/curtisj")

    '''
    while 1:
        x = raw_input("? ")
        fields = string.split(x)
        if len(fields) != 2:
            print("input must have two fields")
            continue
        action, formula = fields
        ok = 0
        try:
            seq = parse(formula,sym2elt)
            ok = 1
        except ValueError, detail:
            print(str(detail))
        if not ok:
            continue
        if action == "molw":
            print("molecular weight", seq.getweight())
        elif action == "syms":
            seq.displaysyms(sym2elt)
        else:
            print("unknown action:", action)
    '''





