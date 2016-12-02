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
import sys,time,os,string,pwd
from Tkinter import *
from sassie.gui.aprogTable import aprogs
from tkMessageBox import *
import tkFont

def user_parameters():

	desired_modules =	{
		"Center Frames"	:	"YES",
		"Align Frames"	:	"YES",
		"Data Interpolation"	:	"YES",
		"Coordinate Tools"	:	"YES",
		"Merge Utilities"	:	"YES",
		"Monomer Monte Carlo"	:	"YES",
		"Complex Monte Carlo"	:	"YES",
		"Structure Minimization"	:	"YES",
		"Structure Open Minimization"	:	"NO",
		"Two-Body Grid"	:	"YES",
		"Xtal2sas"	:	"YES",
		"Cryson"	:	"YES",
		"Crysol"	:	"YES",
		"Scattering Length Density"	:	"NO",
		"Hydropro"	:	"NO",
		"EM to SANS"	:	"YES",
		"Chi-Square Filter"	:	"YES",
		"Density Plot"	:	"YES",
		"APBS"	:	"NO"
				}	
	
	desired_font_sizes = 	{

			'default_font'	: 'YES',
			'button_font'	: ('Helvetica',8,'bold'),	
			'local_font'	: ('Helvetica',8,'bold'),
			'button_pady'	: 6

				}

	enable_saslock = 'YES'

	return desired_modules,desired_font_sizes,enable_saslock

def write_saslock():

	saslockfile = open('.saslock','w')
	saslockfile.write('SASSIE started at : '+time.ctime()+'\n')
	saslockfile.write('IN DIRECTORY : '+os.getcwd()+'\n')
	saslockfile.write('BY USER : '+pwd.getpwuid(os.getuid()).pw_name+'\n')
	saslockfile.close()

	return

class Prog(Frame):
	def __init__(self,parent=None):
		Frame.__init__(self,parent)
		self.config(bg="white")
		self.master.title("SASSIE")
		self.pack()
		st='uname -a'
		val=os.popen(st,'r').readlines()
		sval=string.split(val[0])

		desired_modules,desired_font_sizes,enable_saslock = user_parameters()

		if(enable_saslock == 'YES' or enable_saslock == 'yes'):
			if(os.path.isfile('.saslock')): 
				message_warning =  "\nWARNING: ANOTHER INSTANCE OF SASSIE MAY BE RUNNING IN THIS DIRECTORY"
				print message_warning ; print message_warning ; print message_warning+'\n' ; print os.getcwd()+'\n'
				title = 'SASLOCK FILE FOUND IN CURRENT DIRECTORY'
				txt='\n\nPress "OK" to delete saslock file and continue \nPress cancel to quit now\n'
				ans=Message(title=message_warning,message=title+txt,default='cancel',icon=QUESTION,type='okcancel').show()
				self.update_idletasks()
				if(str(ans)=="ok"):		
					print 'removing old saslock file from current directory'
					rmst='rm -f .saslock'
					os.system(rmst)
					write_saslock()
				elif(str(ans)=="cancel"):		
					self.quitnow()
			else:
				write_saslock()


		license="""
    SASSIE  Copyright (C) 2011 Joseph E. Curtis, Ph.D.
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software. You are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.

    Citation: Comp. Phys. Comm. Volume 183, Issue 2, Pages 382-389 (February 2012).
	"""
		print license

		if(desired_font_sizes['default_font'] != 'YES' and desired_font_sizes['default_font'] != 'yes'):
			button_font = desired_font_sizes['button_font']
			local_font = desired_font_sizes['local_font']
			button_pady = desired_font_sizes['button_pady']
		else:
			root.customFont = tkFont.Font(family = "TkCaptionFont")
			root.customFont.configure(size=12)

			button_font = root.customFont
			local_font = root.customFont
			#button_font = ('Helvetica', 12,'bold')
			local_font = ('Helvetica', 12,'bold')
			button_pady = 8

		self.option_add('*Font', local_font)

		vers='version 0.99_rev_1072 : 11/22/12'
		Label(self,text='---------------',fg='black',bg='white',font=local_font).pack(side=TOP,fill=BOTH)
		Label(self,text='SASSIE',fg='black',bg='white',font=local_font).pack(side=TOP,fill=BOTH)
		Label(self,text=vers,fg='black',bg='white',font=local_font).pack(side=TOP,fill=BOTH)
		Label(self,text='---------------',fg='black',bg='white',font=local_font).pack(side=TOP,fill=BOTH)

		numitems=len(aprogs)
		
		for i in range(numitems):
			lkey=aprogs[i][0]
			if(aprogs[i][1]!='section_header'):
				lvalue=aprogs[i][1]
				if(lkey=='Hydropro' and ('Darwin' in sval)):
					pass
				if(desired_modules[lkey] != 'YES' and desired_modules[lkey] != 'yes'):
					pass
				else:
					this_button=Button(self,text=lkey,command=lvalue,font=button_font,relief=RIDGE,padx=8,pady=button_pady,bg='lightgrey',foreground='black',takefocus='',activeforeground='black',border=6,highlightbackground='lightgrey').pack(side=TOP,fill=BOTH,padx=10,pady=2)
					#this_button=Button(self,text=lkey,command=self.start_thread(lvalue),font=button_font,relief=RIDGE,padx=8,pady=button_pady,bg='lightgrey',foreground='black',takefocus='',activeforeground='black',border=6,highlightbackground='lightgrey').pack(side=TOP,fill=BOTH,padx=10,pady=2)
			else:
				Label(self,text=aprogs[i][0],bg='white',fg='#0198E1',font=local_font).pack()
				
		Label(self,text=' ',fg='black',bg='white').pack()
		Button(self,text="Quit",command=self.quit,font=button_font,relief=RIDGE,padx=12,pady=button_pady,bg='white',fg='red',border=4).pack(side=TOP,fill=BOTH,padx=10)
		ttxt=time.ctime()
		Label(self,text=ttxt,fg='black',bg='white',font=local_font).pack(side=BOTTOM,fill=BOTH)
		self.update()


		root.bind('<Control-s>',lambda event, task=root: self.minus(task))
		root.bind('<Control-b>',lambda event, task=root: self.plus(task))
		root.bind('<Control-q>',lambda event, task=root: self.key_quit(task))

	def start_thread(self,lvalue):

		process=multiprocessing.Process(target=lvalue)

	def quitnow(self):
		self.master.destroy()

	def quit(self):
		os.system('rm -f .saslock')
		self.master.destroy()

	def minus(self,root):
		size = root.customFont['size']
		if(size > 5):
			root.customFont.configure(size=size-2)
		root.update_idletasks()
	
	def plus(self,root):
		size = root.customFont['size']
		root.customFont.configure(size=size+2)
		root.update_idletasks()
	
	def key_quit(self,root):
		os.system('rm -f .saslock')
		root.destroy()


def OnExit():
	os.system('rm -f .saslock')
	root.destroy()
	
if __name__=='__main__': 

	root = Tk()
	app = Prog(root)
	root.protocol("WM_DELETE_WINDOW",OnExit)
	root.mainloop()
