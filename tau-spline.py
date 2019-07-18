#! /usr/bin/python3

# by Silvia Sabotta (silvia@tls-tautenburg.de)

from astropy.io import fits
from tkinter import ttk
from pyraf import iraf, epar
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk

import os
import shutil
import cosmics


class Calibration_Frame(tk.Frame):
    def __init__(self, master=None, calibration_type="Flat"):  
        tk.Frame.__init__(self, master)
        self.LabelUpload = None
            
        self.__init_tabs__(calibration_type)
              
        for i in range(7):
            self.rowconfigure(i, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        
    def __init_tabs__(self, calibration_type):
        if calibration_type == "Flat":
            self.calibration_flat()
        elif calibration_type == "ThAr":
            self.calibration_thar()
        elif calibration_type == "Reduce":
            self.calibration_reduce()
        elif calibration_type == "Prepare":
            self.calibration_prepare()
        elif calibration_type == "Additional":
            self.calibration_additional()
        else:
            print("Something very weird just happened o.O...")       
        
    def get_newname(self, e1):
        self.masterflatname = e1.get()
        print(self.masterflatname)
        
    def call_eparam(self,function):
        toplevel = tk.Toplevel(self)
        toplevel.grid()
        epar.epar(function, parent=toplevel, isChild=1)
        toplevel.withdraw()
    
    def showfiles(self, filelist_new):
     
        filelist_readable = "1. " + filelist_new[0] + "\n" 
        for i in range(1,len(filelist_new)):
            nr = i+1
            filelist_readable = filelist_readable + str(nr) + ". " + filelist_new[i] + "\n"
        toplevel = tk.Toplevel(self)
        toplevel.grid()
        toplevel.configure(bg="white")        
        label1 = tk.Label(toplevel, text=filelist_readable)
        label1.configure(bg="white")
        label1.grid() 
        if len(filelist_new) > 50:
            print(len(filelist_new))
            print(filelist_readable)

    def create_label(self, text, rownr, colnr):
        self.LabelUpload = tk.Label(self)
        self.LabelUpload.grid(row=rownr,column=colnr,sticky='nwse' )
        self.LabelUpload.config(bg = "white")
        self.LabelUpload["text"] = text

    def exit_function(self):
        toplevel = tk.Toplevel(self)
        frame = tk.Frame(toplevel)
        frame.grid(row=0,column=0,sticky='nwse')
        label1 = tk.Label(frame, text="Do you really want to exit?")
        label1.grid(row=0,column=0,sticky='nwse')
        button = tk.Button(frame,  text="Yes",  fg="red",  command=root.destroy)
        button.grid(row=1,column=0,sticky='nwse')
        slogan = tk.Button(frame,text="No", command=toplevel.destroy)
        slogan.grid(row=2,column=0,sticky='nwse')
        
    def create_button(self, text, command, rownr, colnr):

        self.loadfile = tk.Button(self)
        self.loadfile["bg"]   = "white"
        self.loadfile["text"] = text
        self.loadfile["command"] = command
        self.loadfile.grid(row=rownr,column=colnr,sticky='nwse')

   
    def calibration_flat(self):
        self.configure(bg="white")
        self.create_label("Step 1: Upload Flat frames",0,0)
        self.flat_files = []
        self.create_button("Upload Files", lambda: self.load_files(self.flat_files), 1, 0) 
        self.create_button("Show Files", lambda: self.showfiles(self.flat_files), 1, 1)

        self.create_label("Step 2: Choose new Flat frame name",2,0)
        self.masterflatname = "Masterflat.fits"
        e1 = tk.Entry(self)
        e1.grid(row = 3, column= 0, sticky='nwse')
        self.create_button("Send", lambda: self.get_newname(e1), 3,1)

        apflatten_param = {"interac": "no",
                           "find": "no",
                           "recente": "yes",
                           "resize": "no",
                           "edit": "no",
                           "fittrac": "no",
                           "flatten": "yes",
                           "fitspec": "no",
                           "function": "legendre",
                           "order": 40}
        set_iraf_apflatten(apflatten_param)
             
        self.create_label("Step 3: Create Flat frame",4,0)
        self.create_button("Review Parameters for imcombine",\
        lambda: self.call_eparam(iraf.imcombine), 5, 0)
        self.create_button("Review Parameters for apflatten", \
        lambda: self.call_eparam(iraf.apflatten), 5, 1)
        
        self.create_button("GO", lambda: produce_flat(self.flat_files, self.masterflatname), 6, 0)
        self.create_button("Quit", self.exit_function, 6,1)
       
    def calibration_thar(self):
        self.configure(bg="white")
        self.create_label("Step 1: Upload ThAr frames",0,0)
        self.thar_files = []
        self.create_button("Upload Files", lambda: self.load_files(self.thar_files), 1, 0) 
        self.create_button("Show Files", lambda: self.showfiles(self.thar_files), 1, 1)       
        
        apall_param = {"interactive": "no",
                            "find": "no",
                            "recenter": "yes",
                            "trace": "no",
                            "resize": "no",
                            "edit": "no",
                            "fittrac": "no",
                            "extract": "yes",
                            "extras": "no",
                            "review": "no"
                            }     
        set_iraf_apall(apall_param)
                            
        self.thar_reference = "ThArMaster.fits"

        self.create_label("Step 2: Upload different Thar Reference",2,0)
        self.create_button("Upload File", lambda: self.load_files([self.thar_reference]), 3, 0) 

        self.create_label("Step 3: Create ThAr frames",4,0)
        self.create_button("Review Parameters for ecreidentify",\
        lambda: self.call_eparam(iraf.ecreidentify), 5, 0)
        
        self.create_button("GO", lambda: produce_thar(self.thar_files, self.thar_reference), 6, 0)
        self.create_button("Quit", self.exit_function, 6,1)
                 
    def calibration_reduce(self):
        self.configure(bg="white")
        self.create_label("Step 1: Upload Science frames",0,0)
    
        self.calibration_files = []
        self.create_button("Upload Files", lambda: self.load_files(self.calibration_files), 1, 0) 
        self.create_button("Show Files",  lambda: self.showfiles(self.calibration_files), 1, 1)       

        self.create_label("Step 2: Upload calibration frames",2,0)
        self.masterflat_files = []
        self.create_button("Load Flat Frames", lambda: self.load_files(self.masterflat_files), 3, 0) 
        self.masterthar_files = []
        self.create_button("Load ThAr Frames", lambda: self.load_files(self.masterthar_files), 3, 1)    
 
        self.create_label("Step 3: Reduce Data",4,0)         
       
        iraf.fit1d.setParam("interac","no")

        apall_param = {"interactive": "no",
                            "find": "no",
                            "recenter": "yes",
                            "trace": "no",
                            "resize": "no",
                            "edit": "no",
                            "fittrac": "no",
                            "extract": "yes",
                            "extras": "no",
                            "review": "no"}     
        set_iraf_apall(apall_param)

        apscatter_param = {"interactive": "yes",
                           "find": "no",
                           "recenter": "yes",
                           "trace": "no",
                           "resize": "no",
                           "edit": "no",
                           "fittrac": "no",
                           "subtract": "yes",
                           "smooth": "yes",
                           "fitscatter": "yes",
                           "fitsmooth": "yes"}

        apscat1_param = { "high_reject" : "1",
                          "order" : "15"}
        apscat2_param = {"order" : "40"}
        
        set_iraf_apscatter(apscatter_param)
        set_iraf_apscat1(apscat1_param)
        set_iraf_apscat2(apscat2_param)
                        
        self.create_button("Review Parameters for apscatter",\
        lambda: self.call_eparam(iraf.apscatter), 5, 0)
        self.create_button("Review Parameters for apall", \
        lambda: self.call_eparam(iraf.apall), 5, 1)
        self.create_button("Review Parameters for fit1d", \
        lambda: self.call_eparam(iraf.fit1d), 4, 1)
        
        self.create_button("GO", lambda: reduce_data(self.calibration_files, self.masterflat_files, self.masterthar_files), 6, 0)
        self.create_button("Quit", self.exit_function, 6,1)

    def calibration_prepare(self):
        self.configure(bg="white")
        self.create_label("Step 1: Upload files to rename",0,0)
        self.prepare_files = []
        self.create_button("Upload Files", lambda: self.load_files(self.prepare_files), 1, 0) 
        self.create_button("Show Files", lambda: self.showfiles(self.prepare_files), 1, 1)       
        self.create_button("Copy and Give Object Names", lambda: rename(self.prepare_files), 2, 0)
        self.create_button("Make Lists", lambda: makelists(self.prepare_files), 2, 1)
        self.create_button("Quit", self.exit_function, 3,1)

    def calibration_additional(self):
        self.configure(bg="white")
        self.create_label("Step 1: Upload files for additional tasks",0,0)
        self.additional_files = []
        self.create_button("Upload Files", lambda: self.load_files(self.additional_files), 1, 0) 
        self.create_button("Show Files", lambda: self.showfiles(self.additional_files), 1, 1)       
        self.create_button("Apply Barycentric correction", lambda: barycor(self.additional_files), 2, 1)
        self.create_button("Edit header and write Molecfit Parameter file", lambda: headcor(self.additional_files), 2, 0)
        self.create_button("Calculate and Plot RVs", lambda: rv_plot(self.additional_files), 3, 0)
        self.create_button("Just Plot", lambda: justplot(self.additional_files), 3, 1)
        self.create_button("Just Normalize", lambda: normalize_and_merge(self.additional_files), 4, 0)
        self.create_button("Quit", self.exit_function, 4,1)
        
    def load_files(self, file_list):   
        del file_list[:]
        filez = tk.filedialog.askopenfilenames(parent=root, title='Choose a file')
        file_list_dummy = self.tk.splitlist(filez)      
        
        if len(file_list_dummy) == 0:
            return
        
        
        if "lis" in file_list_dummy[0]:
            path = os.path.dirname(file_list_dummy[0])
            with open(file_list_dummy[0]) as fobj:
                for line in fobj:
                    single_file = line.rstrip()
                    file_list.append(os.path.join(path, single_file))

        else:
            for data_file in file_list_dummy:
                file_list.append(data_file)
            
def barycor(filelist_new):
    iraf.reset(obsdb = 'home$obsdb.dat')
    for i in range(len(filelist_new)):
        hdulist = fits.open(filelist_new[i])
        header_time_of_observation = hdulist[0].header['DATE-OBS']
        year_of_observation = int(header_time_of_observation[:4])
        month_of_observation = int(header_time_of_observation[5:7])
        day_of_observation = int(header_time_of_observation[8:10])
        right_ascension = hdulist[0].header['RA']
        declination = hdulist[0].header['DEC']
                        
        try: 
            ut_of_observation = hdulist[0].header['UT']
        except KeyError:
            ut_of_observation = int(header_time_of_observation[11:13]) + int(header_time_of_observation[14:16])/60 + int(header_time_of_observation[17:19])/3600
        exposure_time = hdulist[0].header['EXP_TIME']
    
        output_filename_dummy = filelist_new[i].replace("norm.","norm.dummy.")
        output_filename_dummyI = output_filename_dummy.replace("norm-1","norm-1.dummy.")
        output_filename = output_filename_dummyI.replace("merged","merged.dummy.")

        iraf.scopy(filelist_new[i], output_filename)
        iraf.hedit(images = output_filename, fields = "UT", value = ut_of_observation)
        iraf.hedit(images = output_filename, fields = "EPOCH", value = "2000")
        iraf.hedit(images = output_filename, fields = "EXP-TIME", value = exposure_time)
        iraf.rvcorrect(images = output_filename, year = year_of_observation, month= month_of_observation, day = day_of_observation, ut = ut_of_observation, ra= right_ascension, dec=declination)
        
        output_filename_final = output_filename.replace("dummy.", "rvcorrect.")
        print(output_filename, output_filename_final)
        iraf.dopcor(output_filename, output_filename_final , redshift= "-VHELIO", isvelocity = "yes")
        
        
        os.remove(output_filename)
        hdulist.close()

def headcor(filelist_new):
    bash_script = "/home/silvia/molecfit.sh"

    with open(bash_script, "w") as bashscript:

        for i in range(len(filelist_new)):
            hdulist = fits.open(filelist_new[i])
            seconds = float(hdulist[0].header['TM-START'])
            date_start = hdulist[0].header['JD-START'] - 2450000
            airmass = hdulist[0].header['AIRMASS']
            altitude_angle = np.arcsin(1/airmass)/(2*np.pi)*360.0
    
            right_ascension = hdulist[0].header['RA']
            declination = hdulist[0].header['DEC']
            
            iraf.hedit(images = filelist_new[i], fields = "UTCS", value = seconds)
            iraf.hedit(images = filelist_new[i], fields = "HEIGHT_ANGLE", value = altitude_angle)
            iraf.hedit(images = filelist_new[i], fields = "MJD", value = date_start)
            
            hdulist.close()
    
            parameterfile = "/home/silvia/PhD/A_star_templatefile.par"
            
            basedir = os.path.dirname(filelist_new[0])
        
            filename_dummy = filelist_new[i].replace("fits", "par")
            parfilename = os.path.join(basedir, "molecfit_" + filename_dummy.replace(os.path.join(basedir,""), ""))
            
            with open(parfilename, "w") as parfile:
                with open(parameterfile, "r") as defaultfile:
                    filepath = filelist_new[i]
                    for line in defaultfile:
                        newline = line.replace("REPLACE_THIS", filepath)
                        newline2 = newline.replace("OUTPUT", basedir)

                        parfile.write(newline2)
            bashscript.write("Molecfit/bin/molecfit " +str(parfilename) + "\nMolecfit/bin/calctrans \
" + str(parfilename) + "\n")



def rv_plot(filelist_new):
    rv_directory = os.path.dirname(filelist_new[0])
    os.chdir(rv_directory)
    
    template = "telutemplate.fits"
    basedir = os.path.dirname(filelist_new[0])
    templatefile = os.path.join(basedir, template)

    if os.path.exists("telushift.txt"):
        os.remove("telushift.txt")
    if os.path.exists("rvshift.txt"):
        os.remove("rvshift.txt")      
        
        
    iraf.fxcor(objects=os.path.join(basedir,"*merged.fits"), templates = templatefile, continuum = "both", filter="none", rebin="smallest", \
    function = "gaussian", width = "10.0", height="0.0", minwidth ="3.0", maxwidth = "20.0", window = "200.0", background= "0.2", \
    output = "telushift", verbose = "txtonly", osample = "5880:5975", rsample = "5880:5975", interactive = "no")
    
    rvs = []
    time = []
    with open("telushift.txt", "r") as defaultfile:
        for line in defaultfile:
            if "#" in line:
                pass
            else:
                rvs.append(line[-30:-21])
    
    for i in range(len(filelist_new)):
        hdulist = fits.open(filelist_new[i])
        time.append(hdulist[0].header['JD-START'])


    plt.plot(time, rvs, '.')

    template = "mascara_megatemplate.fits"
    templatefile = os.path.join(basedir, template)

### For A Type Stars
    iraf.fxcor(objects=os.path.join(basedir,"*rvcorrect_TAC.fits"), templates = templatefile, continuum = "both", filter="template", rebin="smallest", \
    function = "center1d", width = "200.0", height="0.0", minwidth ="3.0", maxwidth = "200.0", window = "1000.0", \
    output = "rvshift", verbose = "txtonly", osample = "4610:4830, 4893:5679, 6019:6524, 6598:6860", rsample = "4610:4830, 4893:5679, 6019:6524, 6598:6860", interactive = "no")

### For G Type Stars
#    iraf.fxcor(objects=os.path.join(basedir,"*rvcorrect_TAC.fits"), templates = templatefile, continuum = "both", filter="both", rebin="smallest", \
#    function = "gaussian", width = "20.0", height="0.2", minwidth ="10.0", maxwidth = "200.0", window = "200.0", \
#    output = "rvshift", verbose = "txtonly", osample = "4900:6450", rsample = "4900:6450", interactive = "no")              




    rvs2 = []
    time2 = []
    with open("rvshift.txt", "r") as defaultfile:
        for line in defaultfile:
            if "#" in line:
                pass
            else:
                rvs2.append(line[-30:-21])

#    iraf.fxcor(objects=os.path.join(basedir,"*rvcorrect_TAC.fits"), templates = templatefile, continuum = "both", filter="both", rebin="smallest", \
#    function = "gaussian", width = "20.0", height="0.2", minwidth ="10.0", maxwidth = "200.0", window = "200.0", \
#    output = "rvshift", verbose = "txtonly", osample = "4900:6450", rsample = "4900:6450", interactive = "no")

 #   rvs3 = []
  #  with open("rvshift.txt", "r") as defaultfile:
   #     for line in defaultfile:
    #        if "#" in line:
     #           pass
      #      else:
       #         rvs3.append(line[-30:-21])
    
    for i in range(len(filelist_new)):
        hdulist = fits.open(filelist_new[i])
        time2.append(hdulist[0].header['JD-START'])
 #       iraf.dopcor(filelist_new[i], filelist_new[i].replace("merged_TAC","starvelocity"), redshift= rvs2[i], isvelocity = "yes")

    with open(os.path.join(basedir,"rvshift.dat"), "w") as rvfile:
        with open(os.path.join(basedir,"telushift.dat"), "w") as telufile:
            with open(os.path.join(basedir,"rvtelushift.dat"), "w") as rvtelufile:
                for i in range(len(rvs)):
                    rvfile.write(str(time[i]) + "," + rvs2[i] + "\n")
                    telufile.write(str(time[i]) + "," + rvs[i] + "\n")
                    final_shift = float(rvs2[i]) + float(rvs[i])
                    rvtelufile.write(str(time[i]) + "," + str(final_shift) + "\n")
    
    plt.plot(time2, rvs2, '.')
    plt.show() 

def justplot(filelist_new):
    plot_directory = os.path.dirname(filelist_new[0])
    os.chdir(plot_directory)
    
    for i in range(len(filelist_new)):
        xvalues = np.loadtxt(filelist_new[i],delimiter=",",  usecols=0)
        yvalues = np.loadtxt(filelist_new[i],delimiter=",",  usecols=1)
        plt.plot(xvalues, yvalues, 'o', label= filelist_new[i])
        print(np.mean(yvalues), np.std(yvalues))
     
    plt.title("RV-Standard Shows Long Term Stability of 170 m/s")
    plt.xlabel("Julian Date")
    plt.ylabel("Relative RV")
 #   plt.legend()
    plt.show() 


def rename(filelist_new):    
    for i in range(len(filelist_new)):
        hdulist = fits.open(filelist_new[i])
        objname = hdulist[0].header['OBJECT']
        hdulist.close()
        
        old = filelist_new[i]
     
        path = os.path.dirname(filelist_new[i]) 
        new = filelist_new[i].replace(os.path.join(path,""), os.path.join(path , objname) + "_")
        new2 = new.replace(".fits", ".raw.fits")
        
        iraf.imcopy(old,new2)

        
        path = os.path.dirname(filelist_new[i]) 
        FileName = os.path.join(path , objname) + ".lis"
        filelist = open(FileName,"a")
        filelist.write(new2 + "\n")
        filelist.close()
       
def makelists(filelist_new):
    for i in range(len(filelist_new)):
        hdulist = fits.open(filelist_new[i])
        objname = hdulist[0].header['OBJECT']
        hdulist.close()
        
        path = os.path.dirname(filelist_new[i]) 
        FileName = os.path.join(path , objname) + ".lis"
        filelist = open(FileName,"a")
        filelist.write(filelist_new[i] + "\n")
        filelist.close()

   
def produce_flat(filelist_flat, masterflat_name):
    flat_directory = os.path.dirname(filelist_flat[0])
    
    filelist_text = open("filelist.lis","w")
    for i in range(len(filelist_flat)):
        filelist_text.write(filelist_flat[i] + "\n")

    filelist_text.close()
    
    if len(filelist_flat) > 1:
        iraf.imcombine("@filelist.lis","Flat.raw.fits", combine = "average", reject = "minmax", nlow = 1, nhigh=1, nkeep = 1)
        trim = trim_remove_bias("Flat.raw.fits")    
        iraf.apflatten(trim,masterflat_name, referen = "find_orders_flat", apertures = "3-50", trace = "no", gain=0.1913, readnoise = 3.380)  
    else:
        iraf.apflatten(filelist_flat[0], masterflat_name, referen = "find_orders_flat", apertures = "3-50", trace = "no", gain=0.1913, readnoise = 3.380)    
    
    os.remove("Flat.raw.fits")
    os.remove(trim)
    os.remove("filelist.lis")
    shutil.move(masterflat_name, os.path.join(flat_directory))

    print("Success! Produced file:" + masterflat_name)

def write_reftable(filelist_new, filelist_masterthars):
    reftable_directory = os.path.dirname(filelist_new[0])
    os.chdir(reftable_directory)

    remove_file(reftable_directory,"reftable.dat")
    remove_file(reftable_directory,"science.lis")
    remove_file(reftable_directory,"science_final.lis")
    
    reftable = open("reftable.dat","w")
    science_filelist = open("science.lis","w")
    science_finallist = open("science_final.lis","w")
    
    for j in range(len(filelist_new)):
        hdulist1 = fits.open(filelist_new[j])
    
        date_science = hdulist1[0].header['JD-START']
    
        science_day = int(date_science)
        science_hour = date_science%1
    
        hdulist1.close()
    
        right_one = find_thar(filelist_masterthars, science_day,science_hour)

        try:
            reftable.write(filelist_new[j].replace(".fits",".extracted.fits") + "  " + right_one + "\n")
            science_filelist.write(filelist_new[j].replace(".fits",".extracted.fits")  + "\n")
            science_finallist.write(filelist_new[j].replace(".extracted.fits",".science.fits") + "\n")
        except TypeError:
            print("Cannot find ThAr file for" + filelist_new[j].replace(".fits",".extracted.fits"))
    reftable.close()
    science_filelist.close()
    science_finallist.close() 
    
    print("Success! Produced file: reftable.dat, science.lis, science_final.lis")
    
    
def produce_thar(filelist_new, thar_reference):
    current_directory = os.getcwd()
    thar_directory = os.path.dirname(filelist_new[0])
  
    filelistText = open(os.path.join(thar_directory,"ThAr.extracted.lis"),"w")
    
    for i in range(len(filelist_new)):
        
        trim = trim_remove_bias(filelist_new[i])    
        extracted = trim.replace( ".trim_will_be_removed.fits", ".extracted.fits")        
        
        iraf.apall(trim, output=extracted, references= "find_orders.fits", profiles= "find_orders.fits", apertures = "1-51" )
        os.remove(trim)
        
        extracted_without_path = extracted.replace(os.path.join(thar_directory,""),"")
        
        filelistText.write(extracted_without_path + "\n")

    filelistText.close()
    os.chdir(thar_directory)

    thar_reference_short = thar_reference.replace(".fits","")
    
    iraf.ecreidentify("@ThAr.extracted.lis", reference = thar_reference_short, shift="INDEF")

    os.chdir(current_directory)
    
    print("Success! Produced ThAr reference files.")
       

def find_thar(filelist_masterthars, day,hour):
    for i in range(len(filelist_masterthars)):
        hdulist2 = fits.open(filelist_masterthars[i])
   
        date_thar = hdulist2[0].header['JD-START']
        
        thar_day = int(date_thar)
        thar_hour = date_thar%1 
        
  #      print a2, thar_day, thar_hour, a, science_day, science_hour
        hdulist2.close()        
        if (thar_day == day and hour > thar_hour):
            
            
            directory_name = os.path.dirname(filelist_masterthars[i])
            filename_without_path= filelist_masterthars[i].replace(os.path.join(directory_name,""),"")
 
            return filename_without_path          

def reduce_data(filelist_new, filelist_masterflats,  filelist_masterthars):
    if len(filelist_new) == 0:
        print("Please, choose science files to reduce first!")
    elif len(filelist_masterflats) == 0:
        print("Please, choose masterflat files first!")
    elif len(filelist_masterthars) == 0:
        print("Please, choose Thorium-Argon reference files first!")
    else:
        print("Starting reduction process...")
        extracted_science_files =[]
        current_directory = os.getcwd()
      
        years = ['{:04d}'.format(year) for year in range(2005, date.today().year + 1)]
        month = ['{:02d}'.format(mon) for mon in range(1, 13)]    
        
        for i in range(len(filelist_new)):
            trim = trim_remove_bias(filelist_new[i])
            chosen_masterflat = None
    
            if os.path.exists(os.path.join(current_directory, "dummyI.fits")):
                os.remove(os.path.join(current_directory, "dummyI.fits"))
            if os.path.exists(os.path.join(current_directory, "dummyII.fits")):
                os.remove(os.path.join(current_directory, "dummyII.fits"))
            
            for y in years:
                for mon in month:
                    date_check_1 = y + "_" + mon
                    date_check_2 = y + "-" + mon
                    if ((date_check_1 in trim) or (date_check_2 in trim)) == True:
                        for masterflat in filelist_masterflats:
                            if ((date_check_1 in masterflat) or (date_check_2 in masterflat)) == True:
                                chosen_masterflat = masterflat
                                
            if chosen_masterflat is None :
                print("Masterflat File not found!!")
                continue
                             
            iraf.imarith(trim, "/", chosen_masterflat, "dummyI.fits")
        
            iraf.apscatter("dummyI.fits", output ="dummyII.fits", interactive="No", apertures = "1-51", references= "find_orders.fits")    
        
            clean = trim.replace(".trim_will_be_removed.fits",".clean.fits")
            mask = trim.replace(".trim_will_be_removed.fits",".mask.fits")
                
            mean_image = iraf.imstat("dummyII.fits", Stdout=1, fields="mean", format="no")
            print(mean_image)
    
            if float(mean_image[0]) > 5000:                
                array, header = cosmics.fromfits("dummyII.fits")
            
                sigclip = 50.0
                objlim = 50.0
                
                c = cosmics.cosmicsimage(array, gain=0.368,  readnoise=3.7, sigclip=sigclip, sigfrac=10.0, objlim=objlim)
    
                c.run(maxiter = 0)
                cosmics.tofits(clean, c.cleanarray, header)
                cosmics.tofits(mask, c.mask, header)
                
            elif 5000 > float(mean_image[0]) > 1000:
                array, header = cosmics.fromfits("dummyII.fits")
                sigclip =20.0
                objlim = 20.0
                c = cosmics.cosmicsimage(array, gain=0.368,  readnoise=3.7, sigclip=sigclip, sigfrac=4.0, objlim=objlim)

                c.run(maxiter = 2)
                cosmics.tofits(clean, c.cleanarray, header)
                cosmics.tofits(mask, c.mask, header)   
            else:
                array, header = cosmics.fromfits("dummyII.fits")
                     
                sigclip = 5.0
                objlim = 5.0      
                
                c = cosmics.cosmicsimage(array, gain=0.368,  readnoise=3.7, sigclip=sigclip, sigfrac=1.8, objlim=objlim)
    
                c.run(maxiter = 2)
                cosmics.tofits(clean, c.cleanarray, header)
                cosmics.tofits(mask, c.mask, header)
            
            extract = trim.replace(".trim_will_be_removed.fits",".extracted.fits")
            extracted_science_files.append(extract)        
            
            iraf.apall(clean, output=extract, references= "find_orders", profiles= "find_orders", apertures = "1-51" )
        #    iraf.apall(clean, output=extract, references= "find_orders", profiles= "find_orders", apertures = "3-50" )
        #    iraf.apall(clean, output=extract, references= "find_orders", profiles= "find_orders", apertures = "3-50" )
            os.remove("dummyI.fits")
            os.remove("dummyII.fits")
            os.remove(trim)
            
        write_reftable(filelist_new, filelist_masterthars)    
     
        thar_directory = os.path.dirname(filelist_masterthars[0])
        
        if thar_directory == os.path.dirname(filelist_new[0]):
            pass
        else :
            if os.path.exists(os.path.join(thar_directory, "science.lis")):
                os.remove(os.path.join(thar_directory, "science.lis"))
                
            if os.path.exists(os.path.join(thar_directory, "reftable.dat")):
                os.remove(os.path.join(thar_directory, "reftable.dat"))
        
            shutil.move("science.lis", thar_directory)
            shutil.move("reftable.dat", thar_directory)
        os.chdir(thar_directory)
    
        iraf.refspectra("@science.lis", references="reftable.dat", ignoreaps = "Yes", sort = "", group = "", override = "yes", confirm= "no", assign = "yes")
        iraf.dispcor("@science.lis", output = "@science.lis")

        if thar_directory == os.path.dirname(filelist_new[0]):
            pass
        else:
            remove_file(current_directory, "science.lis")
            remove_file(current_directory, "reftable.dat")
                                                                                               
            shutil.move("science.lis", current_directory)
            shutil.move("reftable.dat", current_directory)
        os.chdir(current_directory)
    
        normalize_and_merge(extracted_science_files)

        print("Success! Reduction process complete.")
    
def set_iraf_apall(params):
    """
    Sets new parameters to iraf.apall.
    
    :param params: 
        dict with the keys as the parameters names and 
        the items as new parameters
    :type params: dict
    """
    for key in params:
        iraf.apall.setParam(key, params[key])


def set_iraf_apscatter(params):
    """
    Sets new parameters to iraf.apscatter.
    
    :param params: 
        dict with the keys as the parameters names and 
        the items as new parameters
    :type params: dict
    """
    for key in params:
        iraf.apscatter.setParam(key, params[key])



def set_iraf_apscat1(params):
    """ 
    Sets new parameters to iraf.apscat1.                                                             
    :param params:                                                                                   
        dict with the keys as the parameters names and 
        the items as new parameters                                                                 
    :type params: dict                                                                               
    """
    for key in params:
        iraf.apscat1.setParam(key, params[key])

def set_iraf_apscat2(params):
    """                                                                  
    Sets new parameters to iraf.apscat2.                                                           \
                                                                                                     
    :param params:                                                                                   
        dict with the keys as the parameters names and                                               
        the items as new parameters                                                                  
    :type params: dict                                                                               
    """
    for key in params:
        iraf.apscat2.setParam(key, params[key])
                                                

def set_iraf_apflatten(params):
    """
    Sets new parameters to iraf.apflatten.
    
    :param params: 
        dict with the keys as the parameters names and 
        the items as new parameters
    :type params: dict
    """
    for key in params:
        iraf.apflatten.setParam(key, params[key])


def remove_file(path, name):
    """
    Removes a file if the files exists.

    :param path: The path to the directory
    :type path: str
    :param name: The name of the file
    :type name: str
    """
    path = os.path.join(path, name)
    if os.path.exists(path):
        os.remove(path)      

def trim_remove_bias(WhatToTrim):
    imstat = WhatToTrim + "[2098:2147,*]"
    mean_overscan = iraf.imstat(imstat, Stdout=1, fields="mean", format="no")
        
    single_file_new =  WhatToTrim.replace(".fits",".-overscan_will_be_removed.fits")
    iraf.imarith(WhatToTrim , "-", mean_overscan[0] , single_file_new)    

    old = single_file_new + "[51:2097,3:2063]"
    new = WhatToTrim.replace(".fits",".trim_will_be_removed.fits")
    
    iraf.imcopy(old,new)

    os.remove(single_file_new)
    
    return new

def place_here(filename):
    current_directory = os.getcwd()
    new_name = os.path.join(current_directory, filename)
    return new_name
    

def normalize_and_merge(reduced_science_files):
    current_directory = os.getcwd()
   
    for k in range(len(reduced_science_files)):
        
        remove_file(current_directory, "norm.dummyI.fits")
        remove_file(current_directory,"norm.dummyIa.fits")
        remove_file(current_directory, "norm.dummyII.fits")
        remove_file(current_directory, "norm.dummyIIa.fits")
        remove_file(current_directory, "norm.dummyIII.fits")
        remove_file(current_directory, "norm.dummyIV.fits")
        remove_file(current_directory, "norm.dummyV.fits")
        remove_file(current_directory, "norm.dummy1.fits")
        remove_file(current_directory, "norm.dummy2.fits")
        remove_file(current_directory, "norm.dummy3.fits")
        remove_file(current_directory, "norm.dummy4.fits")
        remove_file(current_directory, "norm.dummy5.fits")
        remove_file(current_directory, "norm.dummy6.fits")
        remove_file(current_directory, "norm.dummy7.fits")
        remove_file(current_directory, "norm.dummy8.fits")
        
        ranges = np.loadtxt("ranges.lis", delimiter=",")
        aperturlist = open("test.norm.apertures.lis", "w")

        merge = reduced_science_files[k].replace(".extracted.fits",".norm.fits")   

        iraf.scopy(reduced_science_files[k], place_here("norm.dummyI.fits"), apertures ="5:14" ,format = "multispec")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummyII.fits"), apertures ="35:42" ,format = "multispec")

        iraf.fit1d(place_here("norm.dummyI.fits"),  place_here("norm.dummyIa.fits"), naverage = 1, axis= 2, type = "fit", low_rej=1.0,  high_rej = 2.0, order = 2, niterate= 2, func= "spline3", sample = "*")
        iraf.fit1d(place_here("norm.dummyII.fits"), place_here("norm.dummyIIa.fits"), naverage = 1, axis= 2, type = "fit", low_rej=1.0,  high_rej = 2.0, order = 2, niterate= 2, func= "spline3", sample = "*")

        iraf.scopy(reduced_science_files[k], place_here("norm.dummy1.fits"), apertures ="1" ,format = "multispec", w1 = "4544")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummy2.fits"), apertures ="2" ,format = "multispec", w1 = "4569")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummy3.fits"), apertures ="3:7" ,format = "multispec")
        iraf.scopy(place_here("norm.dummyIa.fits"), place_here("norm.dummy4.fits"), apertures ="8:10" ,format = "multispec")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummy5.fits"), apertures ="11:37" ,format = "multispec")
        iraf.scopy(place_here("norm.dummyIIa.fits"), place_here("norm.dummy6.fits"), apertures ="38:40" ,format = "multispec")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummy7.fits"), apertures ="41:50" ,format = "multispec")
        iraf.scopy(reduced_science_files[k], place_here("norm.dummy8.fits"), apertures ="51" ,format = "multispec", w2 = "7591")        
        
        iraf.scombine("@normalization.list.lis", place_here("norm.dummyIII.fits"), group = 'apertures')        
        
        iraf.fit1d(place_here("norm.dummyIII.fits"), place_here("norm.dummyIV.fits"), naverage = 1, axis= 1, type = "fit", low_rej=0.8,  high_rej = 2.0, order = 7, niterate= 4, func= "spline3", sample = "*")

        iraf.sarith(reduced_science_files[k], "/", place_here("norm.dummyIV.fits"), place_here("norm.dummyV.fits"))        
        
        iraf.fit1d(place_here("norm.dummyV.fits"), merge, naverage = 1, axis= 1, type = "ratio", low_rej=0.2,  high_rej = 2.0, order = 1, niterate= 4, func= "chebyshev", sample = "*")

        iraf.hedit(merge, "INSTRUME", 'TLS-echelle')
        
        for i in range(len(ranges)):
            apertures = int(ranges[i][0])       
            
            new_name = merge.replace(".fits","." + str(apertures) + ".fits")

            input1 = os.path.join(current_directory, merge)
            output = os.path.join(current_directory, new_name)

            iraf.scopy(input1, output, w1=ranges[i][1], w2=ranges[i][2], apertur= apertures, format = "multispec")
                    
            aperturlist.write(new_name + "\n")
    
        aperturlist.close()
    
        new_name_merged = reduced_science_files[k].replace(".extracted.fits",".merged.fits")    

        iraf.scombine("@test.norm.apertures.lis", new_name_merged, group = 'all')

        cut_for_ston = new_name_merged.replace("merged","cut")
        iraf.scopy(new_name_merged, cut_for_ston, w1=5603, w2=5612, format = "multispec",apertures = "" )
        stddev = iraf.imstat(cut_for_ston, Stdout=1, fields="stddev", format="no")
        ston = 1/float(stddev[0])
        
        iraf.hedit(new_name, "STON", ston)        
        
        for i in range(len(ranges)):
            apertures = int(ranges[i][0]) 
            os.remove(os.path.join(merge.replace(".fits","." + str(apertures) + ".fits")))
        
        os.remove(os.path.join(current_directory, "test.norm.apertures.lis"))

        remove_file(current_directory, "norm.dummyI.fits")
        remove_file(current_directory,"norm.dummyIa.fits")
        remove_file(current_directory, "norm.dummyII.fits")
        remove_file(current_directory, "norm.dummyIIa.fits")
        remove_file(current_directory, "norm.dummyIII.fits")
        remove_file(current_directory, "norm.dummyIV.fits")
        remove_file(current_directory, "norm.dummyV.fits")
        remove_file(current_directory, "norm.dummy1.fits")
        remove_file(current_directory, "norm.dummy2.fits")
        remove_file(current_directory, "norm.dummy3.fits")
        remove_file(current_directory, "norm.dummy4.fits")
        remove_file(current_directory, "norm.dummy5.fits")
        remove_file(current_directory, "norm.dummy6.fits")
        remove_file(current_directory, "norm.dummy7.fits")
        remove_file(current_directory, "norm.dummy8.fits")       
    print("Success! Produced files *merged.fits, *norm.fits")



if __name__ == '__main__':
    root = tk.Tk()
    root.title("Tautenburg Spectroscopy Pipeline")
    root.config(bg='white')
    root.grid()
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    nb = ttk.Notebook(root)
    
    page1 = Calibration_Frame(nb, "Flat")
    page2 = Calibration_Frame(nb, "ThAr")
    page3 = Calibration_Frame(nb, "Reduce")
    page4 = Calibration_Frame(nb, "Prepare")
    page5 = Calibration_Frame(nb, "Additional")
    
    nb.add(page1, text='Create Flat')
    nb.add(page2, text='Create Thorium Argon')
    nb.add(page3, text='Reduce Data')
    nb.add(page4, text='Prepare Data')
    nb.add(page5, text='Additional tasks')
    
    nb.rowconfigure(0, weight=1)
    nb.columnconfigure(0, weight=1)
    nb.configure(width=300, height=300)
    nb.grid(row=0,column=0,sticky='nwse')  
    root.mainloop()


