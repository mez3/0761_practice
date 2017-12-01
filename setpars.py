import os
import configparser
import re
import shutils
import time

class module:
    file_of_settings = "Settings.ini"
    file_of_det = "det_PF.txt"
    exe = "FDTD_2D.exe"
    savedir = "tests" 
    files = ("Settings.ini",
             "det_PF.txt")
    
    
    def __init__(self, name):
        self.name = name
    
    
    def changepar(self, conf, par, par_value):
        par_r = conf.get(self.name, par)
        string = re.findall(r'\s*;', par_r)[0]
        par_r = re.sub(r'[\w\W]*;', par_value+string, par_r)
        conf.set(self.name, par, par_r)
    
    
    def changesettings(self, dictofargs):
        conf = configparser.RawConfigParser()
        conf.read(self.file_of_settings)
        for key in dictofargs.keys():
            self.changepar(conf, key, str(dictofargs[key]))
        with open(self.file_of_settings, "w") as config:
            conf.write(config)
        config.close()
    
    
    def getvalue(self):
        os.system(self.exe)
        f = open(self.file_of_det)
        result = re.findall(r'\d+', f.read())[0]
        f.close()
        return int(result)
    
    
    def save_settings(self):
        newdir = self.savedir + '\\' + str(time.time())
        os.mkdir(newdir)
        for filename in self.files:
            shutils.copyfile(filename, newdir)
