import os, sys, inspect
import sqlite3 as sql
import pickle as pkl
from datetime import datetime
import re

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import numpy.ma as ma
from numpy.random import uniform, seed

HOME = os.getenv('HOME','/home/meyerson')
BOUT_TOP = os.getenv('BOUT_TOP','/home/meyerson/BOUT')
SCRATCH =  os.getenv('SCRATCH','/tmp')
PWD = os.getenv('PWD','/tmp')
WORK = os.getenv('WORK','/work/01523/meyerson/')
utc = datetime.utcnow()

sys.path.append(HOME+'/local/python')

cmd_folder = WORK+'/BOUT_sims/blob_py'
if cmd_folder not in sys.path:
     sys.path.insert(0, cmd_folder)

print HOME

import numpy as np
#from blob_info import blob_info, Blob2D
from frame import Frame, FrameMovie
from blob_info import blob_info, Blob2D

#we can open 1 file very easily
from paraview.simple import *
from array import *


def vtk_to_array(vtk_array):
    at = vtk_array.GetDataType()
    if at == 11:
        #vtkDoubleArray
        pt='d'
    elif at == 12:
        #vtkIdTypeArray
        pt='l'
    #this is slow. numpy.zeros would be faster.
    r = array(pt, [0]*vtk_array.GetSize())
    vtk_array.ExportToVoidPointer(r)
    return r


#start reading
def get_pos(i):
    print base_dir+all_pvtus[i]
    r.FileName = base_dir+all_pvtus[i]
    r.UpdatePipeline()

    z = servermanager.Fetch(r)

    pos =  z.GetPoints()

    xmin,xmax,ymin,ymax,zmin,zmax = np.round(pos.GetBounds())
        
    print xmin,xmax,ymin,ymax,zmin,zmax
    numPoints = pos.GetNumberOfPoints()

    allx = []
    ally = []
     
   

    for i in range(numPoints):
        x, y, z = pos.GetPoint(i)
        allx.append(x)
        ally.append(y)

    flaty = np.array(sorted(set(ally)))
    flatx = np.array(sorted(set(allx)))
       
    
    x,y = flatx,flaty

    print (x-np.roll(x,1))[1:-2]
    dx = np.mean((x-np.roll(x,1))[1:-2])

    dy = np.mean((y-np.roll(y,1))[1:-2])
    dy = (y.max() - y.min())/np.round((y.max() - y.min())/dy)
    
    y = np.linspace(y.min(),y.max(),1.*np.round((y.max() - y.min())/dy))

    return {'orig':(allx,ally),'new':(x,y)},dx,dy,xmax-xmin,ymax-ymin


def read_data(base_dir,all_files,pos,cached=False):
    
    
    if cached:
        Hist = (np.load('lastX.npy')).item()
        x = Hist['x']
        y = Hist['y']
        out = Hist['n']
        return out

    out = []
        
    for files in all_files:
        print base_dir+files


        r.FileName = base_dir+files
        r.UpdatePipeline()
            
        z = servermanager.Fetch(r)
        
        pdata = z.GetPointData()
 
        a = pdata.GetArray('alpha_0')

        data = vtk_to_array(a)

        x,y = pos

        allx = pos['orig'][0]
        ally = pos['orig'][1]
        
        x = pos['new'][0]
        y = pos['new'][1]
        
        out.append(griddata((allx, ally),data, (x[None,:], y[:,None]), method='cubic'))
        
    out = np.array(out)   

    Hist={'x':x,'y':y,'n':out}
    np.save('lastX',Hist)

    return out

def gamma_theory(ny,dky,mu = 1.0e-2 ,alpha = 3.0e-5,beta = 6.0e-4,
                 Ln = 130.0/4.0, n0 = 10.0):
     allk = dky*np.arange(ny)+(1e-8*dky)

     ii = complex(0,1)
     soln = {}
     soln['freq'] = []
     soln['gamma'] = []
     soln['gammamax'] = []
     soln['freqmax'] = []

     for i,k in enumerate(allk):
          M = np.zeros([2,2],dtype=complex)
     #density
          M[0,0] = -ii*mu*(k**2)
          M[0,1] = k*n0/Ln
     #potential
          M[1,0] = -beta/(n0*k)
          M[1,1] = -ii*(alpha + mu*k**4)/(k**2)
     #M = M.transpose()
          eigsys= np.linalg.eig(M)  
          gamma = (eigsys)[0].imag
          omega =(eigsys)[0].real
          eigvec = eigsys[1]
     #print 'k: ',k
          
          soln['gamma'].append(gamma)
          soln['gammamax'].append(max(gamma))
          where = ((gamma == gamma.max()))
          soln['freqmax'].append(omega[where])
          soln['freq'].append(omega)
     
     #return the analytic solution     
     return soln

#this reader can deal with .vptu files
r = XMLPartitionedUnstructuredGridReader()

#generate some filenames
base_dir = "/work/01523/meyerson/linear_Arcon/density_output/"
field = "density"
all_pvtus = []

for files in os.listdir(base_dir):
    if files.endswith(".pvtu"):
        all_pvtus.append(files)
        

all_pvtus =  sorted(all_pvtus, key = lambda x: int(re.split('-|\.',x)[1]))         
         
pos,dx,dy,Lx,Ly = get_pos(0)

#print all_pvtus[0:1]
data = read_data(base_dir,all_pvtus,pos,cached=True)

nt,nx,ny = data.shape


fftn = np.fft.fft2(data)
Ak = np.sqrt(fftn.conj()*fftn)


pp = PdfPages('summary.pdf')
fig = plt.figure()

frm_data = Frame(data[-1,:,:],meta={'mask':True,'dx':dx,'dy':dy,'title':'n',
                              'stationary':True})
frm_data.render(fig,221)

dkx = 1.
dky = (2.*np.pi)/Ly

power = Frame(np.real(Ak)[-1,0:30,0:30],meta={'mask':True,'dx':dkx,'dy':dky,'title':'n',
                              'stationary':True})
power.render(fig,222)

dt = 20
time = dt*np.arange(nt)
gamma_num = (np.gradient(np.log(np.real(Ak)))[0])/(np.gradient(time)[0])
gamma_ave = np.mean(gamma_num[-80:-20,:,:],axis=0)

analytic_soln =  gamma_theory(ny,dky)
gamma_th = Frame(np.array(analytic_soln['gammamax'][1:ny/3]),
                 meta={'dx':dky,'x0':dky,'stationary':True,'yscale':'linear',
                       'title':r'$\gamma$','fontsz':20,
                       'ylabel':r'$\frac{\omega}{\omega_{ci}}$',
                       'xlabel':r'$k_y$','ticksize':14})


gamma_num = Frame(gamma_num[:,1:ny/3,0],meta={'dx':dky,'xlabel':r'$k_y$',
                          'title':r'$\gamma$',
                          'ylabel':r'$\frac{\omega}{\omega_{ci}}$',
                          'x0':dky,'shareax':False,'style':'ro', 
                                       'stationary':False,'ticksize':14,'fontsz':20})
gamma_num.t = nt -.7*nt
gamma_num.render(fig,223) 
gamma_th.render(fig,223)

print data.shape
namp = abs(data).max(1).max(1) 
namp =  Frame(namp,meta={'ticksize':14})
namp.render(fig,224)

fig.savefig(pp,format='pdf')
plt.close(fig)
pp.close()


pp = PdfPages('gamma.pdf')
fig = plt.figure()
#gamma_num.render(fig,111) 
gamma_th.ax = None
#gamma.t = 20
gamma_num.ax = None
gamma_th.render(fig,111)
gamma_num.render(fig,111)
fig.savefig(pp,format='pdf')
plt.close(fig)
pp.close()


#let's make a movie
fig = plt.figure()
gamma_num.t = 0
gamma_num.ax = fig.add_subplot(111)
gamma_th.ax = gamma_num.ax
#gamma_th.ax = fig.add_subplot(111)
FrameMovie([gamma_num,gamma_th],fast=True,moviename='gamma',fps=6,fig=fig)
