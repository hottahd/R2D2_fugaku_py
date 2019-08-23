import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import R2D2
import config as c
import sys

try:
    caseid
except NameError:
    print("input caseid id (3 digit)")
    caseid = 0
    caseid = input()
    caseid = "d"+caseid.zfill(3)

dir="../run/"+caseid+"/data/"

R2D2.read_init(dir,"3d")
for key in c.p:
    exec('%s = %s%s%s' % (key, 'c.p["',key,'"]'))

try:
    n0
except NameError:
    n0 = 0
if  n0 > c.p["nd"]:
    n0 = c.p["nd"]

print("Maximum time step= ",nd," time ="\
      ,dtout*float(nd)/3600./24.," [day]")

plt.close('all')


rormst = np.zeros((ix,nd-n0+1))
sermst = np.zeros((ix,nd-n0+1))
prrmst = np.zeros((ix,nd-n0+1))
termst = np.zeros((ix,nd-n0+1))

romt = np.zeros((ix,nd-n0+1))
semt = np.zeros((ix,nd-n0+1))
prmt = np.zeros((ix,nd-n0+1))
temt = np.zeros((ix,nd-n0+1))

fet = np.zeros((ix,nd-n0+1))
fkt = np.zeros((ix,nd-n0+1))
frt = np.zeros((ix,nd-n0+1))
ftt = np.zeros((ix,nd-n0+1))

for n in range(n0,nd+1):
#for n in range(0,1):
    print(n)
    ##############################
    # read time
    f = open(dir+"time/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,endian+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")
    
    ##############################
    # read value
    vc = R2D2.read_vc(dir,n)
        
    ##############################    
    fsun = 6.318e10
    fe = np.average(vc["fe"],axis=1)
    fk = np.average(vc["fk"],axis=1)
    fr = np.average(vc["fr"],axis=1)

    ft = fe + fk + fr

    romt[:,n-n0] = np.sqrt(np.average(vc["rom"]**2,axis=1))
    semt[:,n-n0] = np.sqrt(np.average(vc["sem"]**2,axis=1))
    prmt[:,n-n0] = np.sqrt(np.average(vc["prm"]**2,axis=1))
    temt[:,n-n0] = np.sqrt(np.average(vc["tem"]**2,axis=1))

    fet[:,n-n0] = fe
    fkt[:,n-n0] = fk
    frt[:,n-n0] = fr
    ftt[:,n-n0] = ft

    fontsize = 12
    fmax = 3.0
    fmin = -2.0

    plt.rcParams["font.size"] = 15
    fig1 = plt.figure(num=1,figsize=(5,5))
    ax1 = fig1.add_subplot(111)
    ax1.plot(c.p["xr"],fe/fsun,color="red",label="$F_\mathrm{e}$")
    ax1.plot(c.p["xr"],fk/fsun,color="green",label="$F_\mathrm{k}$")
    ax1.plot(c.p["xr"],fr/fsun,color="blue",label="$F_\mathrm{r}$")
    ax1.plot(c.p["xr"],ft/fsun,color="black",label="$F_\mathrm{t}$")
    ax1.set_xlim(c.p["xmin"]/c.p["rsun"],c.p["xmax"]/c.p["rsun"])
    ax1.set_ylim(fmin,fmax)
    ax1.set_xlabel("$x/R_{\odot}$")
    ax1.set_ylabel("$F/F_{\odot}$")
    ax1.set_title("Full convection zone")
    ax1.legend(loc='upper left',prop={'size': 10})
    ax1.annotate(s="t="+"{:.2f}".format(t[0]/3600./24.)+" [day]"\
                     ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)
    
    if n == n0:
        plt.tight_layout()
    plt.pause(0.001)
    #plt.draw()
    
    if n != nd:
        plt.clf() # clear figure

    # loop end
    ###############################################################################
    ###############################################################################
    ###############################################################################

rom = np.average(romt,axis=1)
sem = np.average(semt,axis=1)
prm = np.average(prmt,axis=1)
tem = np.average(temt,axis=1)

fe = np.average(fet,axis=1)
fk = np.average(fkt,axis=1)
fr = np.average(frt,axis=1)
ft = fe + fk + fr
         
plt.rcParams["font.size"] = 15
fig2 = plt.figure(num=2,figsize=(5,5))
ax2 = fig2.add_subplot(111)
ax2.plot(xr,fe/fsun,color="red",label="$F_\mathrm{e}$")
ax2.plot(xr,fk/fsun,color="green",label="$F_\mathrm{k}$")
ax2.plot(xr,fr/fsun,color="blue",label="$F_\mathrm{r}$")
ax2.plot(xr,ft/fsun,color="black",label="$F_\mathrm{t}$")
ax2.set_xlim(xmin/rsun,xmax/rsun)
ax2.set_ylim(fmin,fmax)
ax2.set_xlabel("$x/R_{\odot}$")
ax2.set_ylabel("$F/F_{\odot}$")
ax2.set_title("Full convection zone")
ax2.legend(loc='upper left',prop={'size': 15})
ax2.annotate(s="t="+"{:.2f}".format(t[0]/3600./24.)+" [day]"\
                 ,xy=[0.01,0.01],xycoords="figure fraction",fontsize=18)

fig2.tight_layout()
plt.pause(0.001)
plt.ion()
    
