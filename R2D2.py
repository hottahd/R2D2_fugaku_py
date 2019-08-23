# R2D2.py
######################################################
######################################################
######################################################
def read_init(dir,dimension):
    import numpy as np
    import config as c
    import sys
    c.p = {"a":0}
    f = open(dir+"param/nd.dac","r")
    nn = f.read().split()
    nd = int(nn[0])
    f.close()
    
    c.p['nd'] = nd
    
    R2D2_py_ver = 0.1
    f = open(dir+"param/params.dac","r")
    line = f.readline().split()
    if R2D2_py_ver != float(line[2]):
        print("#######################################################")
        print("#######################################################")
        print("### Current R2D2 Python version is ",R2D2_py_ver,".")
        print("### You use the data from R2D2 version ",float(line[2]),".")
        print("### Please use the same version of fortran R2D2.")
        print("#######################################################")
        print("#######################################################")
        sys.exit()
    
    line = f.readline()
    while line:
        if line.split()[2] == 'i':
            c.p[line.split()[1]] = int(line.split()[0])
        if line.split()[2] == 'd':
            c.p[line.split()[1]] = float(line.split()[0])
        line = f.readline()

    f.close()
            
    if c.p["swap"] == 0: # little
        c.p["endian"] = "<"
    else:
        c.p["endian"] = ">"

    c.p["ix"] = c.p["ix0"]*c.p["nx"]
    c.p["jx"] = c.p["jx0"]*c.p["ny"]
    c.p["kx"] = c.p["kx0"]*c.p["nz"]

    marginx = c.p["margin"]*(c.p["xdcheck"]-1)
    marginy = c.p["margin"]*(c.p["ydcheck"]-1)
    marginz = c.p["margin"]*(c.p["zdcheck"]-1)
   
    ixg = c.p["ix"] + 2*marginx
    jxg = c.p["jx"] + 2*marginy
    kxg = c.p["kx"] + 2*marginz
        
    endian = c.p["endian"]
    dtyp=np.dtype([ \
                    ("head",endian+"i"),\
                    ("x",endian+str(ixg)+"d"),\
                    ("y",endian+str(jxg)+"d"),\
                    ("z",endian+str(kxg)+"d"),\
                    ("ro0",endian+str(ixg)+"d"),\
                    ("te0",endian+str(ixg)+"d"),\
                    ("pr0",endian+str(ixg)+"d"),\
                    ("en0",endian+str(ixg)+"d"),\
                    ("dprdro",endian+str(ixg)+"d"),\
                    ("dprdse",endian+str(ixg)+"d"),\
                    ("dtedro",endian+str(ixg)+"d"),\
                    ("dtedse",endian+str(ixg)+"d"),\
                    ("dendro",endian+str(ixg)+"d"),\
                    ("dendse",endian+str(ixg)+"d"),\
                    ("gx",endian+str(ixg)+"d"),\
                    ("xi",endian+str(ixg)+"d"),\
                    ("fa",endian+str(ixg)+"d"),\
                    ("tail",endian+"i")\
    ])
    f = open(dir+"back.dac",'rb')
    back = np.fromfile(f,dtype=dtyp,count=1)
    f.close()

    for key in back.dtype.names:
        if back[key].size == ixg:
            c.p[key] = back[key].reshape((ixg),order="F")[marginx:ixg-marginx]
        elif back[key].size == jxg:
            c.p[key] = back[key].reshape((jxg),order="F")[marginy:jxg-marginy]
        elif back[key].size == kxg:
            c.p[key] = back[key].reshape((kxg),order="F")[marginz:kxg-marginz]

    c.p["rsun"] = 6.9598947e+10
    c.p["xr"] = c.p["x"]/c.p["rsun"]

    ix0 = c.p['ix0']
    jx0 = c.p['jx0']
    kx0 = c.p['kx0']
    dtyp=np.dtype([ \
                    ("head",endian+"i"),\
                    ("xyz",endian+str(ix0*jx0*kx0*3)+"i"),\
                    ("tail",endian+"i")\
    ])
    f = open(dir+"param/xyz.dac",'rb')
    xyz = np.fromfile(f,dtype=dtyp,count=1)['xyz']
    f.close()
    c.p['xyz'] = xyz.reshape((ix0*jx0*kx0,3),order='F')

    ##############################
    # read value information
    if dimension == "3d":
        f = open(dir+"remap/c.dac","r")
        value = f.read().split('\n')
        c.p["m2da"] = int(value[0])
        del value[0]
        c.p["cl"] = list(map(str.strip,value)) ## strip space from character
        f.close()

        ##############################
        # read mpi information
        npe = c.p["npe"]
        dtyp=np.dtype([ \
                ("iss",endian+str(npe)+"i4"),\
                    ("iee",endian+str(npe)+"i4"),\
                    ("jss",endian+str(npe)+"i4"),\
                    ("jee",endian+str(npe)+"i4"),\
                    ("iixl",endian+str(npe)+"i4"),\
                    ("jjxl",endian+str(npe)+"i4"),\
                    ("np_ijr",endian+str(c.p["ixr"]*c.p["jxr"])+"i4"),\
                    ("ir",endian+str(npe)+"i4"),\
                    ("jr",endian+str(npe)+"i4"),\
                    ("i2ir",endian+str(ixg)+"i4"),\
                    ("j2jr",endian+str(jxg)+"i4"),\
                    ])
        
        f = open(dir+"remap/remap_info.dac",'rb')
        mpi = np.fromfile(f,dtype=dtyp,count=1)    
        f.close()

        for key in mpi.dtype.names:
            if key == "np_ijr":
                c.p[key] = mpi[key].reshape((c.p["ixr"],c.p["jxr"]),order="F")
            else:
                c.p[key] = mpi[key].reshape((mpi[key].size),order="F")
    
        c.p["i2ir"] = c.p["i2ir"][marginx:ixg-marginx]
        c.p["j2jr"] = c.p["j2jr"][marginy:jxg-marginy]

        c.p["iss"] = c.p["iss"] - 1
        c.p["iee"] = c.p["iee"] - 1
        c.p["jss"] = c.p["jss"] - 1
        c.p["jee"] = c.p["jee"] - 1
    
######################################################
def read_qq_original(dir,n):
    import numpy as np
    import config as c

    npe = c.p['npe']
    xyz = c.p['xyz']
    endian = c.p['endian']
    mtype = c.p['mtype']
    ix = c.p['ix']
    jx = c.p['jx']
    kx = c.p['kx']
    nx = c.p['nx']
    ny = c.p['ny']
    nz = c.p['nz']

    margin = c.p['margin']
    
    nxg = c.p['nx'] + 2*margin
    nyg = c.p['ny'] + 2*margin
    nzg = c.p['nz'] + 2*margin

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'qq' in c.q:
        memflag = not c.q['qq'].shape == (mtype,ix,jx,kx)
    if 'qq' not in c.q or memflag:
        print('memory is newly allocated')
        c.q["qq"] = np.zeros((mtype,ix,jx,kx))

    dtyp = np.dtype([('qq',endian+str(mtype*nxg*nyg*nzg)+"d")])
    
    for np0 in range(0,npe):
        f = open(dir+'qq/qq.dac.'+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        ib = xyz[np0,0]
        jb = xyz[np0,1]
        kb = xyz[np0,2]
        c.q['qq'][:,ib*nx:(ib+1)*nx,jb*ny:(jb+1)*ny,kb*nz:(kb+1)*nz] = \
            np.fromfile(f,dtype=dtyp,count=1)['qq'].reshape((mtype,nxg,nyg,nzg),order="F") \
            [:,margin:nxg-margin,margin:nyg-margin,margin:nzg-margin]
        f.close()

    return c.q
        
def read_time(dir,n):
    import numpy as np
    import config as c
    f = open(dir+"time/t.dac."+'{0:08d}'.format(n),"rb")
    t = np.fromfile(f,c.p['endian']+'d',1)
    f.close()    
    t = np.reshape(t,(1),order="F")[0]

    return t    

######################################################
######################################################
######################################################
### read horizontal variable
### prepare array qq = np.zeros((mtype+3,jx,kx))
def read_qq_select(dir,xs,n):
    import numpy as np
    import config as c
    i0 = np.argmin(np.abs(c.p["x"]-xs))
    ir0 = c.p["i2ir"][i0]
    mtype = c.p["mtype"]
    iixl = c.p["iixl"]
    jjxl = c.p["jjxl"]
    jx = c.p["kx"]
    kx = c.p["kx"]
    iss = c.p["iss"]
    jss = c.p["jss"]
    jee = c.p["jee"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in c.q2:
        memflag = not c.q2['ro'].shape == (jx,kx)
    if 'ro' not in c.q2 or memflag:
        print('memory is newly allocated')
        c.q2["ro"] = np.zeros((jx,kx))
        c.q2["vx"] = np.zeros((jx,kx))
        c.q2["vy"] = np.zeros((jx,kx))
        c.q2["vz"] = np.zeros((jx,kx))
        c.q2["se"] = np.zeros((jx,kx))
        c.q2["pr"] = np.zeros((jx,kx))
        c.q2["te"] = np.zeros((jx,kx))
    for jr0 in range(1,c.p["jxr"]+1):
        np0 = c.p["np_ijr"][ir0-1,jr0-1]
        dtyp=np.dtype([ \
                ("qq",c.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("pr",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                ("te",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                    ])
        f = open(dir+"remap/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
        qqq = np.fromfile(f,dtype=dtyp,count=1)
        c.q2["ro"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,i0-iss[np0],:,:]
        c.q2["vx"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,i0-iss[np0],:,:]
        c.q2["vy"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,i0-iss[np0],:,:]
        c.q2["vz"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,i0-iss[np0],:,:]
        c.q2["se"][jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,i0-iss[np0],:,:]
        c.q2["pr"][jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        c.q2["te"][jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")[i0-iss[np0],:,:]
        f.close()
    
    return c.q2

######################################################
######################################################
######################################################
### read horizontal variable
### prepare array qq = np.zeros((mtype+3,ix,jx,kx))
def read_qq(dir,n):
    import numpy as np
    import config as c
    
    mtype = c.p["mtype"]
    iixl = c.p["iixl"]
    jjxl = c.p["jjxl"]
    kx = c.p["kx"]
    iss = c.p["iss"]
    iee = c.p["iee"]
    jss = c.p["jss"]
    jee = c.p["jee"]
    ix = c.p["ix"]
    jx = c.p["jx"]
    kx = c.p["kx"]

    ### Only when memory is not allocated 
    ### and the size of array is different
    ### memory is allocated
    memflag = True
    if 'ro' in c.q3:
        memflag = not c.q3['ro'].shape == (ix,jx,kx)
    if 'ro' not in c.q3 or memflag:
        print('memory is newly allocated')
        c.q3["ro"] = np.zeros((ix,jx,kx))
        c.q3["vx"] = np.zeros((ix,jx,kx))
        c.q3["vy"] = np.zeros((ix,jx,kx))
        c.q3["vz"] = np.zeros((ix,jx,kx))
        c.q3["se"] = np.zeros((ix,jx,kx))
        c.q3["pr"] = np.zeros((ix,jx,kx))
        c.q3["te"] = np.zeros((ix,jx,kx))

    for ir0 in range(1,c.p["ixr"]+1):
        for jr0 in range(1,c.p["jxr"]+1):
            np0 = c.p["np_ijr"][ir0-1,jr0-1]
            dtyp=np.dtype([ \
                    ("qq",c.p["endian"]+str(mtype*iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("pr",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ("te",c.p["endian"]+str(iixl[np0]*jjxl[np0]*kx)+"f"),\
                        ])
            f = open(dir+"remap/qq.dac."+'{0:08d}'.format(n)+"."+'{0:08d}'.format(np0),'rb')
            qqq = np.fromfile(f,dtype=dtyp,count=1)
            c.q3["ro"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[0,:,:,:]
            c.q3["vx"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[1,:,:,:]
            c.q3["vy"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[2,:,:,:]
            c.q3["vz"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[3,:,:,:]
            c.q3["se"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["qq"].reshape((mtype,iixl[np0],jjxl[np0],kx),order="F")[4,:,:,:]
            c.q3["pr"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["pr"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            c.q3["te"][iss[np0]:iee[np0]+1,jss[np0]:jee[np0]+1,:] = qqq["te"].reshape((iixl[np0],jjxl[np0],kx),order="F")
            f.close()
    
    return c.q3

##############################
# read remap_calc variable
def read_vc(dir,n):
    import numpy as np
    import config as c

    f = open(dir+"remap/vla.dac."+'{0:08d}'.format(n),"rb")
    vl0 = np.fromfile(f,c.p["endian"]+'f',c.p['m2da']*c.p['ix']*c.p['jx'])
    f.close()

    vl = np.reshape(vl0,(c.p['ix'],c.p['jx'],c.p['m2da']),order="F")

    vc = {'a':0}
    for m in range(c.p["m2da"]):
        vc[c.p["cl"][m]] = vl[:,:,m]
    
    return vc
