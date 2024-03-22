import logging
import time
import numpy as np
import pandas as pd
import scipy.signal as sg
import scipy.interpolate as sci
import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.backends.backend_pdf import PdfPages
from glob import glob
import h5py
import pickle as pkl
# import matlab.engine
import subprocess as sub
import dask
from dask_jobqueue import SLURMCluster
from dask.distributed import Client

def IsNumeric(x):
    try:
        float(x)
        return True
    except ValueError:
        return False

def AllocWorkers(nclust, ncore):
    cluster = SLURMCluster(cores=ncore, memory="72GiB", walltime="4-00:00:00", job_extra=["--qos=conklin-b"])
    cluster.scale(nclust*ncore)
    client = Client(cluster)
    while ((client.status == "running") and (len(client.scheduler_info()["workers"]) < 1)):
        time.sleep(1.0)
    client.forward_logging()
    return client

# def RunSimEng(args):
#     # logging.info(f"Starting sim {row.name}.")
#     common_args = ['config','nominal.xml','savestats',0]
#     # args = sum(row.to_dict().items(),())
#     eng = matlab.engine.start_matlab()
#     res = pd.Series(eng.TMrel_HG(*common_args,*args,background=False), index=["gapx","gapy","gapz"])
#     return res

def RunSim(args, limit=4*3600):
    logging.basicConfig(
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    common_args = ['config','nominal.xml','savestats','0','savedata','1','run_name','"'+("_".join(map(str, args)))+'"']
    all_args = [arg if IsNumeric(arg) else '\''+arg+'\'' for arg in (common_args + list(args))]
    start = time.time()
    try:
        res = sub.check_output(["./RunSub.sh"] + all_args, timeout=limit)
        stop = time.time()
        gaps = res.splitlines()[-2].decode().split()
    except sub.TimeoutExpired:
        logging.info(f"Run timed out. Args: {args}")
        gaps = [np.nan,np.nan,np.nan,limit]
    if len(gaps) != 4:
        logging.warn(res)
        gaps = [np.nan,np.nan,np.nan,stop-start]
    return pd.Series(gaps, index=["gapx","gapy","gapz","runtime"], dtype=float)

def GetContact(fh, sTM, Lt):
    psi = fh["bTM"][0,:] # angle about x
    the = fh["bTM"][1,:] # angle about y
    phi = fh["bTM"][2,:] # angle about z

    # rotational matrix of TM
    R = np.array([[np.cos(the)*np.cos(phi), np.cos(the)*np.sin(phi), -np.sin(the)],
                  [np.sin(psi)*np.sin(the)*np.cos(phi)-np.cos(psi)*np.sin(phi), np.sin(psi)*np.sin(the)*np.sin(phi)+np.cos(psi)*np.cos(phi), np.sin(psi)*np.cos(the)],
                  [np.cos(psi)*np.sin(the)*np.cos(phi)+np.sin(psi)*np.sin(phi), np.cos(psi)*np.sin(the)*np.sin(phi)-np.sin(psi)*np.cos(phi), np.cos(psi)*np.cos(the)]])

    s1 = np.dot(np.transpose(R,(0,2,1)),np.array([0,sTM/2,0])).T # vector to +Y face
    s2 = np.dot(np.transpose(R,(0,2,1)),np.array([0,-sTM/2,0])).T # vector to -Y face
    s1h = s1/np.sqrt(np.sum(s1*s1,axis=1))[:,None] # unit vector of s1
    s2h = s2/np.sqrt(np.sum(s2*s2,axis=1))[:,None] # unit vector of s2

    # Position and Speed Vectors

    rt1 = np.zeros((fh["xT1"][:].flatten().shape[0],3)) # RT1 position
    rt1[:,1] = fh["xT1"][:]+sTM/2
    rt2 = np.zeros((fh["xT2"][:].flatten().shape[0],3)) # RT2 position
    rt2[:,1] = fh["xT2"][:]-sTM/2
    # vt1 = np.pad(fh["vT1"][:],((1,1),(0,0))) # RT1 speed
    # vt2 = np.pad(fh["vT2"][:],((1,1),(0,0))) # RT2 speed

    rg1 = np.zeros((fh["xG1"][:].flatten().shape[0],3)) # GF1 position
    rg1[:,1] = fh["xG1"][:]+sTM/2
    rg2 = np.zeros((fh["xG2"][:].flatten().shape[0],3)) # GF2 position
    rg2[:,1] = fh["xG2"][:]-sTM/2
    # vg1 = np.array([[0, vG1, 0]]) # GF1 speed
    # vg2 = np.array([[0, vG2, 0]]) # GF2 speed

    rTM = fh["xTM"][:].T # TM position
    omTM = fh["wTM"][:].T # TM angular speed [rad/s]
    vTM = fh["vTM"][:].T # TM linear speed [m/s]

    # Contact Forces

    # GF to RT contact
    Fc1 = (fh["xG1"][:]-fh["xT1"][:]>Lt)
    Fc2 = (fh["xG2"][:]-fh["xT2"][:]>Lt)

    # GF to TM contact
    lg1 = rg1-(rTM+s1) # GF contact point from TM Y-face center
    lg2 = rg2-(rTM+s2)
    dLg1 = np.sum(lg1*s1h,axis=1) # GF to TM face distance along TM face normal
    dLg2 = np.sum(lg2*s2h,axis=1)

    # RT to TM contact
    lt1 = rt1-rTM-s1 # RT contact point from TM Y-face center
    lt2 = rt2-rTM-s2
    dLt1 = np.sum(lt1*s1h,axis=1) # RT to TM face distance along TM face normal
    dLt2 = np.sum(lt2*s2h,axis=1)

    # Check case conditions % MOD
    cG_TM1 = dLg1<0 # RT1 to TM contact
    cG_TM2 = dLg2<0 # RT2 to TM contact
    cRT_TM1 = dLt1<0 # RT1 to TM contact
    cRT_TM2 = dLt2<0 # RT2 to TM contact
    cG_RT1 = Fc1>0 # GF1 to RT1 contact
    cG_RT2 = Fc2>0 # GF2 to RT2 contact
    contact = np.vstack([cG_RT1,cG_RT2,cG_TM1,cG_TM2,cRT_TM1,cRT_TM2,fh["bond1"][:].flatten(),fh["bond2"][:].flatten()]).T
    return contact

def PlotPositions(files, labels):
    dsamp = 50
    dimnames = ["X","Y","Z"]
    fhs = [h5py.File(file, "r") for file in files]
    with PdfPages("SpringScale.pdf") as export_pdf:
        fig, axs = plt.subplots(figsize=(8,8),nrows=3,sharex=True)
        for dim in range(3):
            axs[dim].set_prop_cycle((cycler(color=['r', 'b', 'g']) + cycler(linestyle=['-', '--', ':'])))
            for i in range(len(files)):
                logging.info(f"Plotting {files[i]}")
                t = fhs[i]["t"][0,:]
                axs[dim].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xTM"][dim,::dsamp],"-",label=labels[i])
            axs[dim].set_ylabel(f"TM {dimnames[dim]} pos. [m]")
        plt.legend()
        axs[2].set_xlabel("time [s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

        fig, ax = plt.subplots(figsize=(8,6))
        ax.set_prop_cycle((cycler(color=['r', 'b', 'g']) + cycler(linestyle=['-', '--', ':'])))
        for i in range(len(files)):
            logging.info(f"Plotting {files[i]}")
            t = fhs[i]["t"][0,:]
            tclip = (t>18.5)&(t<18.6)
            ax.plot(fhs[i]["t"][0,tclip],fhs[i]["vTM"][1,tclip],"-",label=labels[i])
        ax.set_ylabel(f"TM y pos. [m]")
        plt.legend()
        ax.set_xlabel("time [s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

def GetRuntime(output):
    if len(output) > 0:
        output = output[0]
    else:
        return np.nan
    with open(output,"r") as fh:
        for line in fh:
            if line.startswith("Simulation runtime"):
                return float(line.split(" = ")[1].split(" ")[0])
    return np.nan

def UnpackVals(fh):
    fh["t"] = fh["t_out"][:]
    fh["xG1"] = fh["y_out"][0,:]
    fh["xG2"] = fh["y_out"][1,:]
    fh["xT1"] = fh["y_out"][2,:]
    fh["vT1"] = fh["y_out"][3,:]
    fh["xT2"] = fh["y_out"][4,:]
    fh["vT2"] = fh["y_out"][5,:]
    fh["xTM"] = fh["y_out"][6:9,:]
    fh["vTM"] = fh["y_out"][9:12,:]
    fh["bTM"] = fh["y_out"][12:15,:]
    fh["wTM"] = fh["y_out"][15:18,:]
    fh["bond1"] = fh["y_out"][18,:]
    fh["bond2"] = fh["y_out"][19,:]

def GetErrors(grp, sTM, Lt, export_pdf=None):
    dsamp = 50
    dimnames = ["X","Y","Z"]
    arg_list = [sum(grp.loc[i].to_dict().items(),()) for i in grp.index]
    grp["files"] = [glob(f"{('_'.join(map(str, args)))}/*.mat")[0] for args in arg_list]
    grp["runtime"] = [GetRuntime(glob(f"{('_'.join(map(str, args)))}/Output*.txt")) for args in arg_list]
    exact = h5py.File(grp["files"][grp["estimateRT"]=="0.0"].iloc[0], "a")
    if "t" not in exact.keys():
        UnpackVals(exact)
    apprx = h5py.File(grp["files"][grp["estimateRT"]=="1.0"].iloc[0], "a")
    if "t" not in apprx.keys():
        UnpackVals(apprx)

    conExact = GetContact(exact, sTM, Lt)
    cngExact = np.where(np.any(np.diff(conExact,axis=0)!=0,axis=1))[0]
    conApprx = GetContact(apprx, sTM, Lt)
    cngApprx = np.where(np.any(np.diff(conApprx,axis=0)!=0,axis=1))[0]
    onApprx = (conApprx[:,4]==1) & (conApprx[:,5]==1) & (conApprx[:,2]==0) & (conApprx[:,3]==0)
    logging.info(onApprx)
    begApprx = np.where(onApprx[:-1] < onApprx[1:])[0]
    endApprx = np.where(onApprx[:-1] > onApprx[1:])[0]
    # begApprx = np.where(np.diff(onApprx)==1)[0]
    # endApprx = np.where(np.diff(onApprx)==-1)[0]
    if begApprx.shape[0] == 0:
        begApprx = np.array([0])
    if endApprx.shape[0] == 0:
        endApprx = np.array([len(conApprx)-1])
    if begApprx[0] > endApprx[0]:
        begApprx = np.insert(begApprx,0,0)
    if begApprx[-1] > endApprx[-1]:
        endApprx = np.append(endApprx, len(conApprx)-1)
    regApprx = np.vstack([begApprx,endApprx]).T

    maxerr = []
    for dim in range(3):
        t0 = exact["t"][0,:]
        t = apprx["t"][0,:]
        x0 = exact["xTM"][dim,:]
        x = sci.InterpolatedUnivariateSpline(t,apprx["xTM"][dim,:])
        maxerr += [np.max(np.abs(x(t0[::dsamp])-x0[::dsamp]))]
    maxerr += [grp["runtime"][grp["estimateRT"]=="0.0"].iloc[0]/grp["runtime"][grp["estimateRT"]=="1.0"].iloc[0]]
    if export_pdf is None:
        return np.array(maxerr)

    fig, axs = plt.subplots(figsize=(8,8),nrows=3,sharex=True)
    for dim in range(3):
        t0 = exact["t"][0,:]
        t = apprx["t"][0,:]
        x0 = exact["xTM"][dim,:]
        x = sci.InterpolatedUnivariateSpline(t,apprx["xTM"][dim,:])
        axs[dim].plot(t0[::dsamp],x(t0[::dsamp])-x0[::dsamp],"-")
        # for i in range(regApprx.shape[0]):
        #     # axs[dim].axvspan(*t0[regApprx[i,:]], alpha=0.5, color="red")
        #     axs[dim].axvspan(t0[begApprx[i]], t0[endApprx[i]], alpha=0.5, color="red")
        # ylim = plt.ylim()
        # plt.vlines(cngExact,*ylim,colors="g",linestyles="dotted")
        # plt.vlines(cngApprx,*ylim,colors="r",linestyles="dotted")
        axs[dim].set_ylabel(f"TM {dimnames[dim]} pos. err. [m]")
    axs[2].set_xlabel("time [s]")
    plt.suptitle(f'kG = {float(grp["kG"].iloc[0]):.3g}, kT = {float(grp["kT"].iloc[0]):.3g}')
    plt.tight_layout()
    export_pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(figsize=(8,6), nrows=2)
    axs[0].plot(exact["t"][0,::dsamp], conExact[::dsamp,2], "-", label="GF1 exact")
    axs[1].plot(exact["t"][0,::dsamp], conExact[::dsamp,3], "-", label="GF2 exact")
    axs[0].plot(exact["t"][0,::dsamp], conExact[::dsamp,4], "-", label="RT1 exact")
    axs[1].plot(exact["t"][0,::dsamp], conExact[::dsamp,5], "-", label="RT2 exact")
    axs[0].plot(apprx["t"][0,::dsamp], conApprx[::dsamp,2], "--", label="GF1 apprx")
    axs[1].plot(apprx["t"][0,::dsamp], conApprx[::dsamp,3], "--", label="GF2 apprx")
    axs[0].plot(apprx["t"][0,::dsamp], conApprx[::dsamp,4], "--", label="RT1 apprx")
    axs[1].plot(apprx["t"][0,::dsamp], conApprx[::dsamp,5], "--", label="RT2 apprx")
    axs[1].set_xlabel("time [s]")
    axs[0].set_ylabel("RT/GF 1 contact")
    axs[1].set_ylabel("RT/GF 2 contact")
    plt.legend()
    plt.tight_layout()
    export_pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(figsize=(8,8),nrows=3,sharex=True)
    for dim in range(3):
        t0 = exact["t"][0,exact["t"][0,:]>17.5]
        t = apprx["t"][0,apprx["t"][0,:]>17.5]
        x0 = exact["xTM"][dim,exact["t"][0,:]>17.5]
        x = sci.InterpolatedUnivariateSpline(t,apprx["xTM"][dim,apprx["t"][0,:]>17.5])
        axs[dim].plot(t0[::dsamp],x(t0[::dsamp])-x0[::dsamp],"-")
        # for i in range(regApprx.shape[0]):
        #     # axs[dim].axvspan(*t0[regApprx[i,:]], alpha=0.5, color="red")
        #     axs[dim].axvspan(t0[begApprx[i]], t0[endApprx[i]], alpha=0.5, color="red")
        # ylim = plt.ylim()
        # plt.vlines(cngExact,*ylim,colors="g",linestyles="dotted")
        # plt.vlines(cngApprx,*ylim,colors="r",linestyles="dotted")
        axs[dim].set_ylabel(f"TM {dimnames[dim]} pos. err. [m]")
    axs[2].set_xlabel("time [s]")
    plt.suptitle(f'kG = {float(grp["kG"].iloc[0]):.3g}, kT = {float(grp["kT"].iloc[0]):.3g}')
    plt.tight_layout()
    export_pdf.savefig()
    plt.close()

    fig, axs = plt.subplots(figsize=(8,6), nrows=2)
    axs[0].plot(exact["t"][:,exact["t"][0,:]>17.5][0,::dsamp], conExact[exact["t"][0,:]>17.5,:][::dsamp,2], "-", label="GF1 exact")
    axs[1].plot(exact["t"][:,exact["t"][0,:]>17.5][0,::dsamp], conExact[exact["t"][0,:]>17.5,:][::dsamp,3], "-", label="GF2 exact")
    axs[0].plot(exact["t"][:,exact["t"][0,:]>17.5][0,::dsamp], conExact[exact["t"][0,:]>17.5,:][::dsamp,4], "-", label="RT1 exact")
    axs[1].plot(exact["t"][:,exact["t"][0,:]>17.5][0,::dsamp], conExact[exact["t"][0,:]>17.5,:][::dsamp,5], "-", label="RT2 exact")
    axs[0].plot(apprx["t"][:,apprx["t"][0,:]>17.5][0,::dsamp], conApprx[apprx["t"][0,:]>17.5,:][::dsamp,2], "--", label="GF1 apprx")
    axs[1].plot(apprx["t"][:,apprx["t"][0,:]>17.5][0,::dsamp], conApprx[apprx["t"][0,:]>17.5,:][::dsamp,3], "--", label="GF2 apprx")
    axs[0].plot(apprx["t"][:,apprx["t"][0,:]>17.5][0,::dsamp], conApprx[apprx["t"][0,:]>17.5,:][::dsamp,4], "--", label="RT1 apprx")
    axs[1].plot(apprx["t"][:,apprx["t"][0,:]>17.5][0,::dsamp], conApprx[apprx["t"][0,:]>17.5,:][::dsamp,5], "--", label="RT2 apprx")
    axs[1].set_xlabel("time [s]")
    axs[0].set_ylabel("RT/GF 1 contact")
    axs[1].set_ylabel("RT/GF 2 contact")
    plt.legend()
    plt.tight_layout()
    export_pdf.savefig()
    plt.close()

    return np.array(maxerr)

if __name__ == "__main__":

    logging.basicConfig(
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    dorun = True
    neng = 10 # Number of Matlab instances
    nval = 10 # Number of vals for each var
    nsamp = 10 # Number of runs
    nclust = 8
    ncore = 4
    if dorun:
        client = AllocWorkers(nclust, ncore)
    fs = 3e6
    dt = 1.0/fs

    sTM = 30e-3
    sEH = 32e-3
    gap = (sEH-sTM)/2
    Lt = 0.5e-4
    tr = 4e-6

    vars = ["kG", "kT", "estimateRT"]
    rngs = [[5,7.5], [5,7.5]]
    # arrs = [np.linspace(*rng, nval) for rng in rngs]
    arrs = [np.logspace(*rng, nval) for rng in rngs]
    cols = np.array([col.flatten() for col in np.meshgrid(*arrs)]).T
    rng = np.random.default_rng()
    cols = rng.choice(cols, nsamp, replace=False)
    cols = np.vstack([np.hstack([cols,np.zeros((cols.shape[0],1))]),
                      np.hstack([cols,np.ones((cols.shape[0],1))])])
    df = pd.DataFrame(cols, columns=vars).astype(str)

    # vecs = ["xTM0"]
    # rngs = [[0,1e-4], [0,0], [0,1e-4]]
    # arrs = [np.linspace(*rng, nval) for rng in rngs]
    # vec_cols = np.array([col.flatten() for col in np.meshgrid(*arrs)]).T
    # vec_df = pd.DataFrame(np.array([pd.DataFrame(vec_cols[:,i:(i+3)], columns=["x","y","z"], dtype="str").agg(
    #     lambda x: "[" + ", ".join(x) + "]", axis=1) for i in range(0,vec_cols.shape[1],3)]).T, columns=vecs)
    # df = df.join(vec_df, how="outer").fillna(0)

    if dorun:
        arg_list = [sum(df.loc[i].to_dict().items(),()) for i in df.index]
        res = client.map(lambda args: RunSim(args, 2*24*3600), arg_list)
        logging.info("Awaiting results.")
        res_list = client.gather(res)
        for res in res_list:
            logging.info(res)
        res_df = pd.concat([df,pd.concat(res_list,axis=1).T],axis=1)
        res_df.to_hdf("Gaps.h5","gaps")

    files = glob("kG_*/*.mat")
    df = pd.concat([pd.DataFrame(np.array(file.split("/")[0].split("_")).reshape((-1,2)).T[1:,:],columns=vars)
                    for file in files], axis=0, ignore_index=True)
    # files = [glob(f"{('_'.join(map(str, args)))}/*.mat")[0] for args in arg_list]
    # labels = [" ".join(sum(df[["kG","kT","estimateRT"]].loc[i].to_dict().items(),())) for i in df.index]#["exact", "est.", "old"]
    # logging.info(f"Found {len(files)} files.")
    # PlotPositions(files, labels)
    errors = pd.DataFrame(np.nan,index=pd.MultiIndex.from_frame(df[["kG","kT"]]).drop_duplicates(),columns=["X","Y","Z","runtime"])
    with PdfPages("EstimateError.pdf") as export_pdf:
        for name, grp in df.groupby(["kG","kT"]):
            if grp.shape[0] != 2:
                continue
            errors.loc[name,:] = GetErrors(grp, sTM, Lt,export_pdf)
        errors = errors.reset_index().astype(np.float64)
        
        fig, ax = plt.subplots(figsize=(8,6))
        plt.scatter(errors["kG"],errors["kT"],c=np.log10(errors["X"]))
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("kG")
        plt.ylabel("kT")
        plt.colorbar(label="log max X err.")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()
        
        fig, ax = plt.subplots(figsize=(8,6))
        plt.scatter(errors["kG"],errors["kT"],c=np.log10(errors["Y"]))
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("kG")
        plt.ylabel("kT")
        plt.colorbar(label="log max Y err.")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()
        
        fig, ax = plt.subplots(figsize=(8,6))
        plt.scatter(errors["kG"],errors["kT"],c=np.log10(errors["Z"]))
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("kG")
        plt.ylabel("kT")
        plt.colorbar(label="log max Z err.")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()
        
        fig, ax = plt.subplots(figsize=(8,6))
        plt.scatter(errors["kG"],errors["kT"],c=errors["runtime"])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("kG")
        plt.ylabel("kT")
        plt.colorbar(label="runtime ratio")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()
    exit()

        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         t = fh["t"][0,:]
        #         TMy = fh["xTM"][1,:]
        #         GF1 = fh["xG1"][0,:]
        #         GF2 = fh["xG2"][0,:]
        #         RT1 = fh["xT1"][0,:]
        #         RT2 = fh["xT2"][0,:]
        #         fig, ax = plt.subplots(figsize=(8,6))
        #         plt.plot(t, TMy, "-", label="TM")
        #         plt.plot(t, GF1, "-", label="GF1")
        #         plt.plot(t, GF2, "-", label="GF2")
        #         plt.plot(t, RT1, "-", label="RT1")
        #         plt.plot(t, RT2, "-", label="RT2")
        #         contact = GetContact(fh, sTM, Lt)
        #         change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
        #         xlim = plt.xlim()
        #         ylim = plt.ylim()
        #         plt.vlines(fh["t"][0,change], *ylim, linestyles="dotted", colors="black")
        #         plt.xlim(xlim)
        #         plt.ylim(ylim)
        #         plt.xlabel("time [s]")
        #         plt.ylabel("y pos. [m]")
        #         plt.legend()
        #         plt.title(files[i])
        #         plt.tight_layout()
        #         export_pdf.savefig()
        #         plt.close()

        # for tstart in range(20):
        #     fig, ax = plt.subplots(figsize=(8,6))
        #     ax.set_prop_cycle((cycler(color=['r', 'g', 'b', 'm', 'c', 'k']) * cycler(linestyle=['-', '--', '-.'])))
        #     for i in range(len(files)):
        #         with h5py.File(files[i], "r") as fh:
        #             trng = (fh["t"][0,:]>=tstart)&(fh["t"][0,:]<tstart+1)
        #             plt.plot(fh["t"][0,trng],fh["xTM"][1,trng],label=labels[i])
        #     plt.legend(fontsize="8")
        #     plt.xlabel("time [s]")
        #     plt.ylabel("TM pos. [m]")
        #     plt.tight_layout()
        #     export_pdf.savefig()
        #     plt.close()
