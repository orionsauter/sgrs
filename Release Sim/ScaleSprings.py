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
    common_args = ['config','nominal.xml','savestats','0','savedata','1','run_name','"'+("_".join(map(str, args)))+'"','kG','1.1e6','kT','1.1e6']
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

    rt1 = np.pad(fh["xT1"][:]+sTM/2,((1,1),(0,0))).T # RT1 position
    rt2 = np.pad(fh["xT2"][:]-sTM/2,((1,1),(0,0))).T # RT2 position
    # vt1 = np.pad(fh["vT1"][:],((1,1),(0,0))) # RT1 speed
    # vt2 = np.pad(fh["vT2"][:],((1,1),(0,0))) # RT2 speed

    rg1 = np.pad(fh["xG1"][:]+sTM/2,((1,1),(0,0))).T # GF1 position
    rg2 = np.pad(fh["xG2"][:]-sTM/2,((1,1),(0,0))).T # GF2 position
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
    contact = np.vstack([cG_RT1,cG_RT2,cG_TM1,cG_TM2,cRT_TM1,cRT_TM2,fh["bond1"][0,:],fh["bond2"][0,:]]).T
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
                t = fhs[i]["t"][0,:]
                # contact = GetContact(fh, sTM, Lt)
                # change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
                axs[dim].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xTM"][dim,::dsamp],"-",label=labels[i])
                # ylim = axs[dim].get_ylim()
                # axs[dim].vlines(t[change], *ylim, linestyles="dotted", colors="black")
                # axs[dim].set_ylim(ylim)
            axs[dim].set_ylabel(f"TM {dimnames[dim]} pos. [m]")
        plt.legend()
        axs[2].set_xlabel("time [s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

        fig, axs = plt.subplots(figsize=(8,8),nrows=2,sharex=True)
        axs[0].set_prop_cycle((cycler(color=['r', 'b', 'g', 'orange']) + cycler(linestyle=['-', '--','-.', ':'])))
        axs[1].set_prop_cycle((cycler(color=['r', 'b', 'g', 'orange']) + cycler(linestyle=['-', '--','-.', ':'])))
        for i in range(len(files)):
            t = fhs[i]["t"][0,:]
            # contact = GetContact(fh, sTM, Lt)
            # change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
            # line, = axs[0].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xG1"][0,::dsamp],"-",label=labels[i])
            axs[0].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xT1"][0,::dsamp],"--",label=labels[i])
            # ylim = axs[0].get_ylim()
            # axs[0].vlines(t[change], *ylim, linestyles="dotted", colors=line.get_color())
            # axs[0].set_ylim(ylim)
            # line, = axs[1].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xG2"][0,::dsamp],"-",label=labels[i])
            axs[1].plot(fhs[i]["t"][0,::dsamp],fhs[i]["xT2"][0,::dsamp],"--",label=labels[i])
            # ylim = axs[1].get_ylim()
            # axs[1].vlines(t[change], *ylim, linestyles="dotted", colors=line.get_color())
            # axs[1].set_ylim(ylim)
        axs[0].set_ylabel(f"GF/RT 1 pos. [m]")
        axs[1].set_ylabel(f"GF/RT 2 pos. [m]")
        plt.legend()
        axs[1].set_xlabel("time [s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

        fig, axs = plt.subplots(figsize=(8,8),nrows=3,sharex=True)
        for dim in range(3):
            t_ex = fhs[0]["t"][0,:]
            t_ap = fhs[1]["t"][0,:]
            x_ex = fhs[0]["xTM"][dim,:]
            x_ap = sci.InterpolatedUnivariateSpline(t_ap,fhs[1]["xTM"][dim,:])
            axs[dim].plot(t_ex[::dsamp],x_ap(t_ex[::dsamp])-x_ex[::dsamp],"-",label=labels[i])
            axs[dim].set_ylabel(f"TM {dimnames[dim]} pos. err. [m]")
        # plt.legend()
        axs[2].set_xlabel("time [s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()


if __name__ == "__main__":

    logging.basicConfig(
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO,
        datefmt='%Y-%m-%d %H:%M:%S')
    
    neng = 10 # Number of Matlab instances
    nval = 3 # Number of vals for each var
    # nsamp = 10 # Number of runs
    nclust = 8
    ncore = 4
    client = AllocWorkers(nclust, ncore)
    fs = 3e6
    dt = 1.0/fs

    sTM = 30e-3
    sEH = 32e-3
    gap = (sEH-sTM)/2
    Lt = 0.5e-4
    tr = 4e-6

    # logging.info("Generating var sets.")
    # vars = []
    # # vars = ["vrsy", "ts1", "ts2", "tp", "tc"]
    # # rngs = [[2e-6,6e-6], [0,10e-5], [0,10e-5], [0,2], [5,7]]#, [0,1e-4], [0,1e-4], [0,1e-2], [0,1e-2]]
    vars = ["posTol", "velTol"]
    rngs = [[-8,-6], [-6.5,-6]]
    # arrs = [np.linspace(*rng, nval) for rng in rngs]
    arrs = [np.logspace(*rng, nval) for rng in rngs]
    cols = np.array([col.flatten() for col in np.meshgrid(*arrs)]).T
    df = pd.DataFrame(cols, columns=vars).astype(str)

    # vecs = ["xTM0"]
    # rngs = [[0,1e-4], [0,0], [0,1e-4]]
    # arrs = [np.linspace(*rng, nval) for rng in rngs]
    # vec_cols = np.array([col.flatten() for col in np.meshgrid(*arrs)]).T
    # vec_df = pd.DataFrame(np.array([pd.DataFrame(vec_cols[:,i:(i+3)], columns=["x","y","z"], dtype="str").agg(
    #     lambda x: "[" + ", ".join(x) + "]", axis=1) for i in range(0,vec_cols.shape[1],3)]).T, columns=vecs)
    # df = df.join(vec_df, how="outer").fillna(0)

    # df = df.sample(nsamp, random_state=42)

    # df = pd.DataFrame({"kG":[5e7,5e7],"kT":[5e7,5e7],"posTol":[1e-8,3e-8,1e-7],"velTol":[1e-8,3e-8,1e-7]})
    arg_list = [sum(df.loc[i].to_dict().items(),()) for i in df.index]
    res = client.map(RunSim, arg_list)
    logging.info("Awaiting results.")
    res_list = client.gather(res)
    for res in res_list:
        logging.info(res)
    res_df = pd.concat([df,pd.concat(res_list,axis=1).T],axis=1)
    # # df[res_df.index,res_df.columns] = res_df
    # logging.info(res_df)
    res_df.to_hdf("Gaps.h5","gaps")
    # res = RunSim(arg_list[-1])
    exit()

    files = glob("Nominal_*/*.mat")
    # files = [glob(f"{('_'.join(map(str, args)))}_*/*.mat")[0] for args in arg_list] + ["WideK/kG_500000.0_kT_500000.0_202412412480/TMrel_202412412480.mat"]
    labels = ["exact", "est.", "old"]
    PlotPositions(files, labels)
    exit()
    # labels = [f"{row['kG']:.0f}, {row['kT']:.0f}" for _, row in df.iterrows()]
    # files = glob(f"*K/*/*.mat")
    # kG = [float(file.split("_")[1]) for file in files]
    # kT = [float(file.split("_")[3]) for file in files]
    # labels = [f"{kG[i]:.0f}, {kT[i]:.0f}" for i in range(len(kG))]
    with PdfPages("SpringScale.pdf") as export_pdf:
        # for dim in range(3):
        # dim = 1
        # fig, ax = plt.subplots(figsize=(8,6))
        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         t = fh["t"][0,:]
        #         amp = fh["vTM"][dim,:]
        #         # Get sign changes & find max in each pos region
        #         idx = np.hstack([0,np.where(np.sign(amp[:-1]) != np.sign(amp[1:]))[0] + 1,-1])
        #         mxs = np.array([idx[i]+np.argmax(amp[idx[i]:idx[i+1]]) for i in range(len(idx)-1) if amp[(idx[i]+idx[i+1])//2]>0])
        #         mts = t[mxs]
        #         f = 1.0/np.diff(mts)
        #         plt.semilogy((mts[:-1]+mts[1:])/2,f,"-",label=labels[i])
        #         contact = GetContact(fh, sTM, Lt)
        #         change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
        #         xlim = plt.xlim()
        #         ylim = plt.ylim()
        #         plt.vlines(t[change], *ylim, linestyles="dotted", colors="black")
        #         plt.xlim(xlim)
        #         plt.ylim(ylim)
        # w0 = np.sqrt(1.317e4/0.5378)/4
        # plt.hlines(np.sqrt(w0*w0-2*(4*1.45547)**2), *xlim, linestyles="dotted", colors="black")
        # plt.xlim(xlim)
        # plt.legend(fontsize="8")
        # plt.xlabel("time [s]")
        # plt.ylabel("TM vel. freq. [Hz]")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # slopes = np.zeros((len(files)))
        # intcps = np.zeros((len(files)))
        # for i in range(len(files)):
        #     fig, axs = plt.subplots(figsize=(8,8),nrows=2,sharex=True)
        #     # axs[0].set_prop_cycle((cycler(color=['r', 'g', 'b', 'm', 'c', 'k']) * cycler(linestyle=['-', '--', '-.'])))
        #     # axs[1].set_prop_cycle((cycler(color=['r', 'g', 'b', 'm', 'c', 'k']) * cycler(linestyle=['-', '--', '-.'])))
        #     with h5py.File(files[i], "r") as fh:
        #         t = fh["t"][0,:]
        #         GF1 = fh["xG1"][0,:]
        #         GF2 = fh["xG2"][0,:]
        #         RT1 = fh["xT1"][0,:]
        #         RT2 = fh["xT2"][0,:]
        #         contact = GetContact(fh, sTM, Lt)
        #         change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
        #         trng = (t>5)&(t<15)
        #         p = np.polyfit(t[trng], RT1[trng], 1)
        #         # m, c = np.linalg.lstsq(np.vstack([t[trng], np.ones(len(t[trng]))]).T, RT1[trng])[0]
        #         slopes[i] = p[0]
        #         intcps[i] = p[1]
        #         estm = 10**(-0.99096642*np.log10(kT[i])-1.640234)
        #         estb = 10**(-0.98271567*np.log10(kT[i])-0.41746246)
        #         estR1 = t*estm+estb
        #         # t = t[t>17.5]
        #         dsamp = 500
        #         axs[0].plot(t[::dsamp],RT1[::dsamp],".",label=labels[i])
        #         # axs[1].plot(t[::dsamp],RT2[::dsamp],".",label=labels[i])
        #         axs[0].plot(t[::dsamp],estm*t[::dsamp]-estb,"--",label=labels[i])
        #         # axs[1].plot(t[::dsamp],-estm*t[::dsamp]+estb,"--",label=labels[i])
        #         axs[1].semilogy(t[::dsamp],np.abs(fh["vTM"][1,::dsamp]),"-")
        #         axs[0].vlines(t[change], *axs[0].get_ylim(), linestyles="dotted", colors="black")
        #         axs[1].vlines(t[change], *axs[1].get_ylim(), linestyles="dotted", colors="black")
        #         axs[0].set_ylim((-np.abs(RT1[0]),np.abs(RT1[0])))
        #         # axs[1].set_ylim((-np.abs(RT2[0]),np.abs(RT2[0])))
        #         # axs[0].plot(t,RT1-RT2,label=labels[i])
        #         # axs[1].plot(t,GF1-RT1,label=labels[i])
        #         # axs[1].plot(t,RT2-GF2,label=labels[i])
        #     # plt.legend(fontsize="8")
        #     axs[1].set_xlabel("time [s]")
        #     axs[0].set_ylabel("RT1 pos [m]")
        #     axs[1].set_ylabel("TM vel. [m/s]")
        #     # axs[1].set_ylabel("RT2 poss [m]")
        #     plt.tight_layout()
        #     export_pdf.savefig()
        #     plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # plt.loglog(kT,slopes,".")
        # p = np.polyfit(np.log10(kT),np.log10(slopes),1)
        # logging.info(f"Slope fit: m = {p[0]:.8g}, b = {p[1]:.8g}")
        # plt.loglog(kT, 10**(p[0]*np.log10(kT)+p[1]),"-")
        # plt.xlabel("kT")
        # plt.ylabel("RT1 pos. slope")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # plt.loglog(kT,-intcps,".")
        # p = np.polyfit(np.log10(kT),np.log10(-intcps),1)
        # logging.info(f"Intcp fit: m = {p[0]:.8g}, b = {p[1]:.8g}")
        # plt.loglog(kT, 10**(p[0]*np.log10(kT)+p[1]),"-")
        # plt.xlabel("kT")
        # plt.ylabel("RT1 pos. slope")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # plt.scatter(kG,slopes,c=np.log10(kT))
        # plt.xlabel("kG")
        # plt.ylabel("RT1 pos. slope")
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.colorbar(label="log kT")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # plt.scatter(kG,-np.array(intcps),c=np.log10(kT))
        # plt.xlabel("kG")
        # plt.ylabel("RT1 pos. intcp.")
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.colorbar(label="log kT")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

            # fig, ax = plt.subplots(figsize=(8,6))
            # for i in range(len(files)):
            #     with h5py.File(files[i], "r") as fh:
            #         t = fh["t"][0,:]
            #         amp = fh["vTM"][dim,:]
            #         amp = amp[(t>3.65)&(t<7.0)]
            #         t = t[(t>3.65)&(t<7.0)]
            #         # Get sign changes & find max in each pos region
            #         idx = np.hstack([0,np.where(np.sign(amp[:-1]) != np.sign(amp[1:]))[0] + 1,-1])
            #         mxs = np.array([idx[i]+np.argmax(amp[idx[i]:idx[i+1]]) for i in range(len(idx)-1) if amp[(idx[i]+idx[i+1])//2]>0])
            #         mts = t[mxs]
            #         f = 1.0/np.diff(mts)
            #         plt.semilogy((mts[:-1]+mts[1:])/2,f,"-",label=labels[i])
            # xlim = plt.xlim()
            # w0 = np.sqrt(1.317e4/0.5378)/4
            # plt.hlines(np.sqrt(w0*w0-2*(4*1.45547)**2), *xlim, linestyles="dotted", colors="black")
            # plt.xlim(xlim)
            # plt.ylim([30,40])
            # plt.legend()
            # plt.xlabel("time [s]")
            # plt.ylabel("TM vel. freq. [Hz]")
            # plt.tight_layout()
            # export_pdf.savefig()
            # plt.close()

            # fig, ax = plt.subplots(figsize=(8,6))
            # for i in range(len(files)):
            #     with h5py.File(files[i], "r") as fh:
            #         t = fh["t"][0,:]
            #         amp = fh["vTM"][dim,:]
            #         # Get sign changes & find max in each pos region
            #         idx = np.hstack([0,np.where(np.sign(amp[:-1]) != np.sign(amp[1:]))[0] + 1,-1])
            #         mxs = np.array([idx[i]+np.argmax(amp[idx[i]:idx[i+1]]) for i in range(len(idx)-1) if amp[(idx[i]+idx[i+1])//2]>0])
            #         # mxs = sg.argrelmax(amp)
            #         mts = t[mxs]
            #         line, = plt.semilogy(mts,amp[mxs],"-",label=labels[i])
            #         contact = GetContact(fh, sTM, Lt)
            #         change = np.where(np.any(np.diff(contact,axis=0)!=0,axis=1))[0]
            #         unbond = np.where(np.any(np.diff(contact[:,-2:],axis=0)!=0,axis=1))[0][0]
            #         # Breaking bond adds impulse to vel. Skip ahead
            #         impendt = t[unbond] + 1e-2
            #         impend = np.argmin(np.abs(t-impendt))
            #         unbt = np.arange(t[impend],t[-1],0.1)
            #         unbv = amp[impend]
            #         xlim = plt.xlim()
            #         ylim = plt.ylim()
            #         plt.vlines(t[change], *ylim, linestyles="dotted", colors=line.get_color())
            #         plt.plot(unbt,unbv*np.exp(-4*1.45547*(unbt-unbt[0])), linestyle="dotted", color=line.get_color())
            #         plt.xlim(xlim)
            #         plt.ylim(ylim)
            # plt.legend()
            # plt.xlabel("time [s]")
            # plt.ylabel("TM vel. amp. [m/s]")
            # plt.tight_layout()
            # export_pdf.savefig()
            # plt.close()

        rtsep1 = np.zeros((len(files),))
        rtsep2 = np.zeros((len(files),))
        rtvel = np.zeros((len(files),))
        fig, ax = plt.subplots(figsize=(8,6))
        for i in range(len(files)):
            with h5py.File(files[i], "r") as fh:
                t = fh["t"][0,:]
                contact = GetContact(fh, sTM, Lt)
                change = np.where(np.any(np.diff(contact[:,4:6],axis=0)!=0,axis=1))[0]
                rtsep1[i] = np.min(t[change])
                rtsep2[i] = np.max(t[change])
                tstart = rtsep2[i]
                trng = (t>=tstart-5e-4)&(t<tstart+3e-4)
                rtvel[i] = np.abs((fh["vTM"][1,trng])[-1])
                plt.semilogy(t[trng]-tstart,np.abs(fh["vTM"][1,trng]),"-",label=labels[i])
        plt.legend(fontsize="8")
        plt.xlabel("time after RT sep. [s]")
        plt.ylabel("TM y velocity [m/s]")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

        # accrise = np.zeros((len(files),))
        # fig, ax = plt.subplots(figsize=(8,6))
        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         t = fh["t"][0,:]
        #         a = np.diff(fh["vTM"][1,:])/np.diff(t)
        #         accrise[i] = np.min(t[:-1][(t[:-1]>17.3)&(np.abs(a)>1e-4)])
        #         # endi = np.max(np.where(a!=0)[0])
        #         # plt.semilogy(t[starti[i]:-1],a[starti[i]:],"-",label=labels[i])
        #         # plt.plot(t[:-1][t[:-1]>rtsep1[i]-2],a[t[:-1]>rtsep1[i]-2],"-",label=labels[i])
        #         plt.plot(t[:-1][(t[:-1]>17.3)&(t[:-1]<17.9)],a[(t[:-1]>17.3)&(t[:-1]<17.9)],"-",label=labels[i])
        # plt.legend()
        # plt.xlabel("time [s]")
        # plt.ylabel("TM y accel. [m/s^2]")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         t = fh["t"][0,:]
        #         # trng = (t>=rtsep2[i]-1)&(t<=rtsep2[i])
        #         yvel = np.abs(fh["vTM"][1,:])
        #         hist, bins = np.histogram(np.log10(yvel[yvel>0]), 200)
        #         plt.step(10**bins[:-1],hist/np.sum(hist),"-",label=labels[i])
        # plt.legend(fontsize="8")
        # plt.xlabel("TM y vel. [m/s]")
        # plt.ylabel("density")
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # kG = [arg_list[i][1] for i in range(len(files))]
        # kT = [arg_list[i][3] for i in range(len(files))]
        tdiff = rtsep2-rtsep1
        fig, ax = plt.subplots(figsize=(8,6))
        plt.loglog(kT,rtsep1,'.',label="1st")
        plt.loglog(kT,rtsep2,'.',label="2nd")
        plt.legend(fontsize="8")
        plt.xlabel("kT")
        plt.ylabel("RT sep. time")
        plt.tight_layout()
        export_pdf.savefig()
        plt.close()

        # fig, ax = plt.subplots(figsize=(8,6))
        # plt.scatter(kT,kG,c=accrise)
        # plt.colorbar(label="acc. rise time")
        # plt.xlabel("kT")
        # plt.ylabel("kG")
        # plt.xscale("log")
        # plt.yscale("log")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

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

        # steps = []
        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         steps += [np.diff(fh["t"][0,:])]
        # fig, ax = plt.subplots(figsize=(8,6))
        # _, bins = np.histogram(np.log10(np.hstack(steps)), 100)
        # ax.hist(steps, 10**bins, histtype="step", stacked=False, fill=False, label=labels)
        # ax.set_xscale("log")
        # ax.set_yscale("log")
        # plt.legend()
        # plt.xlabel("timestep [s]")
        # plt.ylabel("count")
        # plt.tight_layout()
        # export_pdf.savefig()
        # plt.close()

        # for i in range(len(files)):
        #     with h5py.File(files[i], "r") as fh:
        #         # func = sci.CubicSpline(fh["t"][0,:], fh["vTM"][:])
        #         func = sci.UnivariateSpline(fh["t"][0,:], fh["vTM"][1,:])
        #         tunif = np.arange(fh["t"][0,0],fh["t"][0,-1]+dt,dt)
        #         f, t, Sxx = sg.spectrogram(func(tunif), fs, nperseg=int(fs*1e-3))
        #         fig, ax = plt.subplots(figsize=(8,6))
        #         plt.pcolormesh(t, f, np.sqrt(Sxx))
        #         plt.xlabel("time [s]")
        #         plt.ylabel("freq. [Hz]")
        #         plt.ylim([1e5,2e6])
        #         plt.colorbar(label="y vel. ASD [m/s/sqrt(Hz)]")
        #         plt.tight_layout()
        #         export_pdf.savefig()
        #         plt.close()
