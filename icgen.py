# Initial condition generator for GIZMO

import numpy as np
import h5py as h5py

def main():
    ds = np.loadtxt('out.txt')
    ngas = ds.shape[0]
    px_g = ds[:,0]
    py_g = ds[:,1]
    pz_g = ds[:,2]
    vx_g = np.zeros(ngas)
    vy_g = np.zeros(ngas)
    vz_g = np.zeros(ngas)
    uu_g = np.ones(ngas)
    mm_g = np.ones(ngas)/float(ngas)
    id_g = np.ones(ngas,dtype=np.int)
    
    file = h5py.File('ic.hdf5','w') 
    npart = np.array([ngas,0,0,0,0,0]) 
    h = file.create_group("Header")
    h.attrs['NumPart_ThisFile'] = npart
    h.attrs['NumPart_Total'] = npart
    h.attrs['NumPart_Total_HighWord'] = 0*npart
    h.attrs['MassTable'] = np.zeros(6)
    h.attrs['Time'] = 0.0
    h.attrs['Redshift'] = 0.0
    h.attrs['BoxSize'] = 1.0
    h.attrs['NumFilesPerSnapshot'] = 1
    h.attrs['Omega0'] = 1.0
    h.attrs['OmegaLambda'] = 0.0
    h.attrs['HubbleParam'] = 1.0
    h.attrs['Flag_Sfr'] = 0
    h.attrs['Flag_Cooling'] = 0
    h.attrs['Flag_StellarAge'] = 0
    h.attrs['Flag_Metals'] = 0
    h.attrs['Flag_Feedback'] = 0
    h.attrs['Flag_DoublePrecision'] = 0
    h.attrs['Flag_IC_Info'] = 0
    p = file.create_group("PartType0")
    q=np.zeros((ngas,3)); q[:,0]=px_g; q[:,1]=py_g; q[:,2]=pz_g;
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((ngas,3)); q[:,0]=vx_g; q[:,1]=vy_g; q[:,2]=vz_g;
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_g)
    p.create_dataset("Masses",data=mm_g)
    p.create_dataset("InternalEnergy",data=uu_g)
    
    file.close()

if __name__ == "__main__":
    main()
