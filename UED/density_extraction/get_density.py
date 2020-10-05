import sys, os, glob, time
import h5py
from copy import copy as copy
import numpy as np
import scipy as sp
import numpy.random as rnd
from multiprocessing import Pool
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#import jax as jax
#import jax.numpy as jnp

from params.N2O import data_parameters

qScale = 1.03
"""
def dists_(r1, r2):
  R = r1 - r2
  r = jnp.sqrt(jnp.sum(R**2))
  theta = jnp.arccos(R[2]/(r + 1e-20))
  phi   = jnp.arctan2(R[1],R[0])
  return jnp.array([r, theta, phi])

def calc_dists_vm(R):
  return jax.vmap(jax.vmap(dists_, in_axes=(None, 0)),
    in_axes=(0, None))(R,R)
"""

def calc_dists(R):
  r     = np.expand_dims(R, 1) - np.expand_dims(R, 0)
  dR    = np.sqrt(np.sum(r**2, axis=-1))
  theta = np.arccos(r[:,:,2]/(dR + 1e-20))
  phi   = np.arctan2(r[:,:,1], r[:,:,0])

  return np.concatenate([np.expand_dims(dR,-1),\
    np.expand_dims(theta, -1),\
    np.expand_dims(phi, -1)], axis=-1)

def calc_all_dists(R):
  r     = np.expand_dims(R, 2) - np.expand_dims(R, 1)
  dR    = np.sqrt(np.sum(r**2, axis=-1))
  theta = np.arccos(r[:,:,:,2]/(dR + 1e-20))
  phi   = np.arctan2(r[:,:,:,1], r[:,:,:,0])

  return np.concatenate([np.expand_dims(dR,-1),\
    np.expand_dims(theta, -1),\
    np.expand_dims(phi, -1)], axis=-1)



class density_extraction:

  def __init__(self, data_params):
    self.data_params = data_params
    self.atom_names = {
      "H" : "hydrogen",
      "C" : "carbon",
      "O" : "oxygen",
      "N" : "nitrogen",
      "I" : "iodine"
    }

    self.I = 1.
    self.do_rebin = False

    # De Broglie wavelength angs
    self.C_AU = 1./0.0072973525664
    self.eV_to_au = 0.0367493
    self.angs_to_au = 1e-10/5.291772108e-11 
    self.db_lambda = 2*np.pi*self.C_AU/\
        np.sqrt((self.data_params["elEnergy"]*self.eV_to_au + self.C_AU**2)**2\
        - (self.C_AU)**4) #au
    self.db_lambda /= self.angs_to_au  # angs
    self.k0 = 2*np.pi/self.db_lambda
    print("debrog", self.db_lambda)

    if not os.path.exists("./plots/{}".format(self.data_params["molecule"])):
      os.makedirs("./plots/{}".format(self.data_params["molecule"]))

    # Get initial geometry
    self.get_molecule_init_geo()

    # Get data
    self.data_or_sim = True
    if "simulate_data" in self.data_params:
      if self.data_params["simulate_data"]:
        self.data_or_sim = False
    self.get_data()

    # Get scattering amplitudes
    if self.data_params["experiment"] == "UED":
      self.get_scattering_amplitudes()
      self.evaluate_scattering_amplitudes()
    
    dist_inds1 = []
    dist_inds2 = []
    self.dist_sms_scat_amps = []
    for i, a1 in enumerate(self.atom_types):
      for j_, a2 in enumerate(self.atom_types[i+1:]):
        j = j_ + i+1
        dist_inds1.append(i)
        dist_inds2.append(j)
        self.dist_sms_scat_amps.append(
            self.scat_amps[a1]*self.scat_amps[a2]/self.atm_scat)
    self.dist_inds = (np.array(dist_inds1), np.array(dist_inds2))
    self.dist_sms_scat_amps = np.expand_dims(
        np.array(self.dist_sms_scat_amps), axis=0)


    # Simulate data if needed
    if not self.data_or_sim:
      self.input_data_coeffs_var *= self.data_params["simulate_data_scale"]
      self.simulate_data()

    # Prune data in time and dom
    self.prune_data()

    # Rebin Data
    self.do_rebin = False
    if "binning" in self.data_params:
      if self.data_params["binning"] is not None:
        self.do_rebin = True
        if self.data_params["binning"] == "log":
          Nrebin = np.cumsum(np.log2(np.arange(len(self.dom))+1).astype(int)) + 1
          Nrebin[1:] += 1
          prev = 0
          rebin_mat = []
          for Nb in Nrebin:
            rebin_mat.append(np.zeros_like(self.dom))
            rebin_mat[-1][prev:Nb] = 1./(np.min([Nb, len(self.dom)])-prev)
            prev = Nb
            if Nb >= len(self.dom):
              break
          self.rebin_mat = np.transpose(np.array(rebin_mat))
      self.eval_data_coeffs = np.matmul(self.data_coeffs, self.rebin_mat)
      self.eval_data_coeffs_var = np.matmul(self.data_coeffs_var, self.rebin_mat)
      self.eval_dom = np.matmul(self.dom, self.rebin_mat)
    else:
      self.eval_data_coeffs = self.data_coeffs
      self.eval_data_coeffs_var = self.data_coeffs_var
      self.eval_dom = self.dom
    self.eval_data_coeffs -= np.mean(self.eval_data_coeffs[:,:], axis=-1)



  def get_data(self):
    # Get the variance, measurement degree, legendre inds
    with h5py.File(self.data_params["data_fileName"], "r") as h5:
      self.input_data_coeffs_var = h5["fit_coeffs_cov"][:]
      self.input_data_coeffs_var = np.ones_like(self.input_data_coeffs_var)
      self.input_data_coeffs_var *= np.expand_dims(self.input_data_coeffs_var[:,70,:,:], -1)
      self.data_lg = h5["legendre_inds"][:]
      self.dom = h5["fit_axis"][:]*qScale
      if self.data_or_sim:
        self.input_data_coeffs = h5["fit_coeffs"][:]

  
  def simulate_data(self):
    print("SIMMING", self.data_lg, self.atm_scat[70])
    molecule = self.setup_calculations(skip_scale=True, plot=False)
    self.input_data_coeffs = self.calculate_coeffs(molecule)
    tmp = copy(self.input_data_coeffs)
    shift = rnd.normal(size=self.input_data_coeffs.shape)
    shift *= np.sqrt(self.input_data_coeffs_var[:,:,0,0])
    if not self.data_params["isMS"]:
      self.input_data_coeffs *= self.atm_scat
    # TODO debugging
    #self.input_data_coeffs += shift
    self.input_data_coeffs = np.expand_dims(self.input_data_coeffs, -1)
    print(self.input_data_coeffs.shape, self.input_data_coeffs_var.shape)
    plt.plot(self.dom[70:], tmp[1,70:])
    plt.errorbar(self.dom[70:], self.input_data_coeffs[1,70:,0], np.sqrt(self.input_data_coeffs_var[1,70:,0,0]))
    plt.savefig("init_werr.png")
    plt.close()



  def prune_data(self):
    # Devide by atomic scattering
    if not self.data_params["isMS"]:
      atm_scat_ = np.reshape(self.atm_scat, (1,len(self.atm_scat),1))
      self.input_data_coeffs /= atm_scat_
      atm_scat_ = np.expand_dims(atm_scat_, -1)
      print("wtf", self.input_data_coeffs_var.shape, atm_scat_.shape)
      self.input_data_coeffs_var /= atm_scat_**2
    
    # Prune the list of legendre projections
    keep_inds_lg = np.arange(self.input_data_coeffs_var.shape[0])
    if "fit_bases" in self.data_params:
      keep_inds_lg = []
      temp_data_lg = []
      for i,lg in enumerate(self.data_lg):
        if lg in self.data_params["fit_bases"]:
          keep_inds_lg.append(i)
          temp_data_lg.append(lg)
          print("INFO: Will fit Legendre: {}".format(lg)) 
      keep_inds_lg = np.array(keep_inds_lg)
      self.data_lg = np.array(temp_data_lg)


    self.data_coeffs = self.input_data_coeffs[keep_inds_lg,:,0]
    self.data_coeffs_var = self.input_data_coeffs_var[keep_inds_lg,:,0,0]
        
        
    # Prune dom axis
    mask = np.ones_like(self.dom).astype(bool)
    if "fit_range" in self.data_params:
      mask[self.dom<self.data_params["fit_range"][0]] = False
      mask[self.dom>self.data_params["fit_range"][1]] = False
    self.dom = self.dom[mask]
    self.data_coeffs = self.data_coeffs[:,mask]
    self.data_coeffs_var = self.data_coeffs_var[:,mask]
    self.atm_scat = self.atm_scat[mask]
    for atm in self.scat_amps.keys():
      self.scat_amps[atm] = self.scat_amps[atm][mask] 
    self.dist_sms_scat_amps = self.dist_sms_scat_amps[:,:,mask]
    
    for lg in range(len(self.data_lg)):
      handles = []
      handles.append(plt.errorbar(self.dom, self.data_coeffs[lg,:],\
          np.sqrt(self.data_coeffs_var[lg,:]),\
          label="legendre {}".format(self.data_lg[lg])))
      plt.legend(handles=handles)
      plt.savefig("./plots/{}/data_coeffs_lg-{}.png".format(
        self.data_params["molecule"], self.data_lg[lg]))
      plt.close()


  def evaluate_scattering_amplitudes(self):
    self.atm_scat = np.zeros_like(self.dom)
    self.scat_amps = {}
    for atm in self.atom_types:
      if atm not in self.scat_amps:
        self.scat_amps[atm] = self.scat_amps_interp[atm](self.dom)
      self.atm_scat += self.scat_amps[atm]**2


  
  def get_molecule_init_geo(self):
    if not os.path.exists(self.data_params["init_geo_xyz"]):
      print("Cannot find xyz file: " + self.data_params["init_geo_xyz"])
      sys.exit(1)

    self.atom_types      = []
    self.atom_positions  = []
    with open(self.data_params["init_geo_xyz"]) as file:
      for i,ln in enumerate(file):
        if i == 0:
          Natoms = int(ln)
        elif i > 1:
          vals = ln.split()
          print(vals)
          self.atom_types.append(vals[0])
          pos = [float(x) for x in vals[1:]]
          self.atom_positions.append(np.array([pos]))

    self.atom_positions = np.concatenate(self.atom_positions, axis=0)

  
  def get_scattering_amplitudes(self):

    self.scat_amps_interp = {}
    for atm in self.atom_types:
      if atm in self.scat_amps_interp:
        continue

      angStr = []
      sctStr = []
      fName = os.path.join(self.data_params["scat_amps_path"],
          self.atom_names[atm] + "_dcs.dat")
      with open(fName, 'r') as inpFile:
        ind=0
        for line in inpFile:
          if ind < 31:
            ind += 1
            continue

          angStr.append(line[2:11])
          sctStr.append(line[39:50])

      angs = np.array(angStr).astype(np.float64)*np.pi/180
      q = 4*np.pi*np.sin(angs/2.)/self.db_lambda
      scts = np.sqrt(np.array(sctStr).astype(np.float64))

      self.scat_amps_interp[atm] = interp1d(q, scts, 'cubic')

    
  def calculate_coeffs(self, R):
    all_dists = calc_dists(R)
    dists = all_dists[self.dist_inds]

    C = np.complex(0,1)**self.calc_data_lg*np.sqrt(4*np.pi*(2*self.calc_data_lg + 1))
    J = sp.special.spherical_jn(self.calc_data_lg, 
        self.calc_dom*np.expand_dims(dists[:,0], axis=-1))
    Y = sp.special.sph_harm(0, self.calc_data_lg,
        np.expand_dims(np.expand_dims(dists[:,2], axis=0), axis=-1),
        np.expand_dims(np.expand_dims(dists[:,1], axis=0), axis=-1))

    calc_coeffs = np.sum(np.real(self.dist_sms_scat_amps*C*J*Y), axis=1)
    if self.do_rebin:
      calc_coeffs = np.matmul(calc_coeffs, self.rebin_mat)
    calc_coeffs -= np.expand_dims(np.mean(calc_coeffs[:,:], axis=-1), -1)
    calc_coeffs *= self.I

    return calc_coeffs


  def calculate_log_prob(self, R, n=0):

    calc_coeffs = self.calculate_coeffs(R)

    """
    X = np.concatenate(
        [np.expand_dims(calc_coeffs, -1), np.ones(calc_coeffs.shape + (1,))], -1)
    th = np.einsum('abi,ai->ab',
        np.linalg.inv(np.einsum('aib,ai,aic->abc', X, 1./self.data_coeffs_var, X)),
        np.einsum('aib,ai,ai->ab', X, 1./self.data_coeffs_var, self.data_coeffs))

    calc_coeffs *= th[:,0]#scale
    calc_coeffs += th[:,1]
    """
    if n > 500000:
      plt.errorbar(self.dom, self.eval_data_coeffs[0,:], np.sqrt(self.eval_data_coeffs_var[0,:]))
      plt.plot(self.dom, calc_coeffs[0,:])
      plt.plot(self.dom, 3e-5*((self.eval_data_coeffs - calc_coeffs)**2\
          /self.eval_data_coeffs_var)[0,:])
      plt.text(8,0.13,"{}".format(R[0,2], R[1]))
      plt.text(6,0.11,"{}".format(R[2]))
      plt.savefig("test_fit_{}.png".format(n))
      plt.close()

    prob = np.sum(-0.5*(self.eval_data_coeffs - calc_coeffs)**2/self.eval_data_coeffs_var)
    #    + np.log(1/np.sqrt(self.data_coeffs_var)))
   
    return prob
 


  def setup_calculations(self, skip_scale=False, plot=True):
    self.calc_data_lg = np.expand_dims(np.expand_dims(self.data_lg, axis=-1), axis=-1)
    self.calc_dom = np.expand_dims(self.dom, axis=0)
   
    # Calculate Scale Factor
    if not skip_scale:
      self.I = 1.
      calc_coeffs = self.calculate_coeffs(self.atom_positions)
      if "scale" in self.data_params:
        if self.data_params["scale"] is None:
          self.I = np.sum(calc_coeffs/self.eval_data_coeffs_var*self.eval_data_coeffs)\
            /np.sum(calc_coeffs/self.eval_data_coeffs_var*calc_coeffs)
        else:
          self.I = self.data_params["scale"]
      else:
        print(calc_coeffs.shape, self.eval_data_coeffs.shape, self.eval_data_coeffs_var.shape)
        plt.figure()
        plt.plot(self.eval_dom, calc_coeffs[0,:])
        plt.savefig("test.png")
        plt.close()
        self.I = np.sum(calc_coeffs/self.eval_data_coeffs_var*self.eval_data_coeffs)\
          /np.sum(calc_coeffs/self.eval_data_coeffs_var*calc_coeffs)
      print("INFO: I coefficient: {}".format(self.I))
      self.init_fit = self.I*calc_coeffs

    #TODO debugging
    #self.I = 1.
    if plot:
      plt.errorbar(self.eval_dom, self.eval_data_coeffs[0,:], np.sqrt(self.eval_data_coeffs_var[0,:]))
      plt.plot(self.eval_dom, self.I*calc_coeffs[0,:])
      plt.savefig("init_scale_fit_{}-{}.png".format(0,0))#self.data_params["fit_range"][0], self.data_params["fit_range"][1]))
      plt.close()


    return self.atom_positions


  def get_molecule_perturber(self):

    shape = self.atom_positions.shape

    if self.data_params["molecule"] == "N2O":
      def perturb_molecule(prev_molecule):
        perturbation = np.random.uniform(-1,1, 3)
        perturbation *= self.data_params["perturb_range"]

        rNO = np.linalg.norm(prev_molecule[2,:])
        angle = np.arccos(prev_molecule[0,2]*prev_molecule[2,2]\
            /(rNO*prev_molecule[0,2]))
        rNO += perturbation[1] 
        angle += perturbation[2] 
        
        molecule = np.zeros_like(prev_molecule)
        molecule[0,2] = prev_molecule[0,2] + perturbation[0]
        molecule[2,0] = rNO*np.sin(angle)
        molecule[2,2] = rNO*np.cos(angle)

        return molecule

    else:
      def perturb_molecule(molecule):
        perturbation = np.random.uniform(
            -1*self.data_params["perturb_range"],
            self.data_params["perturb_range"],
            shape)

        return molecule + perturbation

    return perturb_molecule


  def save_results(self, probs, distributions):

    fName = os.path.join("output", self.data_params["molecule"],
        "log_probabilities-{}.npy".format(len(probs))) 
    with open(fName, "wb") as file:
      np.save(file, np.array(probs))

    fName = os.path.join("output", self.data_params["molecule"],
        "molecules-{}.npy".format(len(probs))) 
    with open(fName, "wb") as file:
      np.save(file, np.concatenate(distributions))

    fName = os.path.join("output", self.data_params["molecule"],
        "order-{}.npy".format(len(probs))) 
    with open(fName, "wb") as file:
      np.save(file, np.arange(len(probs)))


  def get_results(self):
    files = glob.glob(
        os.path.join("output", self.data_params["molecule"],"molecules*"))
    N = 0
    for fl in files:
      n = int(fl[fl.find("-")+1:-4])
      if n > N:
        N = n
    
    if len(files) > 0:
      fName = os.path.join("output", self.data_params["molecule"],
          "log_probabilities-{}.npy".format(N)) 
      with open(fName, "rb") as file:
        probs = np.load(file)

      fName = os.path.join("output", self.data_params["molecule"],
          "molecules-{}.npy".format(N)) 
      with open(fName, "rb") as file:
        distributions = np.load(file)

      fName = os.path.join("output", self.data_params["molecule"],
          "order-{}.npy".format(N)) 
      with open(fName, "rb") as file:
        order = np.load(file)


      probs = probs[order].tolist()
      distributions = [np.expand_dims(x, axis=0) for x in distributions[order]]

      return probs, distributions

    else:
      return [], []


def main(fit_range=None, simulate_data_scale=None):


  data_params = copy(data_parameters)
  if fit_range is not None:
    data_params["fit_range"] = fit_range
  if simulate_data_scale is not None:
    data_params["simulate_data_scale"] = simulate_data_scale

  print(data_params["fit_range"])
  extraction = density_extraction(data_params)
  molecule = extraction.setup_calculations()
  
  perturb_molecule = extraction.get_molecule_perturber()
  current_log_prob = extraction.calculate_log_prob(molecule)

  #molecule_log_prob,  molecule_distributions = extraction.get_results()
  
  tic = time.time() 
  molecule_distributions = []
  molecule_log_prob = []
  if len(molecule_distributions):
    molecule = molecule_distributions[-1][0]
  while len(molecule_distributions) < 60000:
    new_molecule = perturb_molecule(molecule)
    new_log_prob = extraction.calculate_log_prob(
        new_molecule, len(molecule_distributions))

    if np.exp(new_log_prob - current_log_prob) > np.random.uniform():
      
      molecule_distributions.append(np.expand_dims(copy(new_molecule), axis=0))
      molecule_log_prob.append(new_log_prob)
      molecule = copy(new_molecule)
      current_log_prob = new_log_prob
    #print(len(molecule_distributions), new_log_prob, current_log_prob, new_log_prob - current_log_prob, np.exp(new_log_prob - current_log_prob))

    if len(molecule_distributions) % 1000 == 0 and time.time() - tic > 1:
      print("Finished {}: time since previous {}".format(
          len(molecule_distributions), time.time() - tic))
      tic = time.time()

  extraction.save_results(molecule_log_prob, molecule_distributions)

  #TODO
  """
  #####  Plotting  #####
  Rs = np.concatenate(molecule_distributions, 0)[30000:]
  all_dists = calc_all_dists(Rs)
  angles = np.arccos(np.sum((Rs[:,0,:] - Rs[:,1,:])*(Rs[:,2,:] - Rs[:,1,:]), -1)/(all_dists[:,1,0,0]*all_dists[:,2,1,0]))
  fig, ax = plt.subplots(1,3,figsize=(15,5))
  h, x, y, _ = ax[0].hist2d(all_dists[:,0,1,0], all_dists[:,1,2,0], bins=[100, 100])
  ax[0].set_xlabel("NN Distance [A]")
  ax[0].set_ylabel("NO Distance [A]")
  y, x = np.histogram(angles/np.pi, bins=100)
  ax[1].plot(x[:-1], y)
  ax[1].set_xlim([0.85,1])
  ax[1].set_xlabel("N-N-O Angle [rad/pi]")
  ax[2].errorbar(extraction.eval_dom, extraction.eval_data_coeffs[0,:],
      np.sqrt(extraction.eval_data_coeffs_var[0,:]))
  ax[2].plot(extraction.eval_dom, extraction.init_fit[0,:])
  ax[2].set_xlabel("Q")
  plt.tight_layout()
  fig.savefig("./plots/N2O/constant_error/results_range-{}-{}_scale-{}.png".format(
      data_params["fit_range"][0], data_params["fit_range"][1],
      data_params["simulate_data_scale"]))
  """


def loop_main(i,j):

  fit_ranges = [[3,5], [3,5.5], [3,6], [3,6.5], [3,6.5], [3,7], [3,7.5], [3,8]]
  scales = np.array([0.01, 0.05, 0.1, 0.5, 1, 10])*1e2/(1540**2)
  main(fit_range=fit_ranges[i], simulate_data_scale=scales[j])

def mp_loop_main(i):

  i,j = i//8, i%8
  fit_ranges = [[3,5], [3,5.5], [3,6], [3,6.5], [3,6.5], [3,7], [3,7.5], [3,8]]
  scales = np.array([0.01, 0.05, 0.1, 0.5, 1, 10])*1e2/(1540**2)
  print("INP", i, j, fit_ranges[i], scales[j])
  main(fit_range=fit_ranges[i], simulate_data_scale=scales[j])



if __name__ == "__main__":
  main(fit_range=[3,6], simulate_data_scale=0.01*1e2/(1540**2))

  #for i in range(8):
  #  for j in range(6):
  #    print("Main Loop",i,j)
  #    loop_main(i,j)

  #pool = Pool(processes=6)
  #results = pool.map(mp_loop_main, range(8*6))

