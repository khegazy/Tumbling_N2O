
def get_parameters(run):


  class parameters: 

    def __init__(self, run):
      self.molName     = "N2O"
      self.run         = run
      self.data_format = 'ts'
      self.elEnergy    = 3.7e6
      self.beamline_length = 4
      self.detector_width = 0.04
      self.detector_height = 0.04
      self.probe_FWHM  = 180



      self.QperPix         = 0.028156*(137.0/112.0)
      self.Nlegendres      = 1
      self.NmaxRadBins     = 750
      self.NradAzmBins     = 375 
      self.imgBins         = 895
      self.hasRef          = False
      self.refStagePosCut  = -1

      self.xyz_file            = "N2O.xyz"
      
      self.xyzDir              = "/reg/neh/home/khegazy/analysis/2015/timeBasis_N2O_NO2/simulation/XYZfiles/"
      self.scat_amps_dir       = "/reg/neh/home5/khegazy/simulation/scatteringAmplitudes/3.7MeV/"
      self.simOutputDir        = "/reg/ued/ana/scratch/N2O/simulations/"
      self.preProcOutputDir    = "/reg/ued/ana/scratch/N2O/preProcessing/"
      self.preProcI0OutputDir  = "NULL"
      self.mergeScansOutputDir = "/reg/ued/ana/scratch/N2O/mergeScans/"
      self.scanSearchOutputDir = "/reg/ued/ana/scratch/N2O/scanSearch/"
      self.radialPixelDist     = "/reg/ued/ana/scratch/N2O/radialPixelDist/"
      self.backgroundImage     = "NULL"
      self.backgroundFolder    = "/reg/ued/ana/scratch/N2O/background/"
      self.indexPath           = "/reg/neh/home/khegazy/analysis/radialFitIndices/"
  

      if run == "fullRevival":
        self.backgroundImage = "backgroundImg-20180630_1917.dat"

  return parameters(run)
