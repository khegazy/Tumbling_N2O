This code creates diffraction patterns from saved angular distributions located:
      /reg/neh/home/khegazy/analysis/tumblingN2O/simulation/rotation/output/alignmentPDFs

To run code 
      ./submitClusterJobs.sh folderName Njobs
      folderName = diffraction 
      Njobs = 43 (Number of angular PDFS with name job_#...)

submitClusterJobs.sh calls runDiffraction.sh which then runs n2oDiffSim.exe. 
jobNum is fed into runDiffraction.sh
All other parameters (Nmolecules, pdfName, ...) are passed into n2oDiffSim.exe from runDiffraction.sh


Combining Results
  After the jobs finish they must be merged
  cd output/diffraction
  hadd fullSimulation.root job*
