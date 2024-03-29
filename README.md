# MATLAB/Octave Cast-off Reconstruction

  **Reconstructs stains from cast-off event to reproduce the motion of cast-off.**
  
  >>>
    Scott McCleary
    Email: scott.thomas.mccleary@gmail.com | daniel.attinger@gmail.com
    Phone: (515) 975-5544
    Spatter Stains to Cast-off Reconstruction
    Center for Statistics and Applications in Forensic Evidence
    Department of Mechanical Engineering - Attinger Lab
    Iowa State University
  >>>
  
  _**Last Updated 01/06/2024**_ &#127881;

  ## Note to user:
  
  ### Surface and Stain Orientation Definitions:

  ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/e71faca0-0dfe-495f-82af-293925295205)

  ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/16bed7fb-f392-4541-965e-eb3aa8c83663)

  - _Surface Tangential Vector (t&#770;)_ :
    - _Definition_ : vector tangential to a surface aligned in the direction of gravity (_**-**z-direction_) for all _non-horizonal surfaces_
      - for perfectly horizontal surfaces, e.g. _upward_ surface (ceiling) and _downward_ surface (floor), the -z-direction aligns with the Surface Normal Vector. In these instances, the _Surface Tangential Vector (t&#770;)_ is defined in the _**+**x-direction_
  - _Surface Normal Vector (n&#770;)_ :
    - _Definition_ : vector normal to a surface pointing _**outside**_ the room
  - _**Gamma (&gamma;) :  stain directional angle (relative to impacted surface):**_
    - **_Definition_ : stain glancing angle between the projected trajectory vector onto the impacted surface and the tangent vector to the surface**
  - _Alpha (&alpha;) : stain impact angle_
    - _Definition_ : stain impact angle between the trajectory vector and the impacted surface 

  ## Required Repository Files to run the code:
  
  - Spatter Measurement Data, e.g. `Ink_Trial_INPUT.csv`, `Swineblood_Trial_INPUT.csv`
  - `Castoff_Reconstruction_DRIVER.m`
  - `DRIVER.csv` produces `DRIVER.mat` required for `Castoff_Reconstruction_MAIN.m`
  - `Castoff_Reconstruction_MAIN.m`
  - `Castoff_Reconstruction_FUNC.m`
  - `lineSegmentIntersect.m`
  - `meshVolume.m`
  - `point_to_line.m`
  - `gauss_distribution.m`
  - `CircleFitByPratt.m`
  - `Castoff_Reconstruction_POST.m`
  - `combvec2.m`
  - `inpolyhedron.m` 
  - `triangulateFaces.m` 
  - `linecirc.m`
  - `generate_input.m`
    - `linecirc.m`
    - `plane_line_intersect.m`
    - `rotate_3D.m`
    - `triangulateFaces.m`
    
    ### Licenses:
    
    **All licenses for third party scripts are included and must be kept with provided scripts. If third party materials were not cited within the repository Licenses folder, this was not intentional by the author.**

  ## Required Installation for GNU Octave Compatibility: (Suggested Installation Order)
  
  - Updated [GNU Octave](https://www.gnu.org/software/octave/download.html) (written and tested with 5.2.0)
    - Here is the link for previous [GNU Octave version downloads (5.2.0) for windows](https://mirrors.sarata.com/gnu/octave/windows/)
    - All packages can be updated to the latest version by running:
      > pkg update
  - Updated [Java](https://www.java.com/en/download/win10.jsp)
  - Psychtoolbox-3 Requirements [install prior to Psychtoolbox download](http://psychtoolbox.org/download.html#download-problems)
    - Psychtoolbox is a free set of Matlab and GNU Octave functions for vision and neuroscience research.
    - Subversion 1.7.x command-line client (e.g. Sliksvn) MUST INSTALL COMPLETE SOFTWARE(https://sliksvn.com/download/)
      - SlikSVN is a standalone command-line Subversion client for Windows.
    - 64-Bit GStreamer-1.16.0 MSVC runtime or later versions MUST INSTALL COMPLETE SOFTWARE(https://gstreamer.freedesktop.org/data/pkg/windows/1.16.0/gstreamer-1.0-msvc-x86_64-1.16.0.msi)		 - GStreamer is a library for constructing graphs of media-handling components.
  - Download [Psychtoolbox Toolbox Version 3 (PTB-3)](http://psychtoolbox.org/)
    - Unzip download in new folder C:\toolbox
    - Open Octave, copy and paste the following script to the command window to install Psychtoolbox:
      > cd C:\toolbox
      > DownloadPsychtoolbox _C:\toolbox_
  - Follow command window prompts to finish installation.
  - Prior to running Castoff_Recosntruction_DRIVER.m or Castoff_Reconstruction_MAIN.m the io and signal packages need to be loaded.
    - Load packages on an as needed basis by copying and pasting the following script to the command window:
      > pkg load io
      > pkg load signal
    - Required [packages for download](https://octave.sourceforge.io/packages.php).
    - Automatically load packages at Octave startup with the following:
      - Go to _C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup_
      - Open octaverc with a text editor.
      - Copy and paste the following lines to the end of the document:
      > pkg load io
      > pkg load signal
  - Save file.
    - If octaverc does not exist, copy, paste, and save the following to a text file and save to the directory (C:\Octave\Octave-5.2.0\mingw64\share\octave\site\m\startup):
      > pkg load io
      > pkg load signal

  ## Instructions to Run a first spatter pattern example (`Ink_Trial_INPUT.csv`), with stains measured by hand or with Hemospat:
  
  1. Make sure you have Matlab or Octave installed as per the installation instructions above
  2. Save all files included in the distribution to the same directory.
  3. `Run` (F5) the DRIVER (`Castoff_Reconstruction_DRIVER.m`).
  4. Check that the data name line 64 in the MAIN (`Castoff_Reconstruction_MAIN.m`) reads `INK_Trial_INPUT_DRIVER.mat`
  5. `Run` (F5) the MAIN (`Castoff_Reconstruction_MAIN.m`); producing figure(3) and figure(4+) which simulate Figures 5 and 7 of _**McCleary et al FSI 2021**_, respectively.
  6. Maximize figure 4 to observe the reconstructed swing regions in blue, green and red
  7. MAIN also outputs the `Total Elapsed Cluster Analysis Time:`, in seconds (s), with the total program run time excluding any variable time from user input. See User Inputs: (DRIVER) res for estimated runtimes.
  8. MAIN also save the Resultant Variables as a .mat file denoted by Castoff_Reconstruction.mat
  9. MAIN and Output the Cast-off Reconstruction Results with warnings as a txt-file denoted by `Castoff_Reconstruction_OUTPUT.txt` displayed in Figures(4+) (figure number is dependent on clustering method). Figure(1) shows the Pratt fit automated reference point. Figure(3) outputs each clustered cast-off reconstructed arc.
  10. Run again after changing the Stain Width Measurement Uncertainty (mm) or Cast-off Reconstruction Resolution (cm) in the input file `Ink_Trial_INPUT.csv`
  
  ## Instructions to Run a second spatter pattern example (`FARO_Trial_10_INPUT.csv`), with stains measured with FARO:
  
  1. Make sure you have Matlab or Octave installed as per the installation instructions above
  2. Save all files included in the distribution to the same directory.
  3. Modify the data name line 64 in the MAIN (`Castoff_Reconstruction_MAIN.m`) to `FARO_Trial_10_INPUT_DRIVER.mat`
  4. `Run` (F5) the DRIVER (`Castoff_Reconstruction_DRIVER_FARO.m`).
  5. `Run` (F5) the MAIN (`Castoff_Reconstruction_MAIN.m`); producing figure(3) and figure(4+) which simulate Figures 5 and 6 of _**McCleary et al FSI 2021**_, respectively.
  6. Maximize figure 4 to observe the reconstructed swing regions in blue, green and red
  7. MAIN also outputs the `Total Elapsed Cluster Analysis Time:`, in seconds (s), with the total program run time excluding any variable time from user input. See User Inputs: (DRIVER) res for estimated runtimes.
  8. MAIN also save the Resultant Variables as a .mat file denoted by Castoff_Reconstruction.mat
  9. MAIN and Output the Cast-off Reconstruction Results with warnings as a txt-file denoted by `Castoff_Reconstruction_OUTPUT.txt` displayed in Figures(4+) (figure number is dependent on clustering method). Figure(1) shows the Pratt fit automated reference point. Figure(3) outputs each clustered cast-off reconstructed arc.
  10. Run again after changing the Stain Width Measurement Uncertainty (mm) or Cast-off Reconstruction Resolution (cm) in the input file `Ink_Trial_INPUT.csv`
  
  ## Instructions to Run your own data:
    
  1. Make sure you have Matlab or Octave installed as per the installation instructions above
  2. Enter spatter pattern information in an input `_INPUT.csv` file with same format as either of the examples above
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/8439dd09-a8ae-4bcf-9ef4-128df1dbb1f1)

  3. Modify the data name line 64 in the MAIN (`Castoff_Reconstruction_MAIN.m`) to `input_DRIVER.mat`
  4. Follow same instructions (4-10) as either examples above.

  ## Expected Outputs : `Castoff_Reconstruction_DRIVER.m`

  - creates a `_DRIVER.mat` file as input for `Castoff_Reconstruction_MAIN.m`
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/6ff713dc-376a-4d9d-90e9-20ab2478ed36)

  ## Expected Outputs : `Castoff_Reconstruction_MAIN.m`

  ### Figures:

  - `figure(2)`
    - **inputs overview** : an initial look at the inputted room dimensions, stain impact locations, stain trajectories, automated reference point, and best fit plane
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/1fc9917a-a914-49d9-b143-77a18550afd5)
   
  - `figure(3)`
    - **reconstruction overview** : an active display of clustered stains & trajectories and clustered planes with projected stains reference during spatter reconstruction
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/030f229c-67de-47d1-9b58-ae359c9b99fe)
   
  - `figure(4)`
    - **reconstruction results** : inputs overview overlaid with user-defined three-dimensional region of statistical likelihood
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/0778fa92-9148-4a43-af0b-72ebfc510716)
   
  ### Files:

  - saved simulation workspace and textual overview of `INPUTS` (resolution, user preferences, inputed stain details) and `Cast-off Reconstruction Results:` (High, Medium, and Low Percentile reconstructed XYZ-Vertices & Faces)
    - ![image](https://github.com/scottres/CastoffReconstruction/assets/66079745/90ca379c-8240-40f6-b1bd-79e8417dbf5c)

  ### Note: If actual cast-off motion is known (for research purposes only):
  
  - Enter known cast-off motion path x, y, and z-coordinates into `_MOTION.csv` with same format as `Ink_Trial_MOTION.csv`
  - Follow same instructions (4-10) as either examples above.
  - `Run` (F5) the Post-processor (`Castoff_Reconstruction_POST.m`)

  ### Figure Displaying:
  
  - Show and hide legends from resultant figures using:
    > legend show
    > legend hide
  - Rotate, zoom in/out, and pan figure controls are located on the top of the figures to allow for other viewing angles. Also, for all figures besides Figure(1), view(az,el) can be used for precise three-dimensional viewing angles where az is the azimuth angle (in degrees) and el is the polar (or elevation) angle (in degrees).
