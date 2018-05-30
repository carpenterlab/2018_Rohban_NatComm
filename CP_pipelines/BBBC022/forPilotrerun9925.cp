CellProfiler Pipeline: http://www.cellprofiler.org
Version:1
SVNRevision:9925

LoadData:[module_num:1|svn_version:\'Unknown\'|variable_revision_number:4|show_window:False|notes:\x5B\x5D]
    Input data file location:Default Output Folder\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/PILOT_RERUN_GOODZERNICKES
    Name of the file:allplates_PILOT.csv
    Load images based on this data?:Yes
    Base image location:Elsewhere...\x7C/cbnt/cbimage/HCS/GustafsdottirSigrun
    Process just a range of rows?:No
    Rows to process:1,100000
    Group images by metadata?:No
    Select metadata fields for grouping:

LoadSingleImage:[module_num:2|svn_version:\'9899\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Input image file location:Default Input Folder sub-folder\x7C\\g<PlateID>
    Filename of the image to load (Include the extension, e.g., .tif):illumHoechst.mat
    Name the image that will be loaded:illumHoechst
    Filename of the image to load (Include the extension, e.g., .tif):illumER.mat
    Name the image that will be loaded:illumER
    Filename of the image to load (Include the extension, e.g., .tif):illumSyto.mat
    Name the image that will be loaded:illumSyto
    Filename of the image to load (Include the extension, e.g., .tif):illumPh_golgi.mat
    Name the image that will be loaded:illumPh_golgi
    Filename of the image to load (Include the extension, e.g., .tif):illumMito.mat
    Name the image that will be loaded:illumMito

CorrectIlluminationApply:[module_num:3|svn_version:\'9560\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Select the input image:OrigHoechst
    Name the output image:Hoechst
    Select the illumination function:illumHoechst
    Select how the illumination function is applied:Divide
    Select the rescaling method:No rescaling

CorrectIlluminationApply:[module_num:4|svn_version:\'9560\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Select the input image:OrigER
    Name the output image:ER
    Select the illumination function:illumER
    Select how the illumination function is applied:Divide
    Select the rescaling method:No rescaling

CorrectIlluminationApply:[module_num:5|svn_version:\'9560\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Select the input image:OrigSyto
    Name the output image:Syto
    Select the illumination function:illumSyto
    Select how the illumination function is applied:Divide
    Select the rescaling method:No rescaling

CorrectIlluminationApply:[module_num:6|svn_version:\'9560\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Select the input image:OrigMito
    Name the output image:Mito
    Select the illumination function:illumMito
    Select how the illumination function is applied:Divide
    Select the rescaling method:No rescaling

CorrectIlluminationApply:[module_num:7|svn_version:\'9560\'|variable_revision_number:2|show_window:False|notes:\x5B\x5D]
    Select the input image:OrigPh_golgi
    Name the output image:Ph_golgi
    Select the illumination function:illumPh_golgi
    Select how the illumination function is applied:Divide
    Select the rescaling method:No rescaling

IdentifyPrimaryObjects:[module_num:8|svn_version:\'9856\'|variable_revision_number:7|show_window:False|notes:\x5B\x5D]
    Select the input image:Hoechst
    Name the primary objects to be identified:Nuclei
    Typical diameter of objects, in pixel units (Min,Max):15,150
    Discard objects outside the diameter range?:Yes
    Try to merge too small objects with nearby larger objects?:No
    Discard objects touching the border of the image?:Yes
    Select the thresholding method:Otsu Global
    Threshold correction factor:1.0
    Lower and upper bounds on threshold:0.000000,1.000000
    Approximate fraction of image covered by objects?:0.01
    Method to distinguish clumped objects:Shape
    Method to draw dividing lines between clumped objects:Distance
    Size of smoothing filter:20
    Suppress local maxima that are closer than this minimum allowed distance:20
    Speed up by using lower-resolution image to find local maxima?:Yes
    Name the outline image:PrimaryOutlines
    Fill holes in identified objects?:Yes
    Automatically calculate size of smoothing filter?:No
    Automatically calculate minimum allowed distance between local maxima?:No
    Manual threshold:0.0
    Select binary image:None
    Retain outlines of the identified objects?:No
    Automatically calculate the threshold using the Otsu method?:Yes
    Enter Laplacian of Gaussian threshold:0.5
    Two-class or three-class thresholding?:Two classes
    Minimize the weighted variance or the entropy?:Weighted variance
    Assign pixels in the middle intensity class to the foreground or the background?:Foreground
    Automatically calculate the size of objects for the Laplacian of Gaussian filter?:Yes
    Enter LoG filter diameter:5
    Handling of objects if excessive number of objects identified:Erase
    Maximum number of objects:700
    Select the measurement to threshold with:None

IdentifySecondaryObjects:[module_num:9|svn_version:\'9741\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the input objects:Nuclei
    Name the objects to be identified:Cells
    Select the method to identify the secondary objects:Propagation
    Select the input image:Ph_golgi
    Select the thresholding method:Otsu Global
    Threshold correction factor:0.9
    Lower and upper bounds on threshold:0.003,1.0
    Approximate fraction of image covered by objects?:0.01
    Number of pixels by which to expand the primary objects:10
    Regularization factor:0.05
    Name the outline image:CellOutlines
    Manual threshold:0.0
    Select binary image:None
    Retain outlines of the identified secondary objects?:Yes
    Two-class or three-class thresholding?:Three classes
    Minimize the weighted variance or the entropy?:Weighted variance
    Assign pixels in the middle intensity class to the foreground or the background?:Foreground
    Discard secondary objects that touch the edge of the image?:No
    Discard the associated primary objects?:No
    Name the new primary objects:FilteredNuclei
    Retain outlines of the new primary objects?:No
    Name the new primary object outlines:FilteredNucleiOutlines
    Select the measurement to threshold with:None

IdentifyTertiaryObjects:[module_num:10|svn_version:\'9429\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Select the larger identified objects:Cells
    Select the smaller identified objects:Nuclei
    Name the tertiary objects to be identified:Cytoplasm
    Name the outline image:CytoplasmOutlines
    Retain outlines of the tertiary objects?:No

MeasureObjectIntensity:[module_num:11|svn_version:\'9645\'|variable_revision_number:3|show_window:False|notes:\x5B\x5D]
    Hidden:5
    Select an image to measure:Hoechst
    Select an image to measure:ER
    Select an image to measure:Syto
    Select an image to measure:Ph_golgi
    Select an image to measure:Mito
    Select objects to measure:Nuclei
    Select objects to measure:Cells
    Select objects to measure:Cytoplasm

MeasureObjectNeighbors:[module_num:12|svn_version:\'Unknown\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Select objects to measure:Nuclei
    Method to determine neighbors:Within a specified distance
    Neighbor distance:1
    Retain the image of objects colored by numbers of neighbors for use later in the pipeline (for example, in SaveImages)?:No
    Name the output image:ObjectNeighborCount
    Select colormap:Default
    Retain the image of objects colored by percent of touching pixels for use later in the pipeline (for example, in SaveImages)?:No
    Name the output image:PercentTouching
    Select a colormap:Default

MeasureObjectNeighbors:[module_num:13|svn_version:\'Unknown\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Select objects to measure:Cells
    Method to determine neighbors:Within a specified distance
    Neighbor distance:5
    Retain the image of objects colored by numbers of neighbors for use later in the pipeline (for example, in SaveImages)?:No
    Name the output image:ObjectNeighborCount
    Select colormap:Default
    Retain the image of objects colored by percent of touching pixels for use later in the pipeline (for example, in SaveImages)?:No
    Name the output image:PercentTouching
    Select a colormap:Default

MeasureObjectRadialDistribution:[module_num:14|svn_version:\'9610\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Hidden:4
    Hidden:1
    Hidden:1
    Select an image to measure:ER
    Select an image to measure:Ph_golgi
    Select an image to measure:Mito
    Select an image to measure:Syto
    Select objects to meaasure:Cells
    Object to use as center?:These objects
    Select objects to use as centers:None
    Number of bins:4

MeasureObjectSizeShape:[module_num:15|svn_version:\'1\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Select objects to measure:Nuclei
    Select objects to measure:Cells
    Select objects to measure:Cytoplasm
    Calculate the Zernike features?:Yes

MeasureTexture:[module_num:16|svn_version:\'1\'|variable_revision_number:1|show_window:False|notes:\x5B\x5D]
    Hidden:5
    Hidden:3
    Hidden:2
    Select an image to measure:Hoechst
    Select an image to measure:ER
    Select an image to measure:Syto
    Select an image to measure:Ph_golgi
    Select an image to measure:Mito
    Select objects to measure:Cells
    Select objects to measure:Cytoplasm
    Select objects to measure:Nuclei
    Texture scale to measure:3
    Texture scale to measure:5
    Number of angles to compute for Gabor:4

SaveImages:[module_num:17|svn_version:\'9917\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:CellOutlines
    Select the module display window to save:None
    Select method for constructing file names:Name with metadata
    Select image name for file prefix:None
    Enter file name with metadata:Outline_\\g<CPD_WELL_POSITION>_\\g<Site>
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:png
    Output file location:Default Input Folder sub-folder\x7C\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:Yes
    Select how often to save:Every cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:Yes
    Create subfolders in the output folder?:No

ExportToDatabase:[module_num:18|svn_version:\'9921\'|variable_revision_number:16|show_window:True|notes:\x5B\x5D]
    Database type:MySQL
    Database name:2008_12_04_Imaging_CDRP_for_MLPCN
    Add a prefix to table names?:Yes
    Table prefix:Pilot_Rerun_illumcorr
    SQL file prefix:SQL_
    Output file location:Default Output Folder\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/PILOT_RERUN_GOODZERNICKES
    Create a CellProfiler Analyst properties file?:Yes
    Database host:imgdb01
    Username:cpadmin
    Password:cPus3r
    Name the SQLite database file:DefaultDB.db
    Calculate the per-image mean values of object measurements?:No
    Calculate the per-image median values of object measurements?:No
    Calculate the per-image standard deviation values of object measurements?:No
    Calculate the per-well mean values of object measurements?:No
    Calculate the per-well median values of object measurements?:No
    Calculate the per-well standard deviation values of object measurements?:No
    Export measurements for all objects to the database?:All
    Select the objects:
    Maximum # of characters in a column name:64
    Create one table per object or a single object table?:One table per object type
    Enter an image url prepend if you plan to access your files via http (leave blank if local):
    Write image thumbnails directly to the database?:No
    Select the images you want to save thumbnails of:

CreateBatchFiles:[module_num:19|svn_version:\'9429\'|variable_revision_number:4|show_window:True|notes:\x5B\x5D]
    Store batch files in default output folder?:Yes
    Output folder path:/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/PILOT_RERUN_GOODZERNICKES
    Are the cluster computers running Windows?:No
    Hidden\x3A in batch mode:No
    Hidden\x3A default input folder at time of save:/cbnt/cbimage/HCS/GustafsdottirSigrun
    Hidden\x3A SVN revision number:9655
    Local root path:/cbnt/cbimage/HCS/GustafsdottirSigrun
    Cluster root path:/cbnt/cbimage/HCS/GustafsdottirSigrun
