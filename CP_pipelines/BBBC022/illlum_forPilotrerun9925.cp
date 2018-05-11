CellProfiler Pipeline: http://www.cellprofiler.org
Version:1
SVNRevision:9925

LoadData:[module_num:1|svn_version:\'Unknown\'|variable_revision_number:4|show_window:True|notes:\x5B\x5D]
    Input data file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr
    Name of the file:allplates_PILOT.csv
    Load images based on this data?:Yes
    Base image location:Elsewhere...\x7C/cbnt/cbimage/HCS
    Process just a range of rows?:No
    Rows to process:1,100000
    Group images by metadata?:Yes
    Select metadata fields for grouping:PlateID

CorrectIlluminationCalculate:[module_num:2|svn_version:\'Unknown\'|variable_revision_number:1|show_window:True|notes:\x5B\x5D]
    Select the input image:OrigHoechst
    Name the output image:IllumHoechst
    Select how the illumination function is calculated:Regular
    Dilate objects in the final averaged image?:No
    Dilation radius:1
    Block size:60
    Rescale the illumination function?:Yes
    Calculate function for each image individually, or based on all images?:All
    Smoothing method:Median Filter
    Method to calculate smoothing filter size:Manually
    Approximate object size:10
    Smoothing filter size:100
    Retain the averaged image for use later in the pipeline (for example, in SaveImages)?:No
    Name the averaged image:IllumBlueAvg
    Retain the dilated image for use later in the pipeline (for example, in SaveImages)?:No
    Name the dilated image:IllumBlueDilated

CorrectIlluminationCalculate:[module_num:3|svn_version:\'Unknown\'|variable_revision_number:1|show_window:True|notes:\x5B\x5D]
    Select the input image:OrigER
    Name the output image:IllumER
    Select how the illumination function is calculated:Regular
    Dilate objects in the final averaged image?:No
    Dilation radius:1
    Block size:60
    Rescale the illumination function?:Yes
    Calculate function for each image individually, or based on all images?:All
    Smoothing method:Median Filter
    Method to calculate smoothing filter size:Manually
    Approximate object size:10
    Smoothing filter size:100
    Retain the averaged image for use later in the pipeline (for example, in SaveImages)?:No
    Name the averaged image:IllumBlueAvg
    Retain the dilated image for use later in the pipeline (for example, in SaveImages)?:No
    Name the dilated image:IllumBlueDilated

CorrectIlluminationCalculate:[module_num:4|svn_version:\'Unknown\'|variable_revision_number:1|show_window:True|notes:\x5B\x5D]
    Select the input image:OrigSyto
    Name the output image:IllumSyto
    Select how the illumination function is calculated:Regular
    Dilate objects in the final averaged image?:No
    Dilation radius:1
    Block size:60
    Rescale the illumination function?:Yes
    Calculate function for each image individually, or based on all images?:All
    Smoothing method:Median Filter
    Method to calculate smoothing filter size:Manually
    Approximate object size:10
    Smoothing filter size:100
    Retain the averaged image for use later in the pipeline (for example, in SaveImages)?:No
    Name the averaged image:IllumBlueAvg
    Retain the dilated image for use later in the pipeline (for example, in SaveImages)?:No
    Name the dilated image:IllumBlueDilated

CorrectIlluminationCalculate:[module_num:5|svn_version:\'Unknown\'|variable_revision_number:1|show_window:True|notes:\x5B\x5D]
    Select the input image:OrigPh_golgi
    Name the output image:IllumPh_golgi
    Select how the illumination function is calculated:Regular
    Dilate objects in the final averaged image?:No
    Dilation radius:1
    Block size:60
    Rescale the illumination function?:Yes
    Calculate function for each image individually, or based on all images?:All
    Smoothing method:Median Filter
    Method to calculate smoothing filter size:Manually
    Approximate object size:10
    Smoothing filter size:100
    Retain the averaged image for use later in the pipeline (for example, in SaveImages)?:No
    Name the averaged image:IllumBlueAvg
    Retain the dilated image for use later in the pipeline (for example, in SaveImages)?:No
    Name the dilated image:IllumBlueDilated

CorrectIlluminationCalculate:[module_num:6|svn_version:\'Unknown\'|variable_revision_number:1|show_window:True|notes:\x5B\x5D]
    Select the input image:OrigMito
    Name the output image:IllumMito
    Select how the illumination function is calculated:Regular
    Dilate objects in the final averaged image?:No
    Dilation radius:1
    Block size:60
    Rescale the illumination function?:Yes
    Calculate function for each image individually, or based on all images?:All
    Smoothing method:Median Filter
    Method to calculate smoothing filter size:Manually
    Approximate object size:10
    Smoothing filter size:100
    Retain the averaged image for use later in the pipeline (for example, in SaveImages)?:No
    Name the averaged image:IllumBlueAvg
    Retain the dilated image for use later in the pipeline (for example, in SaveImages)?:No
    Name the dilated image:IllumBlueDilated

SaveImages:[module_num:7|svn_version:\'Unknown\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:IllumHoechst
    Select the module display window to save:None
    Select method for constructing file names:Single name
    Select image name for file prefix:None
    Enter single file name:illumHoechst
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:mat
    Output file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr/\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:No
    Select how often to save:Last cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:No
    Create subfolders in the output folder?:No

SaveImages:[module_num:8|svn_version:\'Unknown\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:IllumER
    Select the module display window to save:None
    Select method for constructing file names:Single name
    Select image name for file prefix:None
    Enter single file name:illumER
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:mat
    Output file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr/\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:No
    Select how often to save:Last cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:No
    Create subfolders in the output folder?:No

SaveImages:[module_num:9|svn_version:\'Unknown\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:IllumSyto
    Select the module display window to save:None
    Select method for constructing file names:Single name
    Select image name for file prefix:None
    Enter single file name:illumSyto
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:mat
    Output file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr/\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:No
    Select how often to save:Last cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:No
    Create subfolders in the output folder?:No

SaveImages:[module_num:10|svn_version:\'Unknown\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:IllumPh_golgi
    Select the module display window to save:None
    Select method for constructing file names:Single name
    Select image name for file prefix:None
    Enter single file name:illumPh_golgi
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:mat
    Output file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr/\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:No
    Select how often to save:Last cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:No
    Create subfolders in the output folder?:No

SaveImages:[module_num:11|svn_version:\'Unknown\'|variable_revision_number:5|show_window:True|notes:\x5B\x5D]
    Select the type of image to save:Image
    Select the image to save:IllumMito
    Select the module display window to save:None
    Select method for constructing file names:Single name
    Select image name for file prefix:None
    Enter single file name:illumMito
    Do you want to add a suffix to the image file name?:No
    Text to append to the image name:Do not use
    Select file format to use:mat
    Output file location:Elsewhere...\x7C/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr/\\g<PlateID>
    Image bit depth:8
    Overwrite existing files without warning?:No
    Select how often to save:Last cycle
    Rescale the images? :No
    Select colormap:gray
    Update file names within CellProfiler?:No
    Create subfolders in the output folder?:No

CreateBatchFiles:[module_num:12|svn_version:\'Unknown\'|variable_revision_number:4|show_window:True|notes:\x5B\x5D]
    Store batch files in default output folder?:Yes
    Output folder path:/imaging/analysis/2008_12_04_Imaging_CDRP_for_MLPCN/Pilot_rerun_illumcorr
    Are the cluster computers running Windows?:No
    Hidden\x3A in batch mode:No
    Hidden\x3A default input folder at time of save:/cbnt/cbimage/HCS
    Hidden\x3A SVN revision number:9881
    Local root path:/cbnt/cbimage/HCS
    Cluster root path:/cbnt/cbimage/HCS
