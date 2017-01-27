#!/bin/sh

echo "Default is wget, if wget is not present, use curl argument. ./update curl"

CMD="curl -kC - -O"  ;
CMD_VER="curl -k" ;
echo $1
if [ $1 == "wget" ] ;
  then 
   echo "using wget" ;
   CMD="wget"
   CMD_VER="wget -qO-"
fi

echo "Backup old file to tar"
tar -cf backup_`date +%F_%s`.tar lib/ bin/xrf_maps VERSION.TXT 

if [ -d "backups" ]; then
  echo "backups dir already exists"
else 
  echo "mkdir backups"
  mkdir backups
fi

mv *.tar backups/

echo "Removing old files"
rm VERSION.TXT
rm -rf bin
rm -rf lib

if [ $1 == "clean" ] ;
   then
   exit;
fi

echo "Getting new files"
$CMD https://jenkins.aps.anl.gov/view/TXM%20Linux64/job/libSZIP-2.1/ws/lib/libsz.so.2
$CMD https://jenkins.aps.anl.gov/view/TXM%20Linux64/job/libHDF5/ws/libHDF5/lib/libhdf5.so.8
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/netcdf_4_linux64/ws/build_lib/lib/libnetcdf.so.11
$CMD https://jenkins.aps.anl.gov/job/gcc_6_2_linux64/ws/build_dir/lib64/libstdc++.so.6
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/ws/reference/henke.xdr
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/ws/reference/xrf_library.csv
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/ws/bin/libxrf_fit.so
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/ws/bin/libxrf_io.so
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/ws/bin/xrf_maps

#get osmesa
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/OSMesa/ws/build/lib/libOSMesa.so.8
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/OSMesa/ws/build/lib/libglapi.so.0
#get vtk libs
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkDomainsChemistryOpenGL2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersFlowPaths-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersGeneric-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersHyperTree-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersParallelImaging-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersPoints-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersProgrammable-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersSMP-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersSelection-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersTexture-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersVerdict-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkverdict-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkGeovisCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkproj4-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOAMR-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOEnSight-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOExodus-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOExport-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingGL2PSOpenGL2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkgl2ps-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOImport-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOInfovis-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtklibxml2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOLSDyna-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOMINC-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOMovie-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkoggtheora-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOPLY-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOParallel-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkjsoncpp-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOParallelXML-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOSQL-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtksqlite-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOTecplotTable-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOVideo-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingMorphological-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingStatistics-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingStencil-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkInteractionImage-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingContextOpenGL2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingImage-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingLOD-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingVolumeOpenGL2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkViewsContext2D-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkViewsInfovis-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkDomainsChemistry-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersAMR-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersParallel-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkexoIIc-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOGeometry-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIONetCDF-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkNetCDF_cxx-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkNetCDF-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkhdf5_hl-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkhdf5-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkParallelCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOLegacy-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingOpenGL2-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkglew-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingMath-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkChartsCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingContext2D-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersImaging-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkInfovisLayout-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkInfovisCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkViewsCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkInteractionWidgets-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersHybrid-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingGeneral-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingSources-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersModeling-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingHybrid-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOImage-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkDICOMParser-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkmetaio-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkpng-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtktiff-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkjpeg-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkInteractionStyle-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersExtraction-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersStatistics-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingFourier-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkalglib-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingAnnotation-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingColor-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingVolume-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkImagingCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOXML-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOXMLParser-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkIOCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkexpat-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingLabel-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingFreeType-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkRenderingCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonColor-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersGeometry-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersSources-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersGeneral-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonComputationalGeometry-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkFiltersCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonExecutionModel-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonDataModel-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonTransforms-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonMisc-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonMath-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonSystem-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkCommonCore-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtksys-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkfreetype-7.1.so.1
$CMD https://jenkins.aps.anl.gov/view/XRF_Maps/job/VTK_7_Linux64/ws/VTK-Release-build/lib/libvtkzlib-7.1.so.1

mkdir lib
chmod +x *.so*
mv *.so* lib

if [ -d "reference" ]; then
  echo "reference dir already exists"
else 
  echo "mkdir reference"
  mkdir reference
fi

mv henke.xdr reference/
mv xrf_library.csv reference/

echo "Generating version info"
echo "XRF-Maps" > VERSION.TXT
$CMD_VER https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/lastBuild/  > VERSION.TXT

chmod +x xrf_maps

mkdir bin
mv xrf_maps bin/

echo "Creating shortcut"
echo "#!/bin/sh" > xrf_maps.sh
echo " ">> xrf_maps.sh
echo "export LD_LIBRARY_PATH=./lib:../lib">> xrf_maps.sh
echo "cd bin">> xrf_maps.sh
echo "./xrf_maps $@">> xrf_maps.sh
echo "cd ..">> xrf_maps.sh

echo "Finished updating."
