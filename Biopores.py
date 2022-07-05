#@ Integer(label="Length filter in voxel", value=60,min=0, max=200, style="slider") t_v
#@ Double(value=0.6, min=0.0, max=1.0, label="Vesselness filter") t_s
#@OUTPUT ImagePlus Roots
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convertService
#@ UIService ui
#@ImageJ ij

# This script segments 3D biopores (based on their tubeness like shape) in 3D CT 
#MLucas 01.2022

import ij.IJ
from net.imagej import Dataset
from ij import ImagePlus
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator
import cmath
from ij import IJ, WindowManager

def Tubeness(data,sigma1, low, high):
	data1 = convertService.convert(data, Dataset)
	slices = data1.getDepth()
	tube = ops.filter().tubeness(data1,sigma1)
	output = ds.create(tube)
	IJ.run("Conversions...", " ");
	tube8bit = ops.convert().uint8(output)
	tube8bit = ds.create(tube8bit)	
	ui.show(tube8bit)
	result2 = IJ.getImage()
	result2.setTitle("result2")
	tube8bit3=IJ.run(result2 ,"Properties...", "channels=1 slices="+str(slices)+" frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");	
	result=IJ.run(tube8bit3, "3D Hysteresis Thresholding", "high="+str(high)+" low="+str(low)+" connectivity");
	result2.close();
	return result
	

def segmentsmallBiopores(scale05):
	results = []
	for i in range(1,5):
		Tubeness(scale05,i,55,65)
		a = IJ.getImage()
		results.append(a)
		
	imp12 = ImageCalculator.run(results[0],results[1], "Max create stack");
	imp13 = ImageCalculator.run(imp12, results[2], "Max create stack");
	imp14 = ImageCalculator.run(imp13, results[3], "Max create stack");
	imp14.show()
	smallBiopores2 = IJ.getImage() 
	smallBiopores2.setTitle("smallBiopores")
	for i in range(0,4):
		results[i].close()
	imp13.close()
	imp12.close()	
	t_s_small=t_s*0.5	
	smallBioporesRemoved= RemoveBloobs(smallBiopores2, t_v, t_s_small, dyn)
	smallBioporesRemoved= IJ.getImage()
	smallBioporesRemoved.setTitle("smallBioporesRemoved");
	i = 1
	while i < 9:
	  IJ.run(smallBioporesRemoved, "Dilate (3D)", "iso=255");
	  i += 1;
	  
	smallBioporesRemoved2 = ImageCalculator.run(smallBioporesRemoved, scale05, "Multiply create stack");
	smallBioporesRemoved2.show();
	smallBioporesRemoved3 = IJ.getImage();
	smallBioporesRemoved4 = smallBioporesRemoved3.resize(total_width,total_height,  total_slices, "none");
	smallBioporesRemoved4.show();
	smallBioporesRemoved4 = IJ.getImage();
	smallBioporesRemoved4.setTitle("smallBioporesRemoved4");	
	return smallBioporesRemoved4

def segmentlargeBiopores(scale02):
	width, height, ch , slices, fr = scale02.getDimensions()
	scale02.setSlice(slices);
	IJ.run(scale02,"Add Slice", "");
	IJ.run(scale02, "Flip Z", "");
	scale02.setSlice(slices+1);
	IJ.run(scale02,"Add Slice", "");
	IJ.run(scale02, "Flip Z", "");
	IJ.run(scale02,"Local Thickness (complete process)", "threshold=255");
	IJ.setThreshold(5, 255);
	IJ.run("Convert to Mask", "method=Default background=Dark black");
	mask = IJ.getImage()
	mask.setTitle("mask");
	results = []
	maxStack = []

	for i in range(2,20):
		Tubeness(mask,i,70,65)
		a = IJ.getImage()
		results.append(a)
	for i in range(20,31):
		Tubeness(mask,i, 60, 55)
		a = IJ.getImage()
		results.append(a)

	a = ImageCalculator.run(results[0], results[1], "Max create stack");	
	maxStack.append(a)
	for i in range(1,28):
		a = ImageCalculator.run(results[i+1],maxStack[i-1], "Max create stack");
		maxStack.append(a)
	maxStack[27].show()
	largeBiopores2 = IJ.getImage() 
	largeBiopores2.setTitle("largeBiopores")
	for i in range(0,29):
		results[i].close()	
	for i in range(0,26):
		maxStack[i].close()
	largeBiopores2.setSlice(1);
	IJ.run(largeBiopores2,"Delete Slice", "");
	largeBiopores2.setSlice(slices+1);
	IJ.run(largeBiopores2,"Delete Slice", "");
	t_s_large=t_s*0.2;
	largeBioporesRemoved= RemoveBloobs(largeBiopores2, t_v, t_s_large, dyn)
	largeBioporesRemoved2= IJ.getImage()
	largeBioporesRemoved2.setTitle("largeBioporesRemoved2");
	i = 1
	while i < 11:
	  IJ.run(largeBioporesRemoved2, "Dilate (3D)", "iso=255");
	  i += 1
	mask.setSlice(1);
	IJ.run(mask,"Delete Slice", "");
	mask.setSlice(slices+1);
	IJ.run(mask,"Delete Slice", "")
	largeBiopores2= ImageCalculator.run(largeBioporesRemoved2, mask, "Multiply create stack");
	largeBiopores2.show();
	largeBiopores3= IJ.getImage();
	largeBioporesfinal = largeBiopores3.resize(total_width,total_height,  total_slices, "none");
	largeBioporesfinal.show();
	largeBioporesfinal1 = IJ.getImage();
	largeBioporesfinal1.setTitle("largeBiioporesfinal1");
	return largeBioporesfinal1

def RemoveBloobs(image, t_v, t_s, dyn):
	image = IJ.getImage()
	image.setTitle("bio");
	smallBioporesRemoved=IJ.run(image,"Distance Transform Watershed 3D", "distances=[Borgefors (3,4,5)] output=[32 bits] normalize dynamic="+str(dyn)+" connectivity=6");
	IJ.run(smallBioporesRemoved,"Analyze Regions 3D", "equivalent_ellipsoid ellipsoid_elongations surface_area_method=[Crofton (13 dirs.)] euler_connectivity=26");	
	IJ.renameResults("biodist-watershed-morpho", "Results");
	rt = ResultsTable.getResultsTable();
	ElliR3 = rt.getColumn(5)
	ElliR2 = rt.getColumn(4)
	ElliR1 = rt.getColumn(3)
	Rb = []
	Ra = []
	Vesselness = []	
	for i in range(len(ElliR3)):
		Rb.append(ElliR3[i]/cmath.sqrt(ElliR2[i]*ElliR1[i])); 
		Ra.append(ElliR2[i]/ElliR1[i]); 
		Vesselness.append(abs(cmath.exp(-Rb[i]*Rb[i])*cmath.exp(-Ra[i]*Ra[i]))); 
		if (Vesselness[i]) >t_v and (ElliR1[i]) > t_s:
			rt.setValue("flag",i ,255) 
		else:
			rt.setValue("flag",i ,0) 	
	rt.updateResults()
	rt.show("flaged")	
	IJ.run("Assign Measure to Label", "results="+"flaged"+" column="+"flag"+" min=0 max=255"); 	
	IJ.run("8-bit", "");
	smallBioporesRemoved = IJ.getImage(smallBioporesRemoved)
	smallBioporesRemoved.setTitle("smallBioporesRemoved")
	return smallBioporesRemoved


thresh = IJ.getImage()
thresh.setTitle("thresh");
total_width, total_height, ch , total_slices, fr = thresh.getDimensions()
scale05 = thresh.resize(total_width/2,total_height/2,  total_slices/2, "none");
dyn =2;
smallBioporesRemoved4 = segmentsmallBiopores(scale05)
scale02 = thresh.resize(total_width/5,total_height/5,  total_slices/5, "none");
largeBioporesfinal2 = segmentlargeBiopores(scale02)
Bioporesfinal = ImageCalculator.run(largeBioporesfinal2,smallBioporesRemoved4, "ADD create stack");
IJ.run("Close All", "")
ui.show(Bioporesfinal)

#print Tube


