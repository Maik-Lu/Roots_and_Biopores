directory="PATH";
ThresholdMethod="Otsu" //"Own value" ;
length_filter =100;
vesselness_filter =0.78;
mask= 0.1;
radius= 1;
threshold_r=3;
/// This bath enables to segment roots in all images (.mhd files) of a folder. It save a binary .mhd file containing the resulting root as foreground (255) 

filelist = getFileList(directory) ;
filelist = Array.sort(filelist);
	for (g = 0; g < lengthOf(filelist); g++) { // each file of a folder
	
		image_title = substring(filelist[g], 0,lengthOf(filelist[g])-1); //cut image title as needed 

		file= directory+ filelist[g]+image_title;

	    	open(directory+ filelist[g]+"nlm_ring_"+image_title+".mhd");	
	    	param = ReadParameters(file+"_2021_Thresh_table.txt");
			print(file);

			run("Roots ", "t_s="+length_filter+" t_v="+vesselness_filter+" mask="+mask+" radius="+radius+" threshold_value="+ThresholdMethod+" upper_threshold="+upper_threshold+" threshold_ratio="+threshold_r+"");
			run("Median 3D...", "x="+1+" y="+1+" z="+1+""); 
  			rename("root_system");
    
 		    run("MHD/MHA ...", "save="+directory+ filelist[g]+ image_title+"_root.mhd");
	    	run("Clear Results");
	    	
		    close("results");
			close("*");
			run("Collect Garbage"); 
	
}	
