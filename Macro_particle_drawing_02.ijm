//************* Drawing spheres in 3D image **********************
//
// Written by Hiroshi Koyama in National Institute for Basic Biology, Japan, 2022/02/03
// Revised by Hiroshi Koyama, 2022/07/21
// hkoyama@nibb.ac.jp

// <Description>
// Read a text file containing xyz-coordinates of multiple particles,
// where the particles are derived from two sample images (e.g. one before fixation, the other after fixation).
// Here, the two sample images called image-before and image-after.
// Each particle is drawn as a sphere with a given radius in two 3D images for the image-before and -after, respectively.

// <Tips for visualization after obtaining the output images>
// If you set pixel sizes of the output images to be the same as the original before-image, you can easily merge the particle image and the before-image in one image.
// Before the merge, you may use ImageJ>Image>Scale in the case that the numbers of y or z pixels/voxels are not equal between the above two images.
// Then, apply ImageJ>Image>Color>Merge Channels...
// For visualization, maximum intensity projection (ImageJ>Image>Stacks>Z Projects...) and 3D viewer (Fiji>Plugins>3D Viewer) may be powerful.
// Good choice of lookup tables improves visualization: ImageJ>Image>Lookup Tables>Green, Magenta, 3-3-2-RGB, Fire, Spectrum, etc.
 
// <Requirement>
// Fiji, version 1.53 or maybe later.
// The name of the text file is output_registered_particles_all.txt with a given format.
// Insted of the above txt file, you can read output_registered_particels_landmarks.txt,
// if you choose that file in "Basic parameters to be defined by users" below.

// <Reference>
// Koyama et al. "An ImageJ-based tool for three-dimensional registration between different types of microscopic images"
// Cite the reference above when you use this macro for any publication.
//
//******************************************************************


//****** Basic parameters to be defined by users ******

print("Please set several values in 'Parameters to be defined by users' on the macro code.\n");

//Choose input text file
input_name = "output_registered_particles_all.txt";
//input_name = "output_registered_particles_landmarks.txt";

image_size_x = 262;		//The pixel number of output images along x axis, which is the same as the original before-image.
image_size_y = 262;		//The pixel number of output images along y axis, which is the same as the original before-image

//Unit of input text file
//If the unit is pixel, um_pixel should be set as 1.0.
//If the unit is um or something, um_pixel should be set so that the values in the input text file will be converted as pixel unit.
///The relationship between um and pixel in the orignal before-image is obtained from ImageJ>Image>Properties...
um_per_pixel = 0.4143204;

//These values should be the same as those used in Macro_3D_particle_registration_220721_06_v2.ijm
x_scale_before = 1.0;
y_scale_before = 1.0;

radius_of_particle=5.0;				//Particle radius to be drawn

modeOFintensity=0;	//0: Intensities between paired particles are the same corresponding to
					//the ID in the before-image.
					//Values other than 0, such as 1: Intensities between paired particles are different
					//and are set according to the original ID.
					//Final intensities are the summation of the above value based on the ID
					//and "foot" defined below as a global variable.
view3D_selection=0; //Select one output image to be 3D-viewed by plugins>3D Viewer.
					//	0 for merge image, 1 for before-image, 2 for after-image
					//However, the 3D Viewer may not work well in Fiji-macro, unfortunately.					
					//Alternatively, you can manually do 3D-viewing for each output three images.
//******************************************************

//*** Global variables ***
var xx1, yy1, zz1, xx2, yy2, zz2, id1, id2;
var p_num, max_id;
var radius, modeOFintensity;
var ix, iy, iz;
var um_pixel;
var foot=100;
var head=13;	//13 means exclusion of header texts in the input text file.

//*** Read input text file ***
dir_input = getDirectory("Choose a Directory ");
input_file = dir_input+input_name;
str1 = File.openAsString(input_file);
//print(str1);
str2 = split(str1, "\t\n\r,");
num = str2.length;
p_num = (num-head)/9;
print("input file name = ", input_file);
print("Particle number = ", p_num);
print("Other values in the input file: ", num, str2[15], str2[num-1]);
print("Basic parameters: ", image_size_x, image_size_y, um_per_pixel, x_scale_before, y_scale_before, radius_of_particle, modeOFintensity, foot);

xx1=newArray(p_num);
yy1=newArray(p_num);
zz1=newArray(p_num);
xx2=newArray(p_num);
yy2=newArray(p_num);
zz2=newArray(p_num);
id1=newArray(p_num);
id2=newArray(p_num);

um_pixel = um_per_pixel;
radius = radius_of_particle;

//*** Interpretation of xyz-coordinates of each particle ***
xyz_coordinate_read();

//*** Scale revision of xyz-coordinates ***
ix = image_size_x;
iy = image_size_y*y_scale_before/x_scale_before;
scale_xyz_revision ();

//*** Sphere drawing for each particles ***
sphere_drawing (modeOFintensity);

//*** save output data and 3D viewing ***
save_3Dview (view3D_selection);

//***************** end of macro ******************************


//***************** functions *********************************

function save_3Dview (view3D_image_selection)
{
	//save output image
	selectWindow("image_before");
	saveAs("Tiff", dir_input+"/image_before.tif");
	selectWindow("image_after");
	saveAs("Tiff", dir_input+"/image_after.tif");

	//merge before and after images and save
	run("Merge Channels...", "c1=image_before.tif c2=image_after.tif create keep");
	saveAs("Tiff", dir_input+"/image_merge_before_after.tif");

	//*** 3D viewer by plugins>3D viewer *** 	
	//3D viewer for merge image
	if (view3D_image_selection==0){
		selectWindow("image_merge_before_after.tif");
		run("3D Viewer");
		call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
		call("ij3d.ImageJ3DViewer.add", "image_merge_before_after.tif", "None", "image_merge_before_after.tif", "0", "true", "true", "true", "2", "0");
	}
	//3D viewer for before image
	selectWindow("image_before.tif");
	run("3-3-2 RGB");	//LUT
	if (view3D_image_selection==1){
		run("3D Viewer");
		call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
		call("ij3d.ImageJ3DViewer.add", "image_before.tif", "None", "image_before.tif", "0", "true", "true", "true", "2", "0");
	}
	//3D viewer for after image
	selectWindow("image_after.tif");
	run("3-3-2 RGB");	//LUT
	if (view3D_image_selection==2){
		run("3D Viewer");
		call("ij3d.ImageJ3DViewer.setCoordinateSystem", "false");
		call("ij3d.ImageJ3DViewer.add", "image_after.tif", "None", "image_after.tif", "0", "true", "true", "true", "2", "0");
	}

	//*** save log file ***
	selectWindow("Log");
	saveAs("Text", dir_input+"/Log_particle_drawing.txt");		
}

function sphere_drawing (modeOFintensity2)
{
	max_intensity = max_id + foot;
	if (max_intensity <256){
		newImage("image_before", "8-bit black", ix, iy, iz);
		newImage("image_after", "8-bit black", ix, iy, iz);
	}
	else{
		newImage("image_before", "16-bit black", ix, iy, iz);
		newImage("image_after", "16-bit black", ix, iy, iz);
	}

	for (ii=0; ii<p_num; ii++){
		sphere_draw_each (ii, modeOFintensity2);	
	}
}

function sphere_draw_each(i, modeOFintensity3)
{
	//drawing of before_image
	selectWindow("image_before");
	z1 = parseInt(zz1[i]-radius);
	if (z1<1){	z1 = 1;	}
	z2 = parseInt(zz1[i]+radius);
	if (z2>iz){	z2 = iz;	}
	
	y1 = parseInt(yy1[i]-radius);
	if (y1<0){	y1 = 0;	}
	y2 = parseInt(yy1[i]+radius);
	if (y2>=iy){	y2 = iy-1;	}	
	
	x1 = parseInt(xx1[i]-radius);
	if (x1<0){	x1 = 0;	}
	x2 = parseInt(xx1[i]+radius);
	if (x2>=ix){	x2 = ix-1;	}	

	for (z=z1; z<=z2; z++){
		setSlice(z);
		for (y=y1; y<y2; y++){
			for (x=x1; x<x2; x++){
				dist2 = pow(xx1[i]-x, 2) + pow(yy1[i]-y, 2) + pow(zz1[i]-z, 2);
				dist = sqrt(dist2);
				if (dist<=radius){
					setPixel(x, y, (id1[i]+foot));
				}
			}
		}
	}
	wait(100);

	//drawing of after_image
	selectWindow("image_after");
	z1 = parseInt(zz2[i]-radius);
	if (z1<1){	z1 = 1;	}
	z2 = parseInt(zz2[i]+radius);
	if (z2>iz){	z2 = iz;	}
	
	y1 = parseInt(yy2[i]-radius);
	if (y1<0){	y1 = 0;	}
	y2 = parseInt(yy2[i]+radius);
	if (y2>=iy){	y2 = iy-1;	}	
	
	x1 = parseInt(xx2[i]-radius);
	if (x1<0){	x1 = 0;	}
	x2 = parseInt(xx2[i]+radius);
	if (x2>=ix){	x2 = ix-1;	}	

	for (z=z1; z<=z2; z++){
		setSlice(z);
		for (y=y1; y<y2; y++){
			for (x=x1; x<x2; x++){
				dist2 = pow(xx2[i]-x, 2) + pow(yy2[i]-y, 2) + pow(zz2[i]-z, 2);
				dist = sqrt(dist2);
				if (dist<=radius){
					if (modeOFintensity3==0){
						setPixel(x, y, (id1[i]+foot));
					}
					else{
						setPixel(x, y, (id2[i]+foot));
					}
				}
			}
		}
	}
	wait(100);
}

function scale_xyz_revision ()
{
	//If the unit in the input text file is um, the unit is converted to pixel.
	for (i=0; i<p_num; i++){
		xx1[i]=xx1[i]/um_pixel;
		yy1[i]=yy1[i]/um_pixel;
		zz1[i]=zz1[i]/um_pixel;
		xx2[i]=xx2[i]/um_pixel;
		yy2[i]=yy2[i]/um_pixel;
		zz2[i]=zz2[i]/um_pixel;		
	}
	max_x=xx1[0];
	max_y=yy1[0];
	max_z=zz1[0];
	min_x=xx1[0];
	min_y=yy1[0];
	min_z=zz1[0];
	for (i=0; i<p_num; i++){
		if (xx1[i]>max_x){	max_x=xx1[i];	}
		if (yy1[i]>max_y){	max_y=yy1[i];	}
		if (zz1[i]>max_z){	max_z=zz1[i];	}
		if (xx1[i]<min_x){	min_x=xx1[i];	}
		if (yy1[i]<min_y){	min_y=yy1[i];	}
		if (zz1[i]<min_z){	min_z=zz1[i];	}
		
		if (xx2[i]>max_x){	max_x=xx2[i];	}
		if (yy2[i]>max_y){	max_y=yy2[i];	}
		if (zz2[i]>max_z){	max_z=zz2[i];	}
		if (xx2[i]<min_x){	min_x=xx2[i];	}
		if (yy2[i]<min_y){	min_y=yy2[i];	}
		if (zz2[i]<min_z){	min_z=zz2[i];	}
	}

	//For parameter setting mistakes detection
	if (min_x<0.0 || min_y<0.0 || min_z<0.0){
		print("There are particles with <0.0 xyz-coordinates, which will be out of output images.");
		print("		Probably, parameter setting of xyz-scales or the original image size may be wrong.");
		print(min_x, min_y, min_z);
	}
	if (max_x > ix || max_y > iy){
		print("There are particles with maxima xy-coordinates, which will be out of output images.");
		print("		Probably, parameter setting of xy-scales or the original image size may be wrong.");
		print(max_x, max_y);
	}

	//image size
	iz = max_z + parseInt(radius*2);
}

function xyz_coordinate_read ()
{
	max_id = 0;
	for (i=0; i<p_num; i++){
		id1[i]=parseInt(str2[i*9+head]);
		id2[i]=parseInt(str2[i*9+1+head]);
		xx1[i]=parseFloat(str2[i*9+3+head]);
		yy1[i]=parseFloat(str2[i*9+4+head]);
		zz1[i]=parseFloat(str2[i*9+5+head]);
		xx2[i]=parseFloat(str2[i*9+6+head]);
		yy2[i]=parseFloat(str2[i*9+7+head]);
		zz2[i]=parseFloat(str2[i*9+8+head]);
		if (id1[i]>max_id){
			max_id = id1[i];
		}
		if (id2[i]>max_id){
			max_id = id2[i];
		}		
	}
}
