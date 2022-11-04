//************* 3D image rotation **********************
//
// Written by Hiroshi Koyama in National Institute for Basic Biology, Japan, 2022/07/22
// hkoyama@nibb.ac.jp

// <Description>
// An input image is 3-dimensionally rotated accoring to three-angles for rotation.
// Read an image file to be rotated: the image is image-after used for Macro_3D_particle_registration_220721_6_v2.ijm.
//		The file name should be "input_after_image.tif".
// Read a text file containing xyz-coordinates of multiple particles, used for Macro_3D_particle_registration_220721_6_v2.ijm.
//		The file name is "input_xyz_registration.csv".
// Manually provide 3D rotation angles obtained from Macro_3D_particle_registration_220721_6_v2.ijm.
//		The angles are written in the last line of "out_rotation_results.txt" obtained from Macro_3D_particle_registation_220721_6_v2.ijm.
// The input_after_image.tif is 3D rotated according to the rotation angles, resulting in image_after_rotation.tif.
// The image_after_rotation.tif file has the same scale as the image-before used for Macro_3D_particle_registration_220721_6_v2.ijm.
// Therefore, we can merge these two images for comparison.
// We can also merge image_after_rotation.tif and images obtained from Macro_particle_drawing_220721_02.ijm.
 
// The procedures in this macro are a little bit complicated.
//		1. First, input_after_image.tif is modified to have the same scale as that of input_xyz_registration.csv.
//		2. Then, the image is rotated.
//		3. Finally, the scale is modified to have the same scale as that of image-before.
// Due to these multiple scale modifications, users must provide precise parameters related to the scales.
//		e.g. x_scale_before, x_pixel_before,,,
// Incorrect parameters cause different scales between image_after_rotation.tif and image-before or cause distorted images.
//	 Therefore, the merge of these two images becomes unsuccessful.
// To avoid the complicated parameter setting and the troubles about the scales,
// I recommend that you set the same scale among image-before, image-after, and input_xyz_registration.csv. 

// <Tips for visualization after obtaining the output images>
// To get merged images, apply ImageJ>Image>Color>Merge Channels...
// For visualization, maximum intensity projection (ImageJ>Image>Stacks>Z Projects...) and 3D viewer (Fiji>Plugins>3D Viewer) may be powerful.
// Good choice of lookup tables improves visualization: ImageJ>Image>Lookup Tables>Green, Magenta, Fire, Spectrum, etc.
 
// <Requirement>
// Fiji, version 1.53 or maybe later.
// The name of the image to be rotated is input_after_image.tif
// The name of the text file is input_xyz_registration.csv with a given format.

// <Running time>
// MacBook air (CPU: Core i7, Dual core, 1.7GHz; Memory: 8 GB 1600 MHz DDR3)
// In the case that range_ave = 1 and image sizes (image_size_x = 262, image_size_y = 262, and image_size_z = 165), the running time was ~15min.
// In the case that range_ave = 0 and image sizes (image_size_x = 262, image_size_y = 262, and image_size_z = 165), the running time was ~2min.
// Theoretically, if the image sizes are 2-fold along x, y, and z, the running time may become ~2^3 = 8-fold.
// Theoretically, if range_ave = 2, the running time may become ~(5x5x5)/(3x3x3) = 4~5-fold compared with range = 1.  

// <Reference>
// Koyama et al. "An ImageJ-based tool for three-dimensional registration between different types of microscopic images"
// Cite the reference above when you use this macro for any publication.
//
//*****************************************************************************************************************

setBatchMode(true);

//********** Parameters to be defined by users ********************************************************************

print("Please set several values in 'Parameters to be defined by users' on the macro code.\n");

//** Scaling of xyz-coordinate: if the values = 1.0, the xyz-coordinates in the input text file (input_xyz_registration.csv) are not modified.
//These values should be the same as those used in Macro_3D_particle_registration_220721_06_v2.ijm.
//for before-image
x_scale_before = 1.0;
y_scale_before = 1.0;
z_scale_before = 1.0;

//for after-image
x_scale_after = 1.0;
y_scale_after = 1.0;	
z_scale_after = 1.0;	

//** For transformation of xyz-coorinates in non-pixel units (e.g. micron) to pixel unit in the input text file.
//A. If the unit of the input text file (input_xyz_registration.csv) is pixel, the following values should be 1.0.
//B. If the unit is not pixel (e.g. micron), the following values should be set other than 1.0 as follows. 
//		The values are obtained from Image>Properties>"Pixel width", "Pixel height", and "Voxel depth".
//		These three parameters correspond to x, y, and z, respectively.
//for before-image
x_pixel_before = 0.4143204;
y_pixel_before = 0.4143204;
z_pixel_before = 0.575;

//for after-image
x_pixel_after = 0.4143204;
y_pixel_after = 0.4143204;
z_pixel_after = 0.625;

//Image sizes of output image: large pixel numbers (e.g. 500~) may cause time consuming.
image_size_x = 262;		//The pixel number of output images along x axis, which is the same as the original before-image.
image_size_y = 262;		//The pixel number of output images along y axis, which is the same as the original before-image.
image_size_z = 165;		//The voxel number of output images along z axis, which is the same as the original before-image.

//** Rotation angles for after-image, which should be obtained from
//output_rotation_results.txt after running the Macro_3D_particle_registration_220721_06_v2.ijm macro.
//The values are written in the last line of the text file, "best rotational angles...".
angle1=5.9529;
angle2=6.2313;
angle3=0.2541;

//** Nomalization of intensity gradient along z-slices
//If your images exhibit severe decay of intensity in deep z-slices, the above normalization may be effective.
z_normalize = 1;	//1 for Yes, 0 for No.

//** Averaging of 3D image drawing
//If range_ave = 0, no averaging. If range_ave = 1, 3x3x3 averaging. If range_ave = 2, 5x5x5 averaging.
// This averaging improves image quality, but time consuming.
range_ave = 1;	

//*******************************************************************************************************************


//***** input text file name and input after-image name *****
dir_input = getDirectory("Choose a Directory ");
input_file = dir_input+"input_xyz_registration.csv";
input_image = dir_input+"input_after_image.tif";

//***** global variables ***********************************************
var c_num=1000;
var rfx, rfy, rfz;
var buf_x, buf_y, buf_z;
var x_scaling, y_scaling, z_scaling;
var x_scaling2, y_scaling2, z_scaling2;
var x_pix, y_pix, z_pix;
var x_pix2, y_pix2, z_pix2;
var x_bf, y_bf, z_bf, z_bf2, x_af, y_af, z_af;
var psi;
var centroid_x1, centroid_y1, centroid_z1, centroid_x2, centroid_y2, centroid_z2;
var range;

var cell_order, cell_num_t;
var cell_share1, cell_pair;
var max_cell_number, cell_number;
var buf_t, buf_t2, buf_id;
var end_t;
var count0, count1, count2;
var landmark_num;

rfx = newArray(2*c_num);
rfy = newArray(2*c_num);
rfz = newArray(2*c_num);
psi = newArray(3);
cell_order = newArray(2*c_num);
cell_num_t = newArray(c_num);
cell_share1 = newArray(2*c_num);
cell_pair = newArray(c_num);

//*************************************************************************

//***** basic parameter setting *****
x_scaling = x_scale_before;
y_scaling = y_scale_before;
z_scaling = z_scale_before;
x_scaling2 = x_scale_after;
y_scaling2 = y_scale_after;
z_scaling2 = z_scale_after;
x_pix = x_pixel_before;
y_pix = y_pixel_before;
z_pix = z_pixel_before;
x_pix2 = x_pixel_after;
y_pix2 = y_pixel_after;
z_pix2 = z_pixel_after;
image_z = parseInt(image_size_z*z_pixel_before/x_pixel_before); //To satisfy the same value between um/pixel along x and um/voxel along z.
x_bf = image_size_x;
y_bf = image_size_y;
z_bf = image_size_z;
z_bf2 = image_z;
psi[0]=angle1;
psi[1]=angle2;
psi[2]=angle3;
range = range_ave;

//***** print parameter setting *****
print("input text file = ", input_file);
print("input image name = ", input_image);
print("xyz_scale before = ", x_scale_before, y_scale_before, z_scale_before);
print("xyz_scale after = ", x_scale_after, y_scale_after, z_scale_after);
print("um/pixel or um/voxel before = ", x_pixel_before, y_pixel_before, z_pixel_before);
print("um/pixel or um/voxel after = ", x_pixel_after, y_pixel_after, z_pixel_after);
print("three angles =", angle1, angle2, angle3);
print("original before-image sizes = ", image_size_x, image_size_y, image_size_z);
print("Averaging for image construction: ", range);
print("Image normalization along z: ", z_normalize);

//***** input text file reading *****
input_file_reading (input_file);

//***** xyz scale revision *****/
//Unit is modified to the same unit as that in Macro_3D_particle_registration_220721_06_v2.ijm.
xyz_scaling_unit ();

//***** extraction of landmark particles (t=0) and (t=1) *****/
landmark_particle_extraction ();

//***** centroid calculation of landmarks *****/
centroid_of_landmarks ();

//***** read after-image *****/
open(input_image);
selectWindow("input_after_image.tif");
x_af = getWidth();
y_af = getHeight();
z_af = nSlices;
if (z_normalize == 1){	//Intensity normalization along z-slices
	run("Enhance Contrast...", "saturated=0.3 normalize process_all");
}
print("Total z_slices = ", z_af);

//***** generation of rotated image for after-image *****/
//The numbers of xy pixels of the output image are set to be the same as the before-image.
//The number of z-slices will be set later in the next section "z-scaling so that...". 
newImage("img_af.tif", "32-bit black", image_size_x, image_size_y, image_z);
newImage("img_af_count.tif", "16-bit black", image_size_x, image_size_y, image_z);
image_rotation_3D ();	

//***** z-scaling so that the output image has the same number of z-slices as the before-image *****
selectWindow("Result of img_af.tif");
run("Scale...", "x=1.0 y=1.0 z=- width="+image_size_x+" height="+image_size_y+" depth="+image_size_z+" interpolation=Bilinear average process create");

//***** save images *****
saveAs("Tiff", dir_input+"/image_after_rotation.tif");
close();
selectWindow("Result of img_af.tif");
//saveAs("Tiff", dir_input+"/Result_img_af.tif");
close();

print("finished");
print("You may manually convert the output image 'image_after_rotation.tif' to 8-bit under appropriate brightness/contrast.");

//*** save log file ***
selectWindow("Log");
saveAs("Text", dir_input+"/Log_image_rotation.txt");

//******************* end of macro ********************************************


//******************* functions ***********************************************

function image_rotation_3D ()
{
	//* rotation angle */
	c1=cos(psi[0]);
    c2=cos(psi[1]);
    c3=cos(psi[2]);
	s1=sin(psi[0]);
    s2=sin(psi[1]);
    s3=sin(psi[2]);
	
	xx=0.1;
	yy=0.1;
	zz=0.1;
	xxx=0.1;
	yyy=0.1;
	zzz=0.1;
	val1=0.1;
	val2=0.1;
	val3=1;
	
	for (z=1; z<=z_af; z++){
		print("	z=", z);
		for (y=0; y<y_af; y++){
			for (x=0; x<x_af; x++){
				selectWindow("input_after_image.tif");
				setSlice(z);
				//Scale units are equivalent to those of the text file after applying xyz_scaling_unit ().
				xx = x * x_pixel_after * x_scaling2 - centroid_x2;
				yy = y * y_pixel_after * y_scaling2 - centroid_y2;
				zz = z * z_pixel_after * z_scaling2 - centroid_z2;
				val1 = parseFloat(getPixel(x,y));

				//* rotation: using rotation matrix with yaw-pitch-roll (not Euler angles) */
				xxx=(xx*c1*c2)+(yy*(c1*s2*s3-s1*c3))+(zz*(c1*s2*c3+s1*s3));
				yyy=(xx*s1*c2)+(yy*(s1*s2*s3+c1*c3))+(zz*(s1*s2*c3-c1*s3));
				zzz=(xx*(-s2))+(yy*c2*s3)+(zz*c2*c3);

				//* xyz-coordinates are transformed according to the scale and the centroid of the before-image.
				xxx = (xxx+centroid_x1)/(x_scaling * x_pixel_before);
				yyy = (yyy+centroid_y1)/(y_scaling * y_pixel_before);
				zzz = (zzz+centroid_z1)/(z_scaling * z_pixel_before);

				//* For initial image construction, the image is set to have the same value between um/pixel along x and um/voxel along z.
				// This process is meaningful for good visualization if the z-resolution in the before-image is extremely lower than xy-resolutions.
				// In other words, if the z-resolution is higher, this process is not important.
				zzz = zzz * z_pixel_before / x_pixel_before;

				//* dipicting a rotated image
				x4 = parseInt(xxx);
				y4 = parseInt(yyy);
				z4 = parseInt(zzz);
				
				//3D drawing: if range = 1, 3x3x3 averaging will be done.
				for (z5 = z4-range; z5<=z4+range; z5++){
					for (y5 = y4-range; y5<=y4+range; y5++){
						for (x5 = x4-range; x5<=x4+range; x5++){
							if (z5>0 && z5<z_bf && x5>=0 && x5<x_bf && y5>=0 && y5<y_bf){	//before-image's size
								selectWindow("img_af.tif");		//summation image
								setSlice(z5);
								val2 = getPixel(x5, y5);
								val2 = val2 + val1;
								setPixel(x5, y5, val2);

								selectWindow("img_af_count.tif");	//count image
								setSlice(z5);
								val3 = getPixel(x5, y5);
								val3++;
								setPixel(x5, y5, val3);
							}	
						}
					}
				}
				
			}
		}
	}

	//Divide the summation image by the count image
	imageCalculator("Divide create stack", "img_af.tif","img_af_count.tif");
	selectWindow("Result of img_af.tif");	//intensity image
	//run("8-bit");							//Conversion to 8-bit may generate an inappropriate brightness/contrast image.

	selectWindow("img_af.tif");
	//saveAs("Tiff", dir_input+"/img_af.tif");
	close();
	selectWindow("img_af_count.tif");
	//saveAs("Tiff", dir_input+"/img_af_count.tif");
	close();
	
}

function input_file_reading (input)
{
	str1 = File.openAsString(input);
	str2 = split(str1, "\t\n\r,");
	s_num = str2.length;
	p_num = s_num/5;
	print("Particle number = ", p_num);
	print("Other values in the input file: ", s_num, str2[0], str2[s_num-1]);

	count0=0;
       	count1=0;
        for (i=0; i<p_num; i++)
        {
		buf_t = parseInt(str2[i*5+3]);

		if (i>0 && buf_t>end_t){
			count0++;
		       	count1=0;
		   	if (buf_t!=end_t+1){
				print("Error. Time is skipped!");
				exit(0);
			}
		}

                buf_id = parseInt(str2[i*5+4]) -1;
                end_t = buf_t;
                rfx[count0*c_num+count1] = parseFloat(str2[i*5]);
	       		rfy[count0*c_num+count1] = parseFloat(str2[i*5+1]);
	       		rfz[count0*c_num+count1] = parseFloat(str2[i*5+2]);;
                cell_order[count0*c_num+count1]=buf_id;
                cell_num_t[count0]=count1;
                count1++;
        }
        print("cell_number = ", cell_num_t[0], "(before)", cell_num_t[1], "(after)");
	max_cell_number=cell_num_t[0];
	if (cell_num_t[1]>max_cell_number){
		max_cell_number=cell_num_t[1];
	}
	if (max_cell_number>=c_num){	
		print("Error in cell number.");
		exit(0);
	}
}

function xyz_scaling_unit ()
{
	//xyz coordinates are modified according to x_scaling, y_scaling, and z_scaling
	//Also, the units of the xyz coordinates are modified to be pixel but not micron etc.

	t=0;
    for (i=0; i<=cell_num_t[t]; i++)
    {
		rfx[t*c_num+i]=rfx[t*c_num+i]*x_scaling;
		rfy[t*c_num+i]=rfy[t*c_num+i]*y_scaling;
		rfz[t*c_num+i]=rfz[t*c_num+i]*z_scaling;
	}

	t=1;
    for (i=0; i<=cell_num_t[t]; i++)
    {		
		rfx[t*c_num+i]=rfx[t*c_num+i]*x_scaling2;
	    rfy[t*c_num+i]=rfy[t*c_num+i]*y_scaling2;
	    rfz[t*c_num+i]=rfz[t*c_num+i]*z_scaling2;
	}
}

function landmark_particle_extraction ()
{
	for (t=0; t<=1; t++){
		for (i=0; i<=cell_num_t[t]; i++){
		       	cell_share1[t*c_num+i]=0;
		}
	}
	t=0;
	count2=-1;
	for (i=0; i<=cell_num_t[t]; i++){
		buf_id=cell_order[t*c_num+i];
		for (j=0; j<=cell_num_t[t+1]; j++){
			if (buf_id==cell_order[(t+1)*c_num+j]){
				cell_share1[t*c_num+i]=1;
				cell_share1[(t+1)*c_num+j]=1;
				cell_pair[i]=j;
				count2++;
			}
		}
	}
	landmark_num = count2;
	print("Number of landmark particles = ", landmark_num);
}

function centroid_of_landmarks ()
{
	t=1;

	//centroid of landmark particles in before-image
	count0=0;
	centroid_x1=0.0;
    centroid_y1=0.0;
    centroid_z1=0.0;
	for (i=0; i<=cell_num_t[t-1]; i++){
		if (cell_share1[(t-1)*c_num+i]==1){	//choose landmark particles
			count0++;
			centroid_x1=centroid_x1+rfx[(t-1)*c_num+i];
			centroid_y1=centroid_y1+rfy[(t-1)*c_num+i];
			centroid_z1=centroid_z1+rfz[(t-1)*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x1=centroid_x1/parseFloat(count0);
		centroid_y1=centroid_y1/parseFloat(count0);
		centroid_z1=centroid_z1/parseFloat(count0);
	}

	//centroid of landmark particles in after-image
	count0=0;
	centroid_x2=0.0;
    centroid_y2=0.0;
    centroid_z2=0.0;
	for (i=0; i<=cell_num_t[t]; i++){
		if (cell_share1[t*c_num+i]==1){	//choose landmark particles
			count0++;
			centroid_x2=centroid_x2+rfx[t*c_num+i];
			centroid_y2=centroid_y2+rfy[t*c_num+i];
			centroid_z2=centroid_z2+rfz[t*c_num+i];
		}
	}
	if (count0!=0){
		centroid_x2=centroid_x2/parseFloat(count0);
		centroid_y2=centroid_y2/parseFloat(count0);
		centroid_z2=centroid_z2/parseFloat(count0);
	}
}


